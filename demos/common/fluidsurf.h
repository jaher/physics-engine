// Screen-space fluid surface renderer. Turns a cloud of SPH particles into a
// smooth liquid surface instead of visible spheres:
//   1. render sphere imposters to an eye-space depth (+ temperature) buffer,
//   2. bilaterally blur that depth into a continuous surface,
//   3. reconstruct per-pixel normals from the smoothed depth,
//   4. shade it (lava emissive ramp, or a glossy Fresnel liquid) and composite
//      over a simply-lit opaque scene, then ACES-tonemap.
// Self-contained: owns its FBOs and shaders, reuses gfx::Mesh/program/writePNG.
#pragma once
#include "gfx.h"
#include <vector>

namespace gfx {

struct FluidSurf {
    int W = 0, H = 0;
    GLuint sceneFBO = 0, sceneColor = 0, sceneLin = 0, sceneDepthRB = 0;
    GLuint fFBO = 0, fTex = 0, fDepthRB = 0;
    GLuint bFBO[2] = {0, 0}, bTex[2] = {0, 0};
    GLuint oFBO = 0, oTex = 0;
    GLuint pScene = 0, pImp = 0, pBlur = 0, pComp = 0, quadVAO = 0;
    glm::mat4 view{1}, proj{1}; glm::vec3 lightDir{0, -1, 0}, camEye{0};

    static GLuint tex(int w, int h, GLint internal, GLenum fmt, GLenum type) {
        GLuint t; glGenTextures(1, &t); glBindTexture(GL_TEXTURE_2D, t);
        glTexImage2D(GL_TEXTURE_2D, 0, internal, w, h, 0, fmt, type, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        return t;
    }

    void init(int w, int h) {
        W = w; H = h;
        glGenVertexArrays(1, &quadVAO);
        // scene: HDR colour + linear eye-depth + depth renderbuffer
        sceneColor = tex(W, H, GL_RGBA16F, GL_RGBA, GL_FLOAT);
        sceneLin = tex(W, H, GL_R32F, GL_RED, GL_FLOAT);
        glGenRenderbuffers(1, &sceneDepthRB); glBindRenderbuffer(GL_RENDERBUFFER, sceneDepthRB);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, W, H);
        glGenFramebuffers(1, &sceneFBO); glBindFramebuffer(GL_FRAMEBUFFER, sceneFBO);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, sceneColor, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, sceneLin, 0);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, sceneDepthRB);
        // fluid: (temp, linDepth, -, coverage) + depth renderbuffer
        fTex = tex(W, H, GL_RGBA16F, GL_RGBA, GL_FLOAT);
        glGenRenderbuffers(1, &fDepthRB); glBindRenderbuffer(GL_RENDERBUFFER, fDepthRB);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, W, H);
        glGenFramebuffers(1, &fFBO); glBindFramebuffer(GL_FRAMEBUFFER, fFBO);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fTex, 0);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, fDepthRB);
        for (int i = 0; i < 2; i++) {
            bTex[i] = tex(W, H, GL_RGBA16F, GL_RGBA, GL_FLOAT);
            glGenFramebuffers(1, &bFBO[i]); glBindFramebuffer(GL_FRAMEBUFFER, bFBO[i]);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, bTex[i], 0);
        }
        oTex = tex(W, H, GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE);
        glGenFramebuffers(1, &oFBO); glBindFramebuffer(GL_FRAMEBUFFER, oFBO);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, oTex, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        compileShaders();
    }

    void setCamera(const glm::mat4& v, const glm::mat4& p, const glm::vec3& eye) { view = v; proj = p; camEye = eye; }
    void setLight(const glm::vec3& d) { lightDir = glm::normalize(d); }

    // ---- opaque scene ----
    void beginScene(glm::vec3 sky) {
        glBindFramebuffer(GL_FRAMEBUFFER, sceneFBO); glViewport(0, 0, W, H);
        GLenum bufs[2] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1}; glDrawBuffers(2, bufs);
        glEnable(GL_DEPTH_TEST); glDepthMask(GL_TRUE);
        glClearColor(sky.r, sky.g, sky.b, 1.0f);
        float farLin[4] = {1e9f, 1e9f, 1e9f, 1e9f}; glClearBufferfv(GL_COLOR, 1, farLin);   // linear-depth = far
        glClear(GL_DEPTH_BUFFER_BIT);
        float skyc[4] = {sky.r, sky.g, sky.b, 1.0f}; glClearBufferfv(GL_COLOR, 0, skyc);
    }
    void drawMesh(const Mesh& m, const glm::mat4& model, glm::vec3 albedo, float rough, int isGround = 0) {
        glUseProgram(pScene);
        setM4(pScene, "uModel", model); setM4(pScene, "uView", view); setM4(pScene, "uProj", proj);
        setV3(pScene, "uAlbedo", albedo); setV3(pScene, "uLightDir", lightDir); setV3(pScene, "uCamPos", camEye);
        setF(pScene, "uRough", rough); setI(pScene, "uGround", isGround);
        m.draw();
    }

    // ---- fluid surface ----  data: [x y z temp] per particle, stride floats
    void surface(const std::vector<float>& data, int stride, float radius) {
        static GLuint vao = 0, vbo = 0; static std::vector<float> quads;
        if (!vao) { glGenVertexArrays(1, &vao); glGenBuffers(1, &vbo); }
        quads.clear();
        const float uv[6][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, -1}, {1, 1}, {-1, 1}};
        int n = (int)data.size() / stride;
        for (int i = 0; i < n; i++) { const float* p = &data[i * stride];
            for (auto& c : uv) quads.insert(quads.end(), {p[0], p[1], p[2], c[0], c[1], p[3]}); }
        glBindVertexArray(vao); glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, quads.size() * sizeof(float), quads.data(), GL_DYNAMIC_DRAW);
        glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(1); glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(2); glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(5 * sizeof(float)));

        glBindFramebuffer(GL_FRAMEBUFFER, fFBO); glViewport(0, 0, W, H);
        GLenum b0 = GL_COLOR_ATTACHMENT0; glDrawBuffers(1, &b0);
        glEnable(GL_DEPTH_TEST); glDepthMask(GL_TRUE);
        float clr[4] = {0, 1e9f, 0, 0}; glClearBufferfv(GL_COLOR, 0, clr); glClear(GL_DEPTH_BUFFER_BIT);
        glUseProgram(pImp);
        setM4(pImp, "uView", view); setM4(pImp, "uProj", proj); setF(pImp, "uRadius", radius);
        setV3(pImp, "uInvRes", glm::vec3(1.0f / W, 1.0f / H, 0));
        glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, sceneLin); setI(pImp, "uSceneLin", 0);
        glBindVertexArray(vao); glDrawArrays(GL_TRIANGLES, 0, n * 6);
    }

    void blur(int iters, float depthSigma) {
        glDisable(GL_DEPTH_TEST); glBindVertexArray(quadVAO); glUseProgram(pBlur);
        setF(pBlur, "uRangeScale", 1.0f / (2.0f * depthSigma * depthSigma));
        GLuint src = fTex;
        for (int i = 0; i < iters; i++) {
            for (int dir = 0; dir < 2; dir++) {
                int dst = (i * 2 + dir) % 2;
                glBindFramebuffer(GL_FRAMEBUFFER, bFBO[dst]); glViewport(0, 0, W, H);
                setV3(pBlur, "uDir", dir == 0 ? glm::vec3(1.0f / W, 0, 0) : glm::vec3(0, 1.0f / H, 0));
                glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, src); setI(pBlur, "uTex", 0);
                glDrawArrays(GL_TRIANGLES, 0, 3);
                src = bTex[dst];
            }
        }
        lastBlur = src;
    }
    GLuint lastBlur = 0;

    // mode 0 = lava (temperature emissive); 1 = glossy liquid whose key channel
    // (the 4th particle float) lerps colour from baseA (0) to baseB (1)
    void composite(int mode, glm::vec3 baseA, glm::vec3 baseB = glm::vec3(0)) {
        glBindFramebuffer(GL_FRAMEBUFFER, oFBO); glViewport(0, 0, W, H);
        glDisable(GL_DEPTH_TEST); glUseProgram(pComp);
        glm::vec3 Lv = glm::normalize(glm::mat3(view) * (-lightDir));
        float tanHalf = 1.0f / proj[1][1], aspect = proj[1][1] / proj[0][0];
        setV3(pComp, "uLightV", Lv); setF(pComp, "uTanHalf", tanHalf); setF(pComp, "uAspect", aspect);
        setV3(pComp, "uInvRes", glm::vec3(1.0f / W, 1.0f / H, 0)); setI(pComp, "uMode", mode);
        setV3(pComp, "uBase", baseA); setV3(pComp, "uBase2", baseB);
        glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, lastBlur ? lastBlur : fTex); setI(pComp, "uFluid", 0);
        glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, sceneColor); setI(pComp, "uScene", 1);
        glBindVertexArray(quadVAO); glDrawArrays(GL_TRIANGLES, 0, 3);
        glEnable(GL_DEPTH_TEST);
    }

    bool screenshot(const char* path) {
        std::vector<unsigned char> px((size_t)W * H * 3);
        glBindFramebuffer(GL_READ_FRAMEBUFFER, oFBO); glReadBuffer(GL_COLOR_ATTACHMENT0);
        glReadPixels(0, 0, W, H, GL_RGB, GL_UNSIGNED_BYTE, px.data());
        return writePNG(path, px.data(), W, H, true);
    }
    void present() {
        glBindFramebuffer(GL_READ_FRAMEBUFFER, oFBO); glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glBlitFramebuffer(0, 0, W, H, 0, 0, W, H, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    }

    void compileShaders() {
        const char* VS_SCN = R"(#version 330 core
layout(location=0) in vec3 aPos; layout(location=1) in vec3 aNormal;
uniform mat4 uModel,uView,uProj; out vec3 vW; out vec3 vN; out vec3 vE;
void main(){ vec4 w=uModel*vec4(aPos,1.0); vW=w.xyz; vN=mat3(uModel)*aNormal;
  vE=(uView*w).xyz; gl_Position=uProj*uView*w; })";
        const char* FS_SCN = R"(#version 330 core
in vec3 vW; in vec3 vN; in vec3 vE; layout(location=0) out vec4 oCol; layout(location=1) out float oLin;
uniform vec3 uAlbedo,uLightDir,uCamPos; uniform float uRough; uniform int uGround;
void main(){ vec3 N=normalize(vN); vec3 L=normalize(-uLightDir); vec3 V=normalize(uCamPos-vW);
  vec3 alb=uAlbedo;
  if(uGround==1){ float c=mod(floor(vW.x*1.4)+floor(vW.z*1.4),2.0); alb*=(0.8+0.2*c); }
  float dif=max(dot(N,L),0.0); vec3 H=normalize(L+V);
  float sp=pow(max(dot(N,H),0.0), mix(8.0,64.0,1.0-uRough))*(1.0-uRough)*0.4;
  vec3 col=alb*(0.25+0.85*dif)+vec3(sp);
  oCol=vec4(col,1.0); oLin=-vE.z; })";

        const char* VS_IMP = R"(#version 330 core
layout(location=0) in vec3 aCenter; layout(location=1) in vec2 aUV; layout(location=2) in float aTemp;
uniform mat4 uView,uProj; uniform float uRadius; out vec2 vUV; out float vTemp; out vec3 vC;
void main(){ vUV=aUV; vTemp=aTemp; vC=(uView*vec4(aCenter,1.0)).xyz;
  vec3 p=vC+vec3(aUV*uRadius,0.0); gl_Position=uProj*vec4(p,1.0); })";
        const char* FS_IMP = R"(#version 330 core
in vec2 vUV; in float vTemp; in vec3 vC; out vec4 frag;
uniform mat4 uProj; uniform float uRadius; uniform vec3 uInvRes; uniform sampler2D uSceneLin;
void main(){ float r2=dot(vUV,vUV); if(r2>1.0) discard;
  vec3 N=vec3(vUV,sqrt(1.0-r2)); vec3 eye=vC+N*uRadius; float lin=-eye.z;
  float sceneLin=texture(uSceneLin, gl_FragCoord.xy*uInvRes.xy).r;
  if(lin>sceneLin) discard;                         // behind opaque scene
  vec4 clip=uProj*vec4(eye,1.0); gl_FragDepth=(clip.z/clip.w)*0.5+0.5;
  frag=vec4(vTemp, lin, 0.0, 1.0); })";

        const char* FS_BLUR = R"(#version 330 core
out vec4 frag; uniform sampler2D uTex; uniform vec3 uDir; uniform float uRangeScale;
void main(){ vec2 uv=gl_FragCoord.xy*vec2(1.0)/vec2(textureSize(uTex,0));
  vec4 c=texture(uTex,uv); if(c.a<0.001){ frag=vec4(0.0,1e9,0.0,0.0); return; }
  float dc=c.g; float wsum=0.0, ts=0.0, ds=0.0, cs=0.0, gwsum=0.0;
  for(int i=-8;i<=8;i++){ vec2 o=uv+uDir.xy*float(i); vec4 s=texture(uTex,o);
    float gw=exp(-float(i*i)/32.0); gwsum+=gw;
    float rw=(s.a>0.001)? exp(-(s.g-dc)*(s.g-dc)*uRangeScale) : 0.0;
    float w=gw*rw; wsum+=w; ts+=w*s.r; ds+=w*s.g; cs+=gw*s.a; }
  frag=vec4(ts/max(wsum,1e-5), ds/max(wsum,1e-5), 0.0, cs/max(gwsum,1e-5)); })";

        const char* FS_COMP = R"(#version 330 core
out vec4 frag; uniform sampler2D uFluid, uScene;
uniform vec3 uLightV, uInvRes, uBase, uBase2; uniform float uTanHalf, uAspect; uniform int uMode;
vec3 aces(vec3 x){return clamp((x*(2.51*x+0.03))/(x*(2.43*x+0.59)+0.14),0.0,1.0);}
vec3 eyePos(vec2 uv,float lin){ return vec3((2.0*uv.x-1.0)*uAspect*uTanHalf*lin, (2.0*uv.y-1.0)*uTanHalf*lin, -lin); }
void main(){ vec2 uv=gl_FragCoord.xy*uInvRes.xy;
  vec3 scene=texture(uScene,uv).rgb;
  vec4 f=texture(uFluid,uv); float cov=f.a;
  if(cov<0.02){ frag=vec4(aces(scene),1.0); return; }
  float lin=f.g; vec3 P=eyePos(uv,lin);
  vec2 px=uInvRes.xy;
  float lr=texture(uFluid,uv+vec2(px.x,0)).g, ll=texture(uFluid,uv-vec2(px.x,0)).g;
  float lu=texture(uFluid,uv+vec2(0,px.y)).g, ld=texture(uFluid,uv-vec2(0,px.y)).g;
  vec3 ddx=(abs(lr-lin)<abs(lin-ll))? eyePos(uv+vec2(px.x,0),lr)-P : P-eyePos(uv-vec2(px.x,0),ll);
  vec3 ddy=(abs(lu-lin)<abs(lin-ld))? eyePos(uv+vec2(0,px.y),lu)-P : P-eyePos(uv-vec2(0,px.y),ld);
  vec3 N=normalize(cross(ddx,ddy)); if(N.z<0.0) N=-N;
  vec3 L=normalize(uLightV); vec3 V=vec3(0,0,1); vec3 Hh=normalize(L+V);
  float dif=max(dot(N,L),0.0); float spec=pow(max(dot(N,Hh),0.0),60.0);
  float fres=pow(1.0-max(dot(N,V),0.0),3.0);
  vec3 col;
  if(uMode==0){                                     // lava: crust + molten emissive
    float T=f.r; vec3 crust=vec3(0.05)+vec3(0.16,0.05,0.025)*T;
    float e=pow(max(T,0.0),2.3)*3.8;
    vec3 emis=vec3(1.0, pow(max(T,0.0),1.7)*0.9, pow(max(T,0.0),4.5)*0.7)*e;
    col=crust*(0.25+0.7*dif)+emis+vec3(spec)*0.4*T+emis*fres*0.6;
  } else {                                          // glossy liquid (key lerps colour)
    vec3 base=mix(uBase, uBase2, clamp(f.r,0.0,1.0));
    col=base*(0.28+0.75*dif)+vec3(spec)*0.9+base*fres*0.7+vec3(fres)*0.15;
  }
  float a=smoothstep(0.15,0.6,cov);
  vec3 outc=mix(scene,col,a);
  frag=vec4(aces(outc),1.0); })";

        pScene = program(VS_SCN, FS_SCN);
        pImp = program(VS_IMP, FS_IMP);
        pBlur = program(VS_FULL, FS_BLUR);
        pComp = program(VS_FULL, FS_COMP);
    }
};

} // namespace gfx
