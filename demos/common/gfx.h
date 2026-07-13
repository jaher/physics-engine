// A compact modern-OpenGL renderer for the demos: HDR + 4x MSAA offscreen target,
// directional-light shadow mapping with PCF, a Cook-Torrance PBR pass (with a
// cloth-sheen and a Kajiya-Kay hair variant), and an ACES tonemap resolve. Plus
// mesh generators, an orbit camera, and PNG screenshots. Uses GLFW/GLEW/GLM.
#pragma once
#define GLFW_INCLUDE_NONE
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "png.h"
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>

namespace gfx {

// ------------------------------------------------------------ shader utilities
inline GLuint compile(GLenum type, const char* src) {
    GLuint s = glCreateShader(type); glShaderSource(s, 1, &src, nullptr); glCompileShader(s);
    GLint ok; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) { char log[4096]; glGetShaderInfoLog(s, 4096, nullptr, log); std::fprintf(stderr, "shader error:\n%s\n", log); }
    return s;
}
inline GLuint program(const char* vs, const char* fs) {
    GLuint p = glCreateProgram(), v = compile(GL_VERTEX_SHADER, vs), f = compile(GL_FRAGMENT_SHADER, fs);
    glAttachShader(p, v); glAttachShader(p, f); glLinkProgram(p);
    GLint ok; glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) { char log[4096]; glGetProgramInfoLog(p, 4096, nullptr, log); std::fprintf(stderr, "link error:\n%s\n", log); }
    glDeleteShader(v); glDeleteShader(f); return p;
}
inline void setM4(GLuint p, const char* n, const glm::mat4& m) { glUniformMatrix4fv(glGetUniformLocation(p, n), 1, GL_FALSE, glm::value_ptr(m)); }
inline void setV3(GLuint p, const char* n, const glm::vec3& v) { glUniform3fv(glGetUniformLocation(p, n), 1, glm::value_ptr(v)); }
inline void setF(GLuint p, const char* n, float f) { glUniform1f(glGetUniformLocation(p, n), f); }
inline void setI(GLuint p, const char* n, int i) { glUniform1i(glGetUniformLocation(p, n), i); }

// ------------------------------------------------------------------- geometry
struct Mesh {
    GLuint vao = 0, vbo = 0, ebo = 0; GLsizei count = 0;
    // interleaved: pos(3) normal(3)  [+ uv(2) if hasUV]
    void upload(const std::vector<float>& verts, const std::vector<unsigned>& idx, bool hasUV = false) {
        int stride = (hasUV ? 8 : 6);
        if (!vao) glGenVertexArrays(1, &vao);
        if (!vbo) glGenBuffers(1, &vbo);
        if (!ebo) glGenBuffers(1, &ebo);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float), verts.data(), GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, idx.size() * sizeof(unsigned), idx.data(), GL_DYNAMIC_DRAW);
        glEnableVertexAttribArray(0); glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride * sizeof(float), (void*)0);
        glEnableVertexAttribArray(1); glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride * sizeof(float), (void*)(3 * sizeof(float)));
        if (hasUV) { glEnableVertexAttribArray(2); glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, stride * sizeof(float), (void*)(6 * sizeof(float))); }
        glBindVertexArray(0);
        count = (GLsizei)idx.size();
    }
    // update only the vertex buffer (for animated meshes with a fixed index set)
    void updateVerts(const std::vector<float>& verts) {
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferSubData(GL_ARRAY_BUFFER, 0, verts.size() * sizeof(float), verts.data());
    }
    void draw() const { glBindVertexArray(vao); glDrawElements(GL_TRIANGLES, count, GL_UNSIGNED_INT, 0); }
};

inline Mesh makeSphere(int seg = 48, int rings = 32) {
    std::vector<float> v; std::vector<unsigned> idx;
    for (int y = 0; y <= rings; y++) for (int x = 0; x <= seg; x++) {
        float u = (float)x / seg, t = (float)y / rings;
        float px = std::cos(u * 2 * M_PI) * std::sin(t * M_PI);
        float py = std::cos(t * M_PI);
        float pz = std::sin(u * 2 * M_PI) * std::sin(t * M_PI);
        v.insert(v.end(), {px, py, pz, px, py, pz});
    }
    int w = seg + 1;
    for (int y = 0; y < rings; y++) for (int x = 0; x < seg; x++) {
        unsigned a = y * w + x, b = a + 1, c = a + w, d = c + 1;
        idx.insert(idx.end(), {a, c, b, b, c, d});
    }
    Mesh m; m.upload(v, idx); return m;
}
inline Mesh makeBox() {
    const float n[6][3] = {{0,0,1},{0,0,-1},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0}};
    const float f[6][12] = {
        {-1,-1,1, 1,-1,1, 1,1,1, -1,1,1}, {1,-1,-1, -1,-1,-1, -1,1,-1, 1,1,-1},
        {1,-1,1, 1,-1,-1, 1,1,-1, 1,1,1}, {-1,-1,-1, -1,-1,1, -1,1,1, -1,1,-1},
        {-1,1,1, 1,1,1, 1,1,-1, -1,1,-1}, {-1,-1,-1, 1,-1,-1, 1,-1,1, -1,-1,1}};
    std::vector<float> v; std::vector<unsigned> idx;
    for (int s = 0; s < 6; s++) { unsigned base = (unsigned)(v.size() / 6);
        for (int k = 0; k < 4; k++) v.insert(v.end(), {f[s][k*3], f[s][k*3+1], f[s][k*3+2], n[s][0], n[s][1], n[s][2]});
        idx.insert(idx.end(), {base, base+1, base+2, base, base+2, base+3}); }
    Mesh m; m.upload(v, idx); return m;
}
inline Mesh makeTorus(float R = 1.0f, float r = 0.28f, int seg = 36, int ring = 18) {
    std::vector<float> v; std::vector<unsigned> idx;
    for (int i = 0; i <= seg; i++) for (int k = 0; k <= ring; k++) {
        float u = 2 * M_PI * i / seg, w = 2 * M_PI * k / ring;
        float cx = std::cos(u) * R, cy = std::sin(u) * R;
        float nx = std::cos(u) * std::cos(w), ny = std::sin(u) * std::cos(w), nz = std::sin(w);
        v.insert(v.end(), {cx + r * nx, cy + r * ny, r * nz, nx, ny, nz});
    }
    int W2 = ring + 1;
    for (int i = 0; i < seg; i++) for (int k = 0; k < ring; k++) {
        unsigned a = i * W2 + k, b = a + 1, c = a + W2, d = c + 1;
        idx.insert(idx.end(), {a, c, b, b, c, d});
    }
    Mesh m; m.upload(v, idx); return m;
}
inline Mesh makePlane(float size, int div = 1) {
    std::vector<float> v; std::vector<unsigned> idx;
    for (int y = 0; y <= div; y++) for (int x = 0; x <= div; x++) {
        float px = (float(x) / div - 0.5f) * 2 * size, pz = (float(y) / div - 0.5f) * 2 * size;
        v.insert(v.end(), {px, 0, pz, 0, 1, 0});
    }
    for (int y = 0; y < div; y++) for (int x = 0; x < div; x++) {
        unsigned a = y * (div + 1) + x, b = a + 1, c = a + div + 1, d = c + 1;
        idx.insert(idx.end(), {a, c, b, b, c, d});
    }
    Mesh m; m.upload(v, idx); return m;
}

// -------------------------------------------------------------------- camera
struct OrbitCamera {
    glm::vec3 target = glm::vec3(0, 1, 0);
    float yaw = 0.6f, pitch = 0.25f, dist = 8.0f;
    glm::vec3 eye() const {
        return target + glm::vec3(std::cos(pitch) * std::sin(yaw), std::sin(pitch), std::cos(pitch) * std::cos(yaw)) * dist;
    }
    glm::mat4 view() const { return glm::lookAt(eye(), target, glm::vec3(0, 1, 0)); }
};

// ------------------------------------------------------------------- shaders
static const char* VS_SCENE = R"(#version 330 core
layout(location=0) in vec3 aPos; layout(location=1) in vec3 aNormal; layout(location=2) in vec2 aUV;
uniform mat4 uModel,uView,uProj,uLightSpace;
out vec3 vWorld; out vec3 vNormal; out vec4 vLS; out vec2 vUV;
void main(){ vec4 wp=uModel*vec4(aPos,1.0); vWorld=wp.xyz; vNormal=mat3(uModel)*aNormal;
  vLS=uLightSpace*wp; vUV=aUV; gl_Position=uProj*uView*wp; })";

static const char* PBR_HEAD = R"(#version 330 core
in vec3 vWorld; in vec3 vNormal; in vec4 vLS; in vec2 vUV; out vec4 frag;
uniform vec3 uCamPos,uLightDir,uLightColor,uAlbedo; uniform float uRough,uMetal,uAO; uniform sampler2D uShadow;
const float PI=3.14159265;
float D_GGX(vec3 N,vec3 H,float a){float a2=a*a;float d=max(dot(N,H),0.0); d=(d*d*(a2-1.0)+1.0); return a2/(PI*d*d);}
float G1(float nx,float k){return nx/(nx*(1.0-k)+k);}
float G_S(vec3 N,vec3 V,vec3 L,float r){float k=(r+1.0);k=k*k/8.0;return G1(max(dot(N,V),0.0),k)*G1(max(dot(N,L),0.0),k);}
vec3 F_S(float ct,vec3 F0){return F0+(1.0-F0)*pow(clamp(1.0-ct,0.0,1.0),5.0);}
float shadow(vec4 ls,vec3 N,vec3 L){ vec3 p=ls.xyz/ls.w*0.5+0.5; if(p.z>1.0)return 1.0;
  float b=max(0.0022*(1.0-dot(N,L)),0.0006); float s=0.0; vec2 tx=1.0/vec2(textureSize(uShadow,0));
  for(int x=-2;x<=2;x++)for(int y=-2;y<=2;y++){float d=texture(uShadow,p.xy+vec2(x,y)*tx).r; s+=(p.z-b>d)?0.0:1.0;} return s/25.0; }
vec3 lighting(vec3 N,vec3 albedo,float rough,float metal,float ao){
  vec3 V=normalize(uCamPos-vWorld), L=normalize(-uLightDir), H=normalize(V+L);
  vec3 F0=mix(vec3(0.04),albedo,metal); float NdL=max(dot(N,L),0.0);
  float a=max(rough,0.04); float D=D_GGX(N,H,a*a); float G=G_S(N,V,L,rough); vec3 F=F_S(max(dot(H,V),0.0),F0);
  vec3 spec=D*G*F/max(4.0*max(dot(N,V),0.0)*NdL,0.001);
  vec3 kd=(1.0-F)*(1.0-metal); float sh=shadow(vLS,N,L);
  vec3 Lo=(kd*albedo/PI+spec)*uLightColor*NdL*sh;
  vec3 sky=vec3(0.5,0.62,0.82), grd=vec3(0.28,0.24,0.20);
  vec3 amb=mix(grd,sky,N.y*0.5+0.5)*albedo*ao;
  vec3 fres=F_S(max(dot(N,V),0.0),F0)*mix(grd,sky,0.7)*0.4*(1.0-rough);
  return amb+fres+Lo;
})";

static const char* FS_SOLID = R"(
void main(){ frag=vec4(lighting(normalize(vNormal),uAlbedo,uRough,uMetal,uAO),1.0); })";

// ground: procedural subtle checker + fade
static const char* FS_GROUND = R"(
void main(){ vec3 N=normalize(vNormal);
  float c=mod(floor(vWorld.x*0.5)+floor(vWorld.z*0.5),2.0);
  vec3 alb=mix(uAlbedo,uAlbedo*0.82,c);
  float r=length(vWorld.xz); alb*=smoothstep(60.0,18.0,r)*0.9+0.1;
  frag=vec4(lighting(N,alb,uRough,uMetal,uAO),1.0); })";

// cloth: two-sided, woven micro-detail, retro-reflective sheen + faint translucency
static const char* FS_CLOTH = R"(
void main(){ vec3 N=normalize(vNormal); if(!gl_FrontFacing) N=-N;
  vec3 V=normalize(uCamPos-vWorld);
  vec2 tu=vUV*110.0;                                        // woven thread pattern
  float weave=0.5+0.5*sin(tu.x*6.2831853)*sin(tu.y*6.2831853);
  float weave2=0.5+0.5*sin((tu.x+tu.y)*3.14159);
  vec3 alb=uAlbedo*(0.86+0.14*weave)*(0.95+0.05*weave2);
  float rough=clamp(uRough+0.10*(weave-0.5),0.3,1.0);
  vec3 base=lighting(N,alb,rough,uMetal,uAO);
  float sheen=pow(1.0-max(dot(N,V),0.0),3.0); base+=sheen*uLightColor*0.12*alb;   // cloth grazing sheen
  vec3 L=normalize(-uLightDir); float back=max(dot(-N,L),0.0); base+=back*alb*0.16*uLightColor; // subsurface hint
  frag=vec4(base,1.0); })";

static const char* VS_DEPTH = R"(#version 330 core
layout(location=0) in vec3 aPos; uniform mat4 uModel,uLightSpace;
void main(){ gl_Position=uLightSpace*uModel*vec4(aPos,1.0); })";
static const char* FS_DEPTH = R"(#version 330 core
void main(){})";

// hair: ribbon geometry with Kajiya-Kay anisotropic shading (two shifted lobes)
static const char* VS_HAIR = R"(#version 330 core
layout(location=0) in vec3 aPos; layout(location=1) in vec3 aTangent; layout(location=2) in vec3 aTV;
uniform mat4 uView,uProj,uLightSpace;
out vec3 vWorld; out vec3 vTan; out vec4 vLS; out float vT; out float vSide; out float vShade;
void main(){ vWorld=aPos; vTan=aTangent; vLS=uLightSpace*vec4(aPos,1.0); vT=aTV.x; vSide=aTV.y; vShade=aTV.z;
  gl_Position=uProj*uView*vec4(aPos,1.0); })";
static const char* FS_HAIR = R"(#version 330 core
in vec3 vWorld; in vec3 vTan; in vec4 vLS; in float vT; in float vSide; in float vShade; out vec4 frag;
uniform vec3 uCamPos,uLightDir,uLightColor; uniform vec3 uRootColor,uTipColor; uniform sampler2D uShadow;
float shadow(vec4 ls){ vec3 p=ls.xyz/ls.w*0.5+0.5; if(p.z>1.0)return 1.0; float b=0.0015; float s=0.0;
  vec2 tx=1.0/vec2(textureSize(uShadow,0)); for(int x=-1;x<=1;x++)for(int y=-1;y<=1;y++){float d=texture(uShadow,p.xy+vec2(x,y)*tx).r; s+=(p.z-b>d)?0.0:1.0;} return s/9.0; }
float strandSpec(vec3 T,vec3 V,vec3 L,float exps,float shift,vec3 N){ vec3 Ts=normalize(T+shift*N); vec3 H=normalize(L+V);
  float d=dot(Ts,H); float s=sqrt(max(1.0-d*d,0.0)); return pow(s,exps); }
void main(){ vec3 T=normalize(vTan); vec3 V=normalize(uCamPos-vWorld); vec3 L=normalize(-uLightDir);
  vec3 Nh=normalize(cross(cross(T,V),T));                 // view-perp component of tangent frame
  float TdL=dot(T,L); float diffuse=sqrt(max(1.0-TdL*TdL,0.0));  // Kajiya-Kay diffuse
  float sp1=strandSpec(T,V,L,110.0,0.10,Nh);              // primary sharp highlight (white)
  float sp2=strandSpec(T,V,L,18.0,-0.07,Nh);              // secondary soft highlight (tinted)
  float sh=shadow(vLS);
  vec3 base=mix(uRootColor,uTipColor,vT)*(1.0+vShade*0.22);   // per-strand colour jitter
  float ao=mix(0.42,1.0,vT);                              // darker/denser near the root
  vec3 ambient=base*vec3(0.34,0.37,0.44);                 // sky fill (NOT scaled by the sun)
  vec3 diff=base*diffuse*uLightColor*0.20*sh;
  vec3 spec=(sp1*(0.5+0.15*vShade)+sp2*base*1.3)*uLightColor*0.22*sh;
  vec3 col=(ambient+diff+spec)*ao;
  float edge=1.0-abs(vSide); float alpha=smoothstep(0.0,0.5,edge)*0.95;   // soft strand edges
  frag=vec4(col,alpha); })";

static const char* VS_FULL = R"(#version 330 core
out vec2 vUV; void main(){ vec2 p=vec2((gl_VertexID<<1)&2, gl_VertexID&2); vUV=p; gl_Position=vec4(p*2.0-1.0,0.0,1.0); })";
static const char* FS_TONE = R"(#version 330 core
in vec2 vUV; out vec4 frag; uniform sampler2D uHDR;
vec3 aces(vec3 x){return clamp((x*(2.51*x+0.03))/(x*(2.43*x+0.59)+0.14),0.0,1.0);}
void main(){ vec3 c=texture(uHDR,vUV).rgb; c=aces(c*1.05); c=pow(c,vec3(1.0/2.2)); frag=vec4(c,1.0); })";

// ------------------------------------------------------------------ renderer
struct Renderer {
    int W, H, shadowRes = 2048;
    GLuint msFBO = 0, msColor = 0, msDepth = 0;      // HDR + MSAA
    GLuint resFBO = 0, resTex = 0;                    // resolved HDR
    GLuint outFBO = 0, outTex = 0;                    // tonemapped RGB8
    GLuint shFBO = 0, shTex = 0;                      // shadow depth
    GLuint pSolid, pGround, pCloth, pDepth, pHair, pTone;
    GLuint emptyVAO = 0;
    glm::mat4 lightSpace;
    glm::vec3 lightDir = glm::normalize(glm::vec3(-0.5f, -1.0f, -0.35f));
    glm::vec3 lightColor = glm::vec3(1.0f, 0.97f, 0.9f) * 3.4f;

    void init(int w, int h) {
        W = w; H = h;
        glEnable(GL_DEPTH_TEST); glEnable(GL_MULTISAMPLE);
        auto full = std::string(PBR_HEAD);
        pSolid = program(VS_SCENE, (full + FS_SOLID).c_str());
        pGround = program(VS_SCENE, (full + FS_GROUND).c_str());
        pCloth = program(VS_SCENE, (full + FS_CLOTH).c_str());
        pDepth = program(VS_DEPTH, FS_DEPTH);
        pHair = program(VS_HAIR, FS_HAIR);
        pTone = program(VS_FULL, FS_TONE);
        glGenVertexArrays(1, &emptyVAO);

        glGenFramebuffers(1, &msFBO); glBindFramebuffer(GL_FRAMEBUFFER, msFBO);
        glGenRenderbuffers(1, &msColor); glBindRenderbuffer(GL_RENDERBUFFER, msColor);
        glRenderbufferStorageMultisample(GL_RENDERBUFFER, 4, GL_RGBA16F, W, H);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, msColor);
        glGenRenderbuffers(1, &msDepth); glBindRenderbuffer(GL_RENDERBUFFER, msDepth);
        glRenderbufferStorageMultisample(GL_RENDERBUFFER, 4, GL_DEPTH_COMPONENT24, W, H);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, msDepth);

        glGenFramebuffers(1, &resFBO); glBindFramebuffer(GL_FRAMEBUFFER, resFBO);
        glGenTextures(1, &resTex); glBindTexture(GL_TEXTURE_2D, resTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, W, H, 0, GL_RGBA, GL_FLOAT, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, resTex, 0);

        glGenFramebuffers(1, &outFBO); glBindFramebuffer(GL_FRAMEBUFFER, outFBO);
        glGenTextures(1, &outTex); glBindTexture(GL_TEXTURE_2D, outTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, W, H, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, outTex, 0);

        glGenFramebuffers(1, &shFBO); glBindFramebuffer(GL_FRAMEBUFFER, shFBO);
        glGenTextures(1, &shTex); glBindTexture(GL_TEXTURE_2D, shTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, shadowRes, shadowRes, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        float border[4] = {1, 1, 1, 1}; glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, border);
        glDrawBuffer(GL_NONE); glReadBuffer(GL_NONE);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, shTex, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }
    void setLightForScene(glm::vec3 centre, float radius) {
        glm::vec3 lpos = centre - lightDir * (radius * 2.2f);
        glm::mat4 lp = glm::ortho(-radius * 1.3f, radius * 1.3f, -radius * 1.3f, radius * 1.3f, 0.1f, radius * 5.0f);
        glm::mat4 lv = glm::lookAt(lpos, centre, glm::vec3(0, 1, 0));
        lightSpace = lp * lv;
    }
    // --- shadow pass ---
    void beginShadow() {
        glBindFramebuffer(GL_FRAMEBUFFER, shFBO); glViewport(0, 0, shadowRes, shadowRes);
        glClear(GL_DEPTH_BUFFER_BIT); glUseProgram(pDepth); setM4(pDepth, "uLightSpace", lightSpace);
        glEnable(GL_POLYGON_OFFSET_FILL); glPolygonOffset(2.0f, 4.0f);
    }
    void shadowDraw(const Mesh& m, const glm::mat4& model) { setM4(pDepth, "uModel", model); m.draw(); }
    void endShadow() { glDisable(GL_POLYGON_OFFSET_FILL); }
    // --- main scene pass ---
    void beginScene(const glm::mat4& view, const glm::mat4& proj, const glm::vec3& camPos) {
        glBindFramebuffer(GL_FRAMEBUFFER, msFBO); glViewport(0, 0, W, H);
        glClearColor(0.46f, 0.56f, 0.72f, 1.0f); glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        curView = view; curProj = proj; curCam = camPos;
    }
    void useProgram(GLuint p) {
        glUseProgram(p);
        setM4(p, "uView", curView); setM4(p, "uProj", curProj); setM4(p, "uLightSpace", lightSpace);
        setV3(p, "uCamPos", curCam); setV3(p, "uLightDir", lightDir); setV3(p, "uLightColor", lightColor);
        glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, shTex); setI(p, "uShadow", 0);
    }
    void drawPBR(GLuint p, const Mesh& m, const glm::mat4& model, glm::vec3 albedo, float rough, float metal, float ao = 1.0f) {
        useProgram(p);
        setM4(p, "uModel", model); setV3(p, "uAlbedo", albedo); setF(p, "uRough", rough); setF(p, "uMetal", metal); setF(p, "uAO", ao);
        m.draw();
    }
    // --- resolve + tonemap + present ---
    void endScene() {
        glBindFramebuffer(GL_READ_FRAMEBUFFER, msFBO); glBindFramebuffer(GL_DRAW_FRAMEBUFFER, resFBO);
        glBlitFramebuffer(0, 0, W, H, 0, 0, W, H, GL_COLOR_BUFFER_BIT, GL_LINEAR);
        glBindFramebuffer(GL_FRAMEBUFFER, outFBO); glViewport(0, 0, W, H);
        glDisable(GL_DEPTH_TEST); glUseProgram(pTone);
        glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, resTex); setI(pTone, "uHDR", 0);
        glBindVertexArray(emptyVAO); glDrawArrays(GL_TRIANGLES, 0, 3);
        glEnable(GL_DEPTH_TEST);
    }
    void present() {   // blit tonemapped output to the default framebuffer (window)
        glBindFramebuffer(GL_READ_FRAMEBUFFER, outFBO); glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glBlitFramebuffer(0, 0, W, H, 0, 0, W, H, GL_COLOR_BUFFER_BIT, GL_LINEAR);
    }
    bool screenshot(const char* path) {
        std::vector<unsigned char> px((size_t)W * H * 3);
        glBindFramebuffer(GL_READ_FRAMEBUFFER, outFBO);
        glPixelStorei(GL_PACK_ALIGNMENT, 1);
        glReadPixels(0, 0, W, H, GL_RGB, GL_UNSIGNED_BYTE, px.data());
        return writePNG(path, px.data(), W, H, true);
    }
    glm::mat4 curView{1}, curProj{1}; glm::vec3 curCam{0};
};

// ---------------------------------------------------------------- app window
struct App {
    GLFWwindow* win = nullptr;
    int W, H; bool headless;
    App(int w, int h, const char* title, bool headless) : W(w), H(h), headless(headless) {
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_SAMPLES, 4);
        glfwWindowHint(GLFW_VISIBLE, headless ? GLFW_FALSE : GLFW_TRUE);
        win = glfwCreateWindow(w, h, title, nullptr, nullptr);
        glfwMakeContextCurrent(win);
        glewExperimental = GL_TRUE; glewInit();
    }
    ~App() { if (win) glfwDestroyWindow(win); glfwTerminate(); }
    bool running() const { return !glfwWindowShouldClose(win); }
    void poll() { glfwPollEvents(); glfwSwapBuffers(win); }
};

} // namespace gfx
