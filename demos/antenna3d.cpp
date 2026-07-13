// A dipole antenna in 3D. The current distribution and far-field are solved by the
// Method of Moments (phys::em::Dipole); the far-field |F(θ)| is drawn as a glowing
// 3-D radiation-pattern surface (the classic dipole "donut"), the wire glows with
// its MoM current, and expanding wavefront shells — modulated by the true pattern —
// show the field radiating outward. The dipole length is swept from 0.5λ to 1.5λ so
// the pattern morphs from the single donut to a multi-lobe pattern, all straight
// from the MoM solution.
//
//   ./antenna3d --video out/f 360   |   ./antenna3d --shot frame.png 90
#include "phys/mom.h"
#include "common/gfx.h"
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstring>
using namespace phys::em;

static const int PT = 91;                       // pattern samples over θ∈[0,π]

// solve the dipole and return a normalized gain profile g(θ)=|F(θ)|/max, θ index 0..PT-1
static void patternProfile(double L, Dipole& dip, std::vector<float>& g) {
    dip.L = L; dip.a = 0.004; dip.N = 51; dip.solve();
    g.assign(PT, 0.f); double mx = 1e-12;
    for (int t = 0; t < PT; t++) { double th = M_PI * t / (PT - 1); double v = dip.pattern(th); g[t] = (float)v; mx = std::max(mx, v); }
    for (auto& x : g) x /= (float)mx;
}
static float sampleG(const std::vector<float>& g, double theta) {
    double f = theta / M_PI * (PT - 1); int i = (int)f; if (i < 0) i = 0; if (i > PT - 2) i = PT - 2;
    double a = f - i; return g[i] * (1 - a) + g[i + 1] * a;
}

// build a parametric radiation-pattern surface  p(θ,φ)=r(θ)·dir, dir along +y axis
static void buildPatternMesh(const std::vector<float>& g, float scale, int NT, int NP,
                             std::vector<float>& verts, std::vector<unsigned>& idx, bool wantIdx) {
    std::vector<glm::vec3> P((NT + 1) * (NP + 1));
    for (int it = 0; it <= NT; it++) {
        double th = M_PI * it / NT; double r = 0.06 + sampleG(g, th) * scale;
        for (int ip = 0; ip <= NP; ip++) {
            double ph = 2 * M_PI * ip / NP;
            glm::vec3 dir((float)(std::sin(th) * std::cos(ph)), (float)std::cos(th), (float)(std::sin(th) * std::sin(ph)));
            P[it * (NP + 1) + ip] = dir * (float)r;
        }
    }
    // smooth normals from neighbouring grid points
    std::vector<glm::vec3> Nrm((NT + 1) * (NP + 1), glm::vec3(0));
    auto at = [&](int a, int b) { return P[a * (NP + 1) + b]; };
    for (int it = 0; it <= NT; it++) for (int ip = 0; ip <= NP; ip++) {
        glm::vec3 du = at(std::min(it + 1, NT), ip) - at(std::max(it - 1, 0), ip);
        glm::vec3 dv = at(it, ip < NP ? ip + 1 : 1) - at(it, ip > 0 ? ip - 1 : NP - 1);
        glm::vec3 n = glm::cross(dv, du); float l = glm::length(n); Nrm[it * (NP + 1) + ip] = l > 1e-9f ? n / l : glm::vec3(0, 1, 0);
    }
    verts.clear();
    for (int k = 0; k < (int)P.size(); k++) verts.insert(verts.end(), { P[k].x, P[k].y, P[k].z, Nrm[k].x, Nrm[k].y, Nrm[k].z });
    if (wantIdx) { idx.clear();
        for (int it = 0; it < NT; it++) for (int ip = 0; ip < NP; ip++) {
            unsigned a = it * (NP + 1) + ip, b = a + 1, c = a + (NP + 1), d = c + 1;
            idx.insert(idx.end(), { a, c, b, b, c, d });
        } }
}

int main(int argc, char** argv) {
    const char* video = nullptr; const char* shot = nullptr; int frames = 360;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--video")) { video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--shot")) { shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }
    const int W = 1280, H = 800;
    gfx::App app(W, H, "antenna3d", true);

    // --- MSAA offscreen target ---
    GLuint msFBO, msC, msD, rFBO, rC;
    glGenFramebuffers(1, &msFBO); glBindFramebuffer(GL_FRAMEBUFFER, msFBO);
    glGenRenderbuffers(1, &msC); glBindRenderbuffer(GL_RENDERBUFFER, msC); glRenderbufferStorageMultisample(GL_RENDERBUFFER, 4, GL_RGBA8, W, H);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, msC);
    glGenRenderbuffers(1, &msD); glBindRenderbuffer(GL_RENDERBUFFER, msD); glRenderbufferStorageMultisample(GL_RENDERBUFFER, 4, GL_DEPTH_COMPONENT24, W, H);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, msD);
    glGenFramebuffers(1, &rFBO); glBindFramebuffer(GL_FRAMEBUFFER, rFBO);
    glGenTextures(1, &rC); glBindTexture(GL_TEXTURE_2D, rC); glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, W, H, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, rC, 0);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // --- shaders ---
    const char* VS = R"(#version 330 core
layout(location=0) in vec3 aPos; layout(location=1) in vec3 aN;
uniform mat4 uM,uV,uP; out vec3 vW; out vec3 vN;
void main(){ vec4 w=uM*vec4(aPos,1.0); vW=w.xyz; vN=mat3(uM)*aN; gl_Position=uP*uV*w; })";

    // radiation-pattern surface: colour by gain, lit + Fresnel rim glow, two-sided
    GLuint pPat = gfx::program(VS, R"(#version 330 core
in vec3 vW; in vec3 vN; out vec4 frag;
uniform vec3 uCam,uLight; uniform float uScale;
vec3 heat(float t){ t=clamp(t,0.0,1.0);
  vec3 c=mix(vec3(0.15,0.05,0.55),vec3(0.1,0.55,0.9),smoothstep(0.0,0.35,t));
  c=mix(c,vec3(0.1,0.85,0.45),smoothstep(0.3,0.6,t));
  c=mix(c,vec3(0.95,0.85,0.2),smoothstep(0.55,0.8,t));
  c=mix(c,vec3(1.0,0.35,0.2),smoothstep(0.78,1.0,t)); return c; }
void main(){ float gain=clamp(length(vW)/uScale,0.0,1.0);
  vec3 base=heat(gain); vec3 N=normalize(vN); vec3 V=normalize(uCam-vW); if(dot(N,V)<0.0) N=-N;
  vec3 L=normalize(uLight); float dif=max(dot(N,L),0.0);
  float rim=pow(1.0-max(dot(N,V),0.0),2.2);
  vec3 col=base*(0.30+0.70*dif)+base*rim*1.5+base*0.30;
  frag=vec4(col,1.0); })");

    // dipole wire: emissive gold, brightness = MoM |I| along its length
    GLuint pRod = gfx::program(VS, R"(#version 330 core
in vec3 vW; in vec3 vN; out vec4 frag;
uniform vec3 uCam; uniform float uHalf; uniform float uCur[64]; uniform int uNC;
void main(){ float f=clamp((vW.y+uHalf)/(2.0*uHalf),0.0,1.0); float x=f*(uNC-1); int i=int(x); i=clamp(i,0,uNC-2);
  float cur=mix(uCur[i],uCur[i+1],x-float(i));
  vec3 gold=vec3(1.0,0.82,0.35); float g=0.35+0.9*cur;
  vec3 V=normalize(uCam-vW); vec3 N=normalize(vN); float rim=pow(1.0-abs(dot(N,V)),2.0);
  frag=vec4(gold*g+vec3(1.0,0.9,0.6)*rim*0.6,1.0); })");

    // expanding wavefront shells: additive, brightness follows the true pattern g(θ),
    // colour alternates crest/trough of the outgoing wave  cos(k r − ω t)
    GLuint pShell = gfx::program(VS, R"(#version 330 core
in vec3 vW; out vec4 frag;
uniform float uR,uPhase,uFade; uniform float uPat[91];
void main(){ vec3 d=normalize(vW); float th=acos(clamp(d.y,-1.0,1.0));
  float fp=th/3.14159265*90.0; int i=int(fp); i=clamp(i,0,89); float gain=mix(uPat[i],uPat[i+1],fp-float(i));
  float c=cos(uPhase);
  float band=pow(max(c,0.0),12.0);          // thin outgoing crest (warm)
  float bandc=pow(max(-c,0.0),12.0);        // thin trough (cool)
  float amp=gain*gain*uFade;
  vec3 col=(vec3(1.0,0.5,0.2)*band + vec3(0.25,0.55,1.0)*bandc)*amp;
  frag=vec4(col,1.0); })");

    // background gradient (fullscreen)
    GLuint pBg = gfx::program(R"(#version 330 core
const vec2 v[3]=vec2[3](vec2(-1,-1),vec2(3,-1),vec2(-1,3)); out vec2 uv;
void main(){ uv=v[gl_VertexID]; gl_Position=vec4(v[gl_VertexID],0,1); })", R"(#version 330 core
in vec2 uv; out vec4 frag;
void main(){ float r=length(uv*vec2(1.1,1.0)); vec3 c=mix(vec3(0.05,0.06,0.10),vec3(0.01,0.01,0.02),clamp(r,0.0,1.0)); frag=vec4(c,1.0); })");

    auto setM = [](GLuint p, const char* n, const glm::mat4& m) { glUniformMatrix4fv(glGetUniformLocation(p, n), 1, GL_FALSE, &m[0][0]); };
    auto setV = [](GLuint p, const char* n, const glm::vec3& v) { glUniform3fv(glGetUniformLocation(p, n), 1, &v[0]); };
    auto setF = [](GLuint p, const char* n, float f) { glUniform1f(glGetUniformLocation(p, n), f); };

    // meshes
    const int NT = 140, NP = 200; const float SCALE = 1.7f;
    Dipole dip; std::vector<float> g, pv; std::vector<unsigned> pidx;
    patternProfile(0.5, dip, g); buildPatternMesh(g, SCALE, NT, NP, pv, pidx, true);
    gfx::Mesh pat; pat.upload(pv, pidx);
    gfx::Mesh rod = gfx::makeBox();
    gfx::Mesh sph = gfx::makeSphere(64, 40);
    GLuint bgVAO; glGenVertexArrays(1, &bgVAO);

    const int NSH = 90; const float Rin = 0.35f, Rout = 5.6f;
    const float kWave = 2.0f * M_PI * 5.0f / (Rout - Rin);   // ~5 crests across the span (densely sampled)
    glm::mat4 proj = glm::perspective(glm::radians(37.0f), (float)W / H, 0.05f, 100.0f);
    gfx::OrbitCamera cam; cam.target = glm::vec3(0, 0, 0); cam.dist = 8.2f; cam.pitch = 0.40f;

    auto renderFrame = [&](int f, float tphase, double L) {
        // update MoM pattern + current for this frame's length
        patternProfile(L, dip, g);
        buildPatternMesh(g, SCALE, NT, NP, pv, pidx, false); pat.updateVerts(pv);
        float cur[64]; int NC = std::min(64, dip.N); double Imax = dip.maxCurrent();
        for (int i = 0; i < NC; i++) cur[i] = (float)(std::abs(dip.I[i * dip.N / NC]) / Imax);
        float half = 0.9f * (float)L * 2.0f;                 // wire half-length (visual)

        cam.yaw = 0.5f + f * 0.0075f;
        glm::mat4 view = cam.view(); glm::vec3 eye = cam.eye();
        glm::vec3 light = glm::normalize(glm::vec3(0.5f, 0.8f, 0.6f));

        glBindFramebuffer(GL_FRAMEBUFFER, msFBO); glViewport(0, 0, W, H);
        glClearColor(0, 0, 0, 1); glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        // background
        glDisable(GL_DEPTH_TEST); glUseProgram(pBg); glBindVertexArray(bgVAO); glDrawArrays(GL_TRIANGLES, 0, 3);
        glEnable(GL_DEPTH_TEST); glDepthMask(GL_TRUE); glDisable(GL_BLEND); glDisable(GL_CULL_FACE);
        glm::mat4 I(1);
        // pattern donut
        glUseProgram(pPat); setM(pPat, "uM", I); setM(pPat, "uV", view); setM(pPat, "uP", proj);
        setV(pPat, "uCam", eye); setV(pPat, "uLight", light); setF(pPat, "uScale", SCALE);
        pat.draw();
        // dipole wire (thin tall box along y)
        glm::mat4 rodM = glm::scale(glm::mat4(1), glm::vec3(0.018f, half, 0.018f));
        glUseProgram(pRod); setM(pRod, "uM", rodM); setM(pRod, "uV", view); setM(pRod, "uP", proj);
        setV(pRod, "uCam", eye); setF(pRod, "uHalf", half);
        glUniform1fv(glGetUniformLocation(pRod, "uCur"), NC, cur); glUniform1i(glGetUniformLocation(pRod, "uNC"), NC);
        rod.draw();
        // radiating wavefront shells (additive, pattern-modulated)
        glEnable(GL_BLEND); glBlendFunc(GL_ONE, GL_ONE); glDepthMask(GL_FALSE);
        glUseProgram(pShell); setM(pShell, "uV", view); setM(pShell, "uP", proj);
        glUniform1fv(glGetUniformLocation(pShell, "uPat"), PT, g.data());
        for (int s = 0; s < NSH; s++) {
            float r = Rin + (s + 0.5f) / NSH * (Rout - Rin);
            float env = std::sin(M_PI * (s + 0.5f) / NSH);   // fade in/out at the ends
            float fade = env * (Rin / r) * 0.22f;             // 1/r amplitude decay (dim ambience)
            glm::mat4 M = glm::scale(glm::mat4(1), glm::vec3(r));
            setM(pShell, "uM", M); setF(pShell, "uR", r); setF(pShell, "uFade", fade);
            setF(pShell, "uPhase", kWave * r - tphase);
            sph.draw();
        }
        glDepthMask(GL_TRUE); glDisable(GL_BLEND);

        // resolve MSAA and read back
        glBindFramebuffer(GL_READ_FRAMEBUFFER, msFBO); glBindFramebuffer(GL_DRAW_FRAMEBUFFER, rFBO);
        glBlitFramebuffer(0, 0, W, H, 0, 0, W, H, GL_COLOR_BUFFER_BIT, GL_LINEAR);
    };
    auto save = [&](const char* path) {
        std::vector<unsigned char> px((size_t)W * H * 3);
        glBindFramebuffer(GL_READ_FRAMEBUFFER, rFBO); glPixelStorei(GL_PACK_ALIGNMENT, 1);
        glReadPixels(0, 0, W, H, GL_RGB, GL_UNSIGNED_BYTE, px.data());
        gfx::writePNG(path, px.data(), W, H, true);
    };
    // dipole length schedule: hold half-wave, sweep to 1.5λ, hold
    auto lenAt = [&](int f) {
        double u = (double)f / std::max(1, frames - 1);
        if (u < 0.35) return 0.5;
        if (u < 0.75) return 0.5 + (u - 0.35) / 0.40 * 1.0;   // 0.5 → 1.5
        return 1.5;
    };

    if (video) {
        for (int f = 0; f < frames; f++) { renderFrame(f, f * 0.5f, lenAt(f)); char q[512]; std::snprintf(q, sizeof(q), "%s_%04d.png", video, f); save(q); }
        std::printf("wrote %d frames\n", frames); return 0;
    }
    int fr = shot ? frames : 90;
    renderFrame(fr, fr * 0.5f, lenAt(fr)); save(shot ? shot : "antenna3d.png");
    std::printf("wrote shot at L=%.2f\n", lenAt(fr));
    return 0;
}
