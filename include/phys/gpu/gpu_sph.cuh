// GPU execution — a CUDA massively-parallel SPH fluid solver, the compute-backend
// parity item this header-only CPU engine was missing (cf. Brax/MJX, PhysX-GPU).
// Standard uniform-grid pipeline (Green, "Particle Simulation using CUDA"): hash
// particles into cells, radix-sort by cell (thrust), find per-cell ranges, then
// per-particle density → pressure → force → integrate, each a coalesced kernel over
// the sorted arrays. Runs hundreds of thousands of particles on the GPU in real time.
// Compile with nvcc (PTX target JITs to any newer GPU, e.g. sm_120 Blackwell).
#pragma once
#include <cuda_runtime.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <vector>
#include <cstdio>
#include <cmath>

namespace phys { namespace gpu {

__host__ __device__ inline float3 operator+(float3 a, float3 b) { return make_float3(a.x + b.x, a.y + b.y, a.z + b.z); }
__host__ __device__ inline float3 operator-(float3 a, float3 b) { return make_float3(a.x - b.x, a.y - b.y, a.z - b.z); }
__host__ __device__ inline float3 operator*(float3 a, float s)  { return make_float3(a.x * s, a.y * s, a.z * s); }
__host__ __device__ inline float  dot3(float3 a, float3 b)      { return a.x * b.x + a.y * b.y + a.z * b.z; }
__host__ __device__ inline float3 cross3(float3 a, float3 b)    { return make_float3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }

// up to this many static axis-aligned box obstacles (cliff, ledge, pool floor, …)
#define GPU_SPH_MAXBOX 8

struct Params {
    int n; float h, h2, mass, rho0, stiff, visc, xsph;
    float cPoly6, cSpiky, cVisc;
    float surfTens;                 // cohesion / surface tension (Becker–Teschner): holds free-falling streams & droplets together
    float3 gmin, gmax, grav; float cell; int3 grid; float rest, wallDamp;
    // --- solid box obstacles the fluid collides with (0 boxes → disabled) ---
    int nBox; float3 boxMin[GPU_SPH_MAXBOX], boxMax[GPU_SPH_MAXBOX]; float boxDamp, boxFric;
    // --- recirculating emitter: a particle that sinks past recycleY is teleported
    //     back into the emit AABB with emitVel (+jitter), so a fixed N sustains a
    //     continuous inflow — a waterfall source/sink without realloc (emitOn=0 → off) ---
    int emitOn; float3 emitMin, emitMax, emitVel; float recycleY;
    // --- rotating reference frame: Coriolis acceleration a = -2·omega×v (0 → inertial
    //     frame). Geophysical "f-plane": gravity stays vertical, no centrifugal term ---
    float3 omega;
};
__constant__ Params P;

// deterministic per-particle hash → 3 uniforms in [0,1) (seeded per step so respawns spread)
__device__ inline void hash3(unsigned s, float& a, float& b, float& c) {
    s ^= s >> 16; s *= 0x7feb352du; s ^= s >> 15; s *= 0x846ca68bu; s ^= s >> 16;
    a = (s & 0xffff) / 65535.0f;
    unsigned t = s * 2654435761u; b = (t & 0xffff) / 65535.0f;
    unsigned u = t * 40503u;      c = (u & 0xffff) / 65535.0f;
}
// push a particle out of a solid box along its least-penetration face; reflect the
// normal velocity (damped) and shear the tangential (friction)
__device__ inline void collideBox(float3& p, float3& v, float3 bmn, float3 bmx, float damp, float fric) {
    if (p.x <= bmn.x || p.x >= bmx.x || p.y <= bmn.y || p.y >= bmx.y || p.z <= bmn.z || p.z >= bmx.z) return;
    float dxl = p.x - bmn.x, dxh = bmx.x - p.x, dyl = p.y - bmn.y, dyh = bmx.y - p.y, dzl = p.z - bmn.z, dzh = bmx.z - p.z;
    float m = fminf(fminf(fminf(dxl, dxh), fminf(dyl, dyh)), fminf(dzl, dzh));
    if      (m == dxl) { p.x = bmn.x; v.x = -v.x * damp; v.y *= fric; v.z *= fric; }
    else if (m == dxh) { p.x = bmx.x; v.x = -v.x * damp; v.y *= fric; v.z *= fric; }
    else if (m == dyl) { p.y = bmn.y; v.y = -v.y * damp; v.x *= fric; v.z *= fric; }
    else if (m == dyh) { p.y = bmx.y; v.y = -v.y * damp; v.x *= fric; v.z *= fric; }
    else if (m == dzl) { p.z = bmn.z; v.z = -v.z * damp; v.x *= fric; v.y *= fric; }
    else               { p.z = bmx.z; v.z = -v.z * damp; v.x *= fric; v.y *= fric; }
}

__device__ inline int3 cellOf(float3 p) {
    return make_int3(min(max(int((p.x - P.gmin.x) / P.cell), 0), P.grid.x - 1),
                     min(max(int((p.y - P.gmin.y) / P.cell), 0), P.grid.y - 1),
                     min(max(int((p.z - P.gmin.z) / P.cell), 0), P.grid.z - 1));
}
__device__ inline int hashOf(int3 c) { return (c.z * P.grid.y + c.y) * P.grid.x + c.x; }

__global__ void calcHash(const float4* pos, int* hash, int* idx) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; if (i >= P.n) return;
    hash[i] = hashOf(cellOf(make_float3(pos[i].x, pos[i].y, pos[i].z))); idx[i] = i;
}
__global__ void reorder(const int* hash, const int* idx, const float4* oPos, const float4* oVel,
                        float4* sPos, float4* sVel, int* cellStart, int* cellEnd) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; if (i >= P.n) return;
    int h = hash[i];
    if (i == 0 || h != hash[i - 1]) cellStart[h] = i;
    if (i == P.n - 1 || h != hash[i + 1]) cellEnd[h] = i + 1;
    int j = idx[i]; sPos[i] = oPos[j]; sVel[i] = oVel[j];
}
__global__ void densityK(const float4* sPos, const int* cellStart, const int* cellEnd, float* dens, float* pres) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; if (i >= P.n) return;
    float3 pi = make_float3(sPos[i].x, sPos[i].y, sPos[i].z); int3 c = cellOf(pi); float d = 0;
    for (int dz = -1; dz <= 1; dz++) for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
        int3 cc = make_int3(c.x + dx, c.y + dy, c.z + dz);
        if (cc.x < 0 || cc.y < 0 || cc.z < 0 || cc.x >= P.grid.x || cc.y >= P.grid.y || cc.z >= P.grid.z) continue;
        int h = hashOf(cc), s = cellStart[h]; if (s == 0x7fffffff) continue; int e = cellEnd[h];
        for (int j = s; j < e; j++) { float3 r = pi - make_float3(sPos[j].x, sPos[j].y, sPos[j].z);
            float r2 = dot3(r, r); if (r2 < P.h2) { float t = P.h2 - r2; d += P.mass * P.cPoly6 * t * t * t; } }
    }
    dens[i] = d > 1e-6f ? d : 1e-6f;
    float p = P.stiff * (dens[i] - P.rho0); pres[i] = p > 0 ? p : 0;      // clamp ≥0 (no tensile clumping)
}
__global__ void forceK(const float4* sPos, const float4* sVel, const float* dens, const float* pres,
                       const int* cellStart, const int* cellEnd, float4* acc) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; if (i >= P.n) return;
    float3 pi = make_float3(sPos[i].x, sPos[i].y, sPos[i].z), vi = make_float3(sVel[i].x, sVel[i].y, sVel[i].z);
    int3 c = cellOf(pi); float3 fp = make_float3(0, 0, 0), fv = make_float3(0, 0, 0), fc = make_float3(0, 0, 0);
    for (int dz = -1; dz <= 1; dz++) for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) {
        int3 cc = make_int3(c.x + dx, c.y + dy, c.z + dz);
        if (cc.x < 0 || cc.y < 0 || cc.z < 0 || cc.x >= P.grid.x || cc.y >= P.grid.y || cc.z >= P.grid.z) continue;
        int h = hashOf(cc), s = cellStart[h]; if (s == 0x7fffffff) continue; int e = cellEnd[h];
        for (int j = s; j < e; j++) { if (j == i) continue;
            float3 r = pi - make_float3(sPos[j].x, sPos[j].y, sPos[j].z); float rl = sqrtf(dot3(r, r));
            if (rl >= P.h || rl < 1e-7f) continue;
            float3 rn = r * (1.0f / rl);
            fp = fp + rn * (-P.mass * (pres[i] + pres[j]) / (2 * dens[j]) * P.cSpiky * (P.h - rl) * (P.h - rl));
            float3 dv = make_float3(sVel[j].x - vi.x, sVel[j].y - vi.y, sVel[j].z - vi.z);
            fv = fv + dv * (P.visc * P.mass / dens[j] * P.cVisc * (P.h - rl));
            if (P.surfTens > 0) { float hr = P.h2 - rl * rl; float w6 = P.cPoly6 * hr * hr * hr;   // cohesion: attract toward neighbours (poly6)
                fc = fc + r * (-P.surfTens * P.mass * w6); }
        }
    }
    float3 a = (fp + fv + fc) * (1.0f / dens[i]) + P.grav;
    acc[i] = make_float4(a.x, a.y, a.z, 0);
}
__global__ void integrateK(float4* sPos, float4* sVel, const float4* acc, float dt, unsigned seed) {
    int i = blockIdx.x * blockDim.x + threadIdx.x; if (i >= P.n) return;
    float dye = sPos[i].w;                                               // advected passive scalar (dye/tracer), rides in pos.w
    float3 v = make_float3(sVel[i].x + acc[i].x * dt, sVel[i].y + acc[i].y * dt, sVel[i].z + acc[i].z * dt);
    float ow = sqrtf(dot3(P.omega, P.omega));                            // Coriolis: exact rotation of v about the omega axis (magnitude-preserving)
    if (ow > 1e-8f) { float3 k = P.omega * (1.0f / ow); float th = -2.0f * ow * dt, c = cosf(th), s = sinf(th);
        v = v * c + cross3(k, v) * s + k * (dot3(k, v) * (1.0f - c)); }
    float3 p = make_float3(sPos[i].x + v.x * dt, sPos[i].y + v.y * dt, sPos[i].z + v.z * dt);
    for (int b = 0; b < P.nBox; b++) collideBox(p, v, P.boxMin[b], P.boxMax[b], P.boxDamp, P.boxFric);   // solid cliff/ledge/floor
    float* pp = &p.x; float* vv = &v.x; const float* lo = &P.gmin.x; const float* hi = &P.gmax.x;
    for (int a = 0; a < 3; a++) {                                        // domain bounds with damping
        if (pp[a] < lo[a]) { pp[a] = lo[a]; vv[a] = -vv[a] * P.wallDamp; }
        if (pp[a] > hi[a]) { pp[a] = hi[a]; vv[a] = -vv[a] * P.wallDamp; }
    }
    if (P.emitOn && p.y < P.recycleY) {                                 // drained → respawn at the source (fixed-N recirculation)
        float a0, a1, a2; hash3((unsigned)i * 2654435761u ^ seed, a0, a1, a2);
        p = make_float3(P.emitMin.x + (P.emitMax.x - P.emitMin.x) * a0,
                        P.emitMin.y + (P.emitMax.y - P.emitMin.y) * a1,
                        P.emitMin.z + (P.emitMax.z - P.emitMin.z) * a2);
        v = make_float3(P.emitVel.x + (a0 - 0.5f) * 0.25f, P.emitVel.y, P.emitVel.z + (a2 - 0.5f) * 0.25f);
    }
    sPos[i] = make_float4(p.x, p.y, p.z, dye); sVel[i] = make_float4(v.x, v.y, v.z, 0);   // keep the dye scalar with the particle
}

class GpuSPH {
public:
    int n = 0; Params p{};
    float4 *dPos = 0, *dVel = 0, *dSort = 0, *dSortV = 0, *dAcc = 0;
    float *dDens = 0, *dPres = 0; int *dHash = 0, *dIdx = 0, *dCellStart = 0, *dCellEnd = 0;
    int numCells = 0;

    void init(const std::vector<float3>& pos, Params params, const std::vector<float>* dye = nullptr) {
        n = (int)pos.size(); p = params; p.n = n;
        p.h2 = p.h * p.h;
        p.cPoly6 = 315.0f / (64.0f * 3.14159265f * powf(p.h, 9));
        p.cSpiky = -45.0f / (3.14159265f * powf(p.h, 6));
        p.cVisc = 45.0f / (3.14159265f * powf(p.h, 6));
        p.cell = p.h;
        p.grid = make_int3(int((p.gmax.x - p.gmin.x) / p.cell) + 1, int((p.gmax.y - p.gmin.y) / p.cell) + 1, int((p.gmax.z - p.gmin.z) / p.cell) + 1);
        numCells = p.grid.x * p.grid.y * p.grid.z;
        std::vector<float4> h4(n); for (int i = 0; i < n; i++) h4[i] = make_float4(pos[i].x, pos[i].y, pos[i].z, dye ? (*dye)[i] : 0.0f);
        cudaMalloc(&dPos, n * sizeof(float4)); cudaMalloc(&dVel, n * sizeof(float4));
        cudaMalloc(&dSort, n * sizeof(float4)); cudaMalloc(&dSortV, n * sizeof(float4)); cudaMalloc(&dAcc, n * sizeof(float4));
        cudaMalloc(&dDens, n * sizeof(float)); cudaMalloc(&dPres, n * sizeof(float));
        cudaMalloc(&dHash, n * sizeof(int)); cudaMalloc(&dIdx, n * sizeof(int));
        cudaMalloc(&dCellStart, numCells * sizeof(int)); cudaMalloc(&dCellEnd, numCells * sizeof(int));
        cudaMemcpy(dPos, h4.data(), n * sizeof(float4), cudaMemcpyHostToDevice);
        cudaMemset(dVel, 0, n * sizeof(float4));
        // measure rho0 from the initial packing so the fluid starts near zero pressure
        p.rho0 = 1.0f; cudaMemcpyToSymbol(P, &p, sizeof(Params));
        buildGrid(); densityK<<<blocks(), 256>>>(dSort, dCellStart, dCellEnd, dDens, dPres);
        std::vector<float> hd(n); cudaMemcpy(hd.data(), dDens, n * sizeof(float), cudaMemcpyDeviceToHost);
        double s = 0; for (float d : hd) s += d; p.rho0 = (float)(s / n);
        cudaMemcpyToSymbol(P, &p, sizeof(Params));
    }
    int blocks() const { return (n + 255) / 256; }
    void buildGrid() {
        calcHash<<<blocks(), 256>>>(dPos, dHash, dIdx);
        thrust::sort_by_key(thrust::device_ptr<int>(dHash), thrust::device_ptr<int>(dHash + n), thrust::device_ptr<int>(dIdx));
        cudaMemset(dCellStart, 0x7f, numCells * sizeof(int));
        reorder<<<blocks(), 256>>>(dHash, dIdx, dPos, dVel, dSort, dSortV, dCellStart, dCellEnd);
    }
    unsigned frame = 0;
    void step(float dt) {
        cudaMemcpyToSymbol(P, &p, sizeof(Params));                       // this instance's params → supports several GpuSPH sharing __constant__ P
        buildGrid();
        densityK<<<blocks(), 256>>>(dSort, dCellStart, dCellEnd, dDens, dPres);
        forceK<<<blocks(), 256>>>(dSort, dSortV, dDens, dPres, dCellStart, dCellEnd, dAcc);
        integrateK<<<blocks(), 256>>>(dSort, dSortV, dAcc, dt, frame++ * 2654435761u + 12345u);
        std::swap(dPos, dSort); std::swap(dVel, dSortV);                 // sorted arrays become current
    }
    void computeDensityNow() {                                          // grid + density on the current positions
        buildGrid(); densityK<<<blocks(), 256>>>(dSort, dCellStart, dCellEnd, dDens, dPres); cudaDeviceSynchronize();
    }
    double densitySum() { std::vector<float> h(n); cudaMemcpy(h.data(), dDens, n * sizeof(float), cudaMemcpyDeviceToHost); double s = 0; for (float d : h) s += d; return s; }
    double densityMean() { return densitySum() / n; }
    void download(std::vector<float3>& out) {
        std::vector<float4> h4(n); cudaMemcpy(h4.data(), dPos, n * sizeof(float4), cudaMemcpyDeviceToHost);
        out.resize(n); for (int i = 0; i < n; i++) out[i] = make_float3(h4[i].x, h4[i].y, h4[i].z);
    }
    // positions + speed (for shading whitewater vs. calm pool). Speed via |vel|.
    void downloadPosSpeed(std::vector<float4>& out) {
        std::vector<float4> hp(n), hv(n);
        cudaMemcpy(hp.data(), dPos, n * sizeof(float4), cudaMemcpyDeviceToHost);
        cudaMemcpy(hv.data(), dVel, n * sizeof(float4), cudaMemcpyDeviceToHost);
        out.resize(n);
        for (int i = 0; i < n; i++) out[i] = make_float4(hp[i].x, hp[i].y, hp[i].z, sqrtf(hv[i].x * hv[i].x + hv[i].y * hv[i].y + hv[i].z * hv[i].z));
    }
    // per-particle SPH density (low ⇒ airborne/aerated spray, high ⇒ dense bulk water)
    void downloadDens(std::vector<float>& out) { out.resize(n); cudaMemcpy(out.data(), dDens, n * sizeof(float), cudaMemcpyDeviceToHost); }
    // positions + advected dye scalar (pos.w) — for Lagrangian tracer visualisation
    void downloadPosDye(std::vector<float4>& out) { out.resize(n); cudaMemcpy(out.data(), dPos, n * sizeof(float4), cudaMemcpyDeviceToHost); }
    void downloadVel(std::vector<float3>& out) {
        std::vector<float4> h(n); cudaMemcpy(h.data(), dVel, n * sizeof(float4), cudaMemcpyDeviceToHost);
        out.resize(n); for (int i = 0; i < n; i++) out[i] = make_float3(h[i].x, h[i].y, h[i].z);
    }
    void setVel(const std::vector<float3>& vel) {                       // seed initial velocities (init zeroes them)
        std::vector<float4> h(n); for (int i = 0; i < n; i++) h[i] = make_float4(vel[i].x, vel[i].y, vel[i].z, 0);
        cudaMemcpy(dVel, h.data(), n * sizeof(float4), cudaMemcpyHostToDevice);
    }
    // CPU reference poly6 density at a host-space point (for correctness tests)
    static float cpuDensity(const std::vector<float3>& pos, int i, float h, float mass) {
        float h2 = h * h, c = 315.0f / (64.0f * 3.14159265f * powf(h, 9)), d = 0;
        for (size_t j = 0; j < pos.size(); j++) { float3 r = pos[i] - pos[j]; float r2 = dot3(r, r);
            if (r2 < h2) { float t = h2 - r2; d += mass * c * t * t * t; } }
        return d;
    }
    void free() {
        for (void* q : {(void*)dPos, (void*)dVel, (void*)dSort, (void*)dSortV, (void*)dAcc, (void*)dDens, (void*)dPres, (void*)dHash, (void*)dIdx, (void*)dCellStart, (void*)dCellEnd}) cudaFree(q);
    }
};

}} // namespace phys::gpu
