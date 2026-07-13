// GPU SPH: the CUDA solver's grid density must match a CPU brute-force reference,
// the sim must stay finite/in-bounds and near rest density, and it must scale to
// hundreds of thousands of particles on the GPU.  Build: nvcc.
#include "phys/gpu/gpu_sph.cuh"
#include "check.h"
#include <vector>
#include <cmath>
#include <chrono>
using namespace phys::gpu;

static std::vector<float3> lattice(float3 lo, float3 hi, float s) {
    std::vector<float3> v;
    for (float x = lo.x; x <= hi.x; x += s) for (float y = lo.y; y <= hi.y; y += s) for (float z = lo.z; z <= hi.z; z += s)
        v.push_back(make_float3(x, y, z));
    return v;
}
static Params baseParams(float h) {
    Params p{}; p.h = h; p.mass = 0.02f; p.rho0 = 1.f; p.stiff = 250.f; p.visc = 4.f; p.xsph = 0;
    p.gmin = make_float3(-0.6f, 0.f, -0.4f); p.gmax = make_float3(0.6f, 1.4f, 0.4f);
    p.grav = make_float3(0, -9.81f, 0); p.wallDamp = 0.3f; return p;
}

int main() {
    int dev; if (cudaGetDevice(&dev) != cudaSuccess) { std::printf("no CUDA device\n"); return 0; }
    cudaDeviceProp prop; cudaGetDeviceProperties(&prop, dev);

    // A) GPU grid density == CPU brute-force density (same neighbour set within h)
    {
        float h = 0.10f;
        auto pos = lattice(make_float3(-0.3f, 0.1f, -0.2f), make_float3(0.3f, 0.5f, 0.2f), 0.05f);
        GpuSPH sph; sph.init(pos, baseParams(h));
        sph.computeDensityNow();
        double gpuSum = sph.densitySum();
        double cpuSum = 0; for (size_t i = 0; i < pos.size(); i++) cpuSum += GpuSPH::cpuDensity(pos, (int)i, h, sph.p.mass);
        CHECK(pos.size() > 500);
        CHECK(std::isfinite(gpuSum) && gpuSum > 0);
        CHECK_NEAR(gpuSum, cpuSum, cpuSum * 1e-3);                       // GPU matches CPU to float tolerance
        sph.free();
    }

    // B) a dropped block settles: finite, in-bounds, particle count conserved, near rest density
    {
        float h = 0.055f;
        auto pos = lattice(make_float3(-0.5f, 0.05f, -0.35f), make_float3(0.5f, 0.7f, 0.35f), 0.03f);
        int n0 = (int)pos.size();
        GpuSPH sph; sph.init(pos, baseParams(h));
        for (int s = 0; s < 400; s++) sph.step(1.2e-3f);
        std::vector<float3> out; sph.download(out);
        CHECK((int)out.size() == n0);
        bool finite = true, inside = true;
        for (auto& q : out) {
            if (!std::isfinite(q.x) || !std::isfinite(q.y) || !std::isfinite(q.z)) finite = false;
            if (q.x < sph.p.gmin.x - 1e-3f || q.x > sph.p.gmax.x + 1e-3f || q.y < sph.p.gmin.y - 1e-3f || q.y > sph.p.gmax.y + 1e-3f) inside = false;
        }
        CHECK(finite);
        CHECK(inside);
        sph.computeDensityNow();
        double mean = sph.densityMean();
        CHECK(mean > 0.8 * sph.p.rho0 && mean < 1.25 * sph.p.rho0);      // stays near rest density
        sph.free();
    }

    // C) GPU scale: hundreds of thousands of particles, real-time throughput
    {
        float h = 0.02f;
        Params p = baseParams(h);
        p.gmin = make_float3(-1.0f, 0.f, -1.0f); p.gmax = make_float3(1.0f, 2.0f, 1.0f);
        auto pos = lattice(make_float3(-0.9f, 0.05f, -0.9f), make_float3(0.9f, 0.9f, 0.9f), 0.02f);
        GpuSPH sph; sph.init(pos, p);
        cudaDeviceSynchronize();
        auto t0 = std::chrono::high_resolution_clock::now();
        const int STEPS = 60; for (int s = 0; s < STEPS; s++) sph.step(1.0e-3f);
        cudaDeviceSynchronize();
        double ms = std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - t0).count();
        double perStep = ms / STEPS;
        std::printf("  [gpu] %s (%d SMs): %d particles, %.2f ms/step, %.1fM particle-steps/s\n",
                    prop.name, prop.multiProcessorCount, sph.n, perStep, sph.n * STEPS / (ms / 1000.0) / 1e6);
        CHECK(sph.n > 200000);                                          // genuinely large scale
        CHECK(perStep < 60.0);                                          // interactive on the GPU
        std::vector<float3> out; sph.download(out); bool ok = true;
        for (auto& q : out) if (!std::isfinite(q.y)) ok = false;
        CHECK(ok);
        sph.free();
    }

    return test::report("gpu");
}
