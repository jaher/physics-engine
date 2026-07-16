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

    // D) solid box obstacle: fluid dropped onto a raised platform settles ON it, never inside it
    {
        Params p = baseParams(0.05f);
        p.gmin = make_float3(-0.5f, 0.f, -0.5f); p.gmax = make_float3(0.5f, 1.2f, 0.5f);
        p.nBox = 1; p.boxDamp = 0.2f; p.boxFric = 0.5f;
        p.boxMin[0] = make_float3(-0.25f, 0.20f, -0.25f); p.boxMax[0] = make_float3(0.25f, 0.50f, 0.25f);   // rock platform
        auto pos = lattice(make_float3(-0.20f, 0.55f, -0.20f), make_float3(0.20f, 0.95f, 0.20f), 0.035f);
        GpuSPH sph; sph.init(pos, p);
        for (int s = 0; s < 400; s++) sph.step(1.0e-3f);
        std::vector<float3> out; sph.download(out);
        int inside = 0; bool finite = true;
        for (auto& q : out) {
            if (!std::isfinite(q.y)) finite = false;
            if (q.x > p.boxMin[0].x + 1e-3f && q.x < p.boxMax[0].x - 1e-3f && q.y > p.boxMin[0].y + 3e-3f &&
                q.y < p.boxMax[0].y - 3e-3f && q.z > p.boxMin[0].z + 1e-3f && q.z < p.boxMax[0].z - 1e-3f) inside++;
        }
        CHECK(finite);
        CHECK(inside == 0);                                             // the box is solid — nothing penetrates it
        sph.free();
    }

    // E) recirculating emitter: particles that sink past recycleY respawn in the emit AABB,
    //    so the count is conserved and none linger below the drain plane
    {
        Params p = baseParams(0.05f);
        p.gmin = make_float3(-0.4f, -0.5f, -0.4f); p.gmax = make_float3(0.4f, 1.1f, 0.4f);
        p.emitOn = 1; p.recycleY = -0.20f;
        p.emitMin = make_float3(-0.30f, 0.85f, -0.30f); p.emitMax = make_float3(0.30f, 1.00f, 0.30f);
        p.emitVel = make_float3(0.f, 0.f, 0.f);
        auto pos = lattice(make_float3(-0.30f, 0.05f, -0.30f), make_float3(0.30f, 0.45f, 0.30f), 0.045f);
        int n0 = (int)pos.size();
        GpuSPH sph; sph.init(pos, p);                                   // no floor → water falls, drains, respawns up top
        for (int s = 0; s < 500; s++) sph.step(1.5e-3f);
        std::vector<float3> out; sph.download(out);
        CHECK((int)out.size() == n0);                                   // fixed-N recirculation conserves count
        int below = 0, inEmit = 0; bool finite = true;
        for (auto& q : out) {
            if (!std::isfinite(q.y)) finite = false;
            if (q.y < p.recycleY - 1e-2f) below++;
            if (q.y > p.emitMin.y - 0.05f && q.x > p.emitMin.x - 0.05f && q.x < p.emitMax.x + 0.05f) inEmit++;
        }
        CHECK(finite);
        CHECK(below == 0);                                              // drained particles are teleported, never left below
        CHECK(inEmit > 0);                                              // some have respawned at the source
        sph.free();
    }

    // F) cohesion (surface tension): a free blob in zero-g holds together more tightly with
    //    surfTens>0 — its radius of gyration is smaller than the same run without cohesion.
    {
        Params p = baseParams(0.10f);
        p.gmin = make_float3(-2, -2, -2); p.gmax = make_float3(2, 2, 2);
        p.grav = make_float3(0, 0, 0); p.visc = 1.f;
        auto blob = lattice(make_float3(-0.12f, -0.12f, -0.12f), make_float3(0.12f, 0.12f, 0.12f), 0.04f);
        auto radiusGyration = [&](float st) {
            Params pp = p; pp.surfTens = st; GpuSPH s; s.init(blob, pp);
            for (int k = 0; k < 90; k++) s.step(8.0e-4f);
            std::vector<float3> o; s.download(o); s.free();
            float3 c = make_float3(0, 0, 0); for (auto& q : o) c = c + q; c = c * (1.0f / o.size());
            double rg = 0; for (auto& q : o) { float3 d = q - c; rg += dot3(d, d); } return std::sqrt(rg / o.size());
        };
        double rgNo = radiusGyration(0.f), rgCoh = radiusGyration(400.f);
        CHECK(rgNo > 0);
        CHECK(rgCoh < rgNo - 0.003);                                    // cohesion keeps the blob more compact
    }

    // G) Coriolis (rotating frame): a lone particle with no other force conserves speed
    //    and its velocity vector rotates at rate 2·|omega| (an inertial oscillation),
    //    while with omega=0 it flies straight.
    {
        Params p = baseParams(0.10f);
        p.gmin = make_float3(-3, -3, -3); p.gmax = make_float3(3, 3, 3);
        p.grav = make_float3(0, 0, 0);
        const float Om = 3.0f; const float dt = 1.0e-3f; const int N = 200;
        std::vector<float3> pos = {make_float3(0, 0, 0)};
        auto velAfter = [&](float3 om) {
            Params pp = p; pp.omega = om; GpuSPH s; s.init(pos, pp); s.setVel({make_float3(0.5f, 0, 0)});
            for (int k = 0; k < N; k++) s.step(dt);
            std::vector<float3> v; s.downloadVel(v); s.free(); return v[0];
        };
        float3 vR = velAfter(make_float3(0, Om, 0));                     // rotating frame
        float speed = std::sqrt(vR.x * vR.x + vR.y * vR.y + vR.z * vR.z);
        float ang = std::atan2(vR.z, vR.x);
        CHECK_NEAR(speed, 0.5, 0.03);                                    // Coriolis does no work → |v| conserved
        CHECK_NEAR(ang, 2.0 * Om * (N * dt), 0.12);                      // velocity rotated by 2·omega·T
        float3 vI = velAfter(make_float3(0, 0, 0));                      // inertial frame
        CHECK_NEAR(vI.x, 0.5, 0.02); CHECK(std::fabs(vI.z) < 0.02);      // flies straight, no deflection
    }

    return test::report("gpu");
}
