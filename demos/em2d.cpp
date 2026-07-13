// Maxwell + Method of Moments, visualized. A center-fed dipole antenna is solved
// in the frequency domain by the Method of Moments (phys::em::Dipole) for its
// current distribution and radiation pattern; that current then drives a 2-D FDTD
// time-domain solution of Maxwell's equations (phys::em::FdtdTE), so we watch the
// electric and magnetic fields actually radiate away from the antenna.
//
//   left  : magnetic field  H_z  (diverging colour, the expanding wavefronts) with
//           the in-plane electric field  E = (E_x,E_y)  drawn as flow arrows,
//   right : the MoM radiation pattern (polar) and current distribution |I(z)|.
//
//   ./em2d --video out/f 300   |   ./em2d --shot frame.png 220
#include "phys/fdtd.h"
#include "phys/mom.h"
#include "png.h"
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <algorithm>
using namespace phys::em;

static const int W = 1200, H = 760;
static std::vector<unsigned char> img(W * H * 3);

static inline void blend(int x, int y, double r, double g, double b, double a) {
    if (x < 0 || y < 0 || x >= W || y >= H || a <= 0) return; if (a > 1) a = 1;
    unsigned char* p = &img[(y * W + x) * 3];
    p[0] = (unsigned char)std::min(255.0, p[0] * (1 - a) + r * a);
    p[1] = (unsigned char)std::min(255.0, p[1] * (1 - a) + g * a);
    p[2] = (unsigned char)std::min(255.0, p[2] * (1 - a) + b * a);
}
static void line(double x0, double y0, double x1, double y1, double r, double g, double b, double a, double w = 1.0) {
    double dx = x1 - x0, dy = y1 - y0; int n = (int)(std::max(std::fabs(dx), std::fabs(dy)) + 1);
    for (int i = 0; i <= n; i++) { double t = (double)i / n, cx = x0 + dx * t, cy = y0 + dy * t;
        for (double ox = -w; ox <= w; ox += 1.0) for (double oy = -w; oy <= w; oy += 1.0) {
            double fall = 1.0 - std::sqrt(ox * ox + oy * oy) / (w + 0.5); if (fall > 0) blend((int)(cx + ox), (int)(cy + oy), r, g, b, a * fall); } }
}

// diverging colormap for the (signed) magnetic field on a dark ground. Brightness
// tracks |field| (a gentle power, not sqrt) so the amplitude decay and the dipole's
// directivity read as light/dark instead of being flattened by saturation.
static void hcolor(double s, double& r, double& g, double& b) {
    s = std::max(-1.0, std::min(1.0, s));
    double a = std::fabs(s), m = std::pow(a, 0.72);
    double R = 8, G = 9, B = 14;                              // base
    if (s >= 0) { R += m * 245; G += m * 105 + m * m * 130; B += m * 25; }   // warm: ember→gold→white
    else        { R += m * 30 + m * m * 60; G += m * 120; B += m * 240; }    // cool: indigo→sky
    r = R; g = G; b = B;
}

int main(int argc, char** argv) {
    const char* video = nullptr; const char* shot = nullptr; int frames = 300;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--video")) { video = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
        else if (!strcmp(argv[i], "--shot")) { shot = argv[++i]; if (i + 1 < argc) frames = atoi(argv[i + 1]); }
    }

    // --- Method of Moments: solve the dipole current & pattern once ---
    Dipole dip; dip.L = 0.5; dip.a = 0.005; dip.N = 41; dip.solve();
    double Imax = dip.maxCurrent();
    double patMax = 0; for (int t = 0; t <= 180; t++) patMax = std::max(patMax, dip.pattern(t * M_PI / 180));

    // --- FDTD grid & antenna geometry ---
    const int nx = 380, ny = 380, bd = 34;
    const double lam = 44.0;                                  // wavelength in cells
    const double omega = 2 * M_PI / lam;                     // c = 1
    FdtdTE fd(nx, ny, 0.5); fd.border = bd;
    const int ic = nx / 2;                                    // antenna column
    const double half = 0.5 * dip.L * lam;                   // half-length in cells
    const int jlo = (int)std::round(ny / 2 - half), jhi = (int)std::round(ny / 2 + half);
    // MoM current sampled onto each antenna cell (magnitude + phase), for a soft-started source
    auto Jat = [&](int i, int j, double time) -> double {
        if (i != ic || j < jlo || j > jhi) return 0.0;
        double frac = (double)(j - jlo) / (jhi - jlo);        // 0..1 along the wire
        double zn = (frac - 0.5) * dip.L;                     // wavelengths from center
        // nearest MoM segment
        int seg = std::min(dip.N - 1, std::max(0, (int)std::round((zn + 0.5 * dip.L) / (dip.L / dip.N) - 0.5)));
        std::complex<double> Iph = dip.I[seg];
        double env = 1.0 - std::exp(-time / (3 * lam));       // smooth ramp-up
        return env * (std::abs(Iph) / Imax) * std::cos(omega * time + std::arg(Iph));
    };

    // --- field-panel placement: FDTD interior (minus absorbing border) → left square ---
    const int PS = H;                                         // panel size (square, full height)
    const int gi0 = bd, gi1 = nx - bd, gj0 = bd, gj1 = ny - bd;
    auto sampleHz = [&](double gi, double gj) { int i = (int)gi, j = (int)gj; i = std::max(0, std::min(nx - 2, i)); j = std::max(0, std::min(ny - 2, j));
        double fi = gi - i, fj = gj - j; auto v = [&](int a, int b) { return fd.Hz[fd.id(a, b)]; };
        return v(i, j) * (1 - fi) * (1 - fj) + v(i + 1, j) * fi * (1 - fj) + v(i, j + 1) * (1 - fi) * fj + v(i + 1, j + 1) * fi * fj; };

    auto renderFrame = [&]() {
        // background
        for (int k = 0; k < W * H; k++) { img[k * 3] = 8; img[k * 3 + 1] = 9; img[k * 3 + 2] = 14; }
        // running field scale so the near field saturates and the radiation shows
        double vmax = 1e-9; for (int i = gi0; i < gi1; i++) for (int j = gj0; j < gj1; j++) vmax = std::max(vmax, std::fabs(fd.Hz[fd.id(i, j)]));
        double scale = 0.5 * vmax + 1e-6;                    // higher scale ⇒ far field dims, near field clips (shows decay)
        // Hz colour field
        for (int py = 0; py < PS; py++) for (int px = 0; px < PS; px++) {
            double gi = gi0 + (double)px / PS * (gi1 - gi0);
            double gj = gj0 + (double)py / PS * (gj1 - gj0);
            double r, g, b; hcolor(sampleHz(gi, gj) / scale, r, g, b);
            unsigned char* p = &img[((py) * W + px) * 3]; p[0] = (unsigned char)std::min(255.0, r); p[1] = (unsigned char)std::min(255.0, g); p[2] = (unsigned char)std::min(255.0, b);
        }
        // E-field flow arrows on a regular lattice
        double emax = 1e-9; for (int i = gi0; i < gi1; i++) for (int j = gj0; j < gj1; j++) { double e = std::hypot(fd.Ex[fd.id(i, j)], fd.Ey[fd.id(i, j)]); emax = std::max(emax, e); }
        int stepc = 20; double cellpx = (double)PS / (gi1 - gi0);
        for (int gj = gj0 + stepc / 2; gj < gj1; gj += stepc) for (int gi = gi0 + stepc / 2; gi < gi1; gi += stepc) {
            double ex = fd.Ex[fd.id(gi, gj)], ey = fd.Ey[fd.id(gi, gj)]; double mag = std::hypot(ex, ey); if (mag < 1e-9) continue;
            double s = std::tanh(3.0 * mag / emax); if (s < 0.05) continue;
            double px = (gi - gi0) * cellpx, py = (gj - gj0) * cellpx;
            double L = stepc * cellpx * 0.42 * s, ux = ex / mag, uy = ey / mag;
            double x1 = px + ux * L, y1 = py + uy * L, x0 = px - ux * L, y0 = py - uy * L;
            double a = 0.25 + 0.55 * s;
            line(x0, y0, x1, y1, 230, 240, 255, a, 0.7);
            line(x1, y1, x1 - (ux * 0.5 + uy * 0.5) * L, y1 - (uy * 0.5 - ux * 0.5) * L, 230, 240, 255, a, 0.7);   // arrowhead
            line(x1, y1, x1 - (ux * 0.5 - uy * 0.5) * L, y1 - (uy * 0.5 + ux * 0.5) * L, 230, 240, 255, a, 0.7);
        }
        // glowing antenna rod
        double ax = (ic - gi0) * cellpx;
        for (int j = jlo; j <= jhi; j++) { double ay = (j - gj0) * cellpx;
            for (int gx = -6; gx <= 6; gx++) { double fall = std::exp(-gx * gx / 6.0); blend((int)(ax + gx), (int)ay, 255, 235, 150, 0.9 * fall); } }
        line(ax, (jlo - gj0) * cellpx, ax, (jhi - gj0) * cellpx, 255, 245, 200, 0.95, 1.4);

        // ---- right panels ----
        int rx = PS + 24, rw = W - rx - 24;
        // radiation pattern (polar), broadside pointing +x
        int cxp = rx + rw / 2, cyp = 190, Rp = 150;
        for (int k = 0; k < 4; k++) { double rr = Rp * (k + 1) / 4.0;                 // range rings
            for (int t = 0; t < 360; t += 3) line(cxp + rr * cos(t * M_PI / 180), cyp + rr * sin(t * M_PI / 180), cxp + rr * cos((t + 3) * M_PI / 180), cyp + rr * sin((t + 3) * M_PI / 180), 60, 66, 80, 0.5, 0.5); }
        line(cxp - Rp, cyp, cxp + Rp, cyp, 70, 76, 92, 0.6, 0.5); line(cxp, cyp - Rp, cxp, cyp + Rp, 70, 76, 92, 0.6, 0.5);
        double pxprev = 0, pyprev = 0; bool have = false;
        for (int a = 0; a <= 360; a += 2) {
            double th = a * M_PI / 180;                        // pattern angle from the wire axis (vertical)
            double val = dip.pattern(th) / patMax;             // radius
            // place: axis is vertical (antenna along y); broadside (θ=90°) → horizontal lobes
            double vx = cxp + Rp * val * std::sin(th), vy = cyp - Rp * val * std::cos(th);
            if (have) line(pxprev, pyprev, vx, vy, 255, 110, 60, 0.95, 1.3);
            pxprev = vx; pyprev = vy; have = true;
        }
        // current distribution |I(z)| along a vertical axis matching the antenna
        int cx0 = rx + 40, cy0 = 430, cyh = 260, cw = rw - 80;
        line(cx0, cy0, cx0, cy0 + cyh, 70, 76, 92, 0.6, 0.6);                          // axis
        line(cx0, cy0 + cyh / 2, cx0 + cw, cy0 + cyh / 2, 50, 56, 70, 0.5, 0.5);
        double ppx = 0, ppy = 0; bool h2 = false;
        for (int n = 0; n < dip.N; n++) {
            double yy = cy0 + (double)n / (dip.N - 1) * cyh;
            double xx = cx0 + std::abs(dip.I[n]) / Imax * cw;
            if (h2) line(ppx, ppy, xx, yy, 120, 220, 255, 0.95, 1.3);
            ppx = xx; ppy = yy; h2 = true;
        }
        line(cx0, cy0, cx0 + cw, cy0, 40, 44, 56, 0.6, 0.4); line(cx0, cy0 + cyh, cx0 + cw, cy0 + cyh, 40, 44, 56, 0.6, 0.4);
    };

    auto advance = [&](int steps) { for (int s = 0; s < steps; s++) fd.step([&](int i, int j) { return Jat(i, j, fd.t); }); };

    if (video) {
        for (int f = 0; f < frames; f++) { advance(4); renderFrame(); char q[512]; std::snprintf(q, sizeof(q), "%s_%04d.png", video, f); gfx::writePNG(q, img.data(), W, H, false); }
        std::printf("wrote %d frames (%d MoM segments, Zin=%.1f%+.1fj)\n", frames, dip.N, dip.Zin.real(), dip.Zin.imag());
        return 0;
    }
    int target = shot ? frames : 220;
    for (int f = 0; f < target; f++) advance(4);
    renderFrame(); gfx::writePNG(shot ? shot : "em2d.png", img.data(), W, H, false);
    std::printf("wrote shot (Zin=%.1f%+.1fj)\n", dip.Zin.real(), dip.Zin.imag());
    return 0;
}
