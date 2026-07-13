// Method of Moments (MoM) — the frequency-domain integral-equation workhorse of
// antenna engineering, here for a center-fed thin-wire dipole. We solve the
// electric-field integral equation (Pocklington's form) for the current
// distribution along the wire, then read off the input impedance and the
// far-field radiation pattern.
//
//   Pocklington:  ∫ I(z') [ ∂²/∂z² + k² ] G(z,z') dz' = -jωε · E_z^inc(z)
//   with the thin-wire kernel  G(z,z') = e^{-jkR}/(4πR),  R = √((z−z')² + a²)
//
// Discretized with pulse (rectangular) basis functions and point matching, the
// ∂²/∂z² becomes a central difference across neighbouring match points, giving a
// dense complex linear system Z·I = V that we solve directly. A delta-gap source
// models the feed. Everything is in wavelength units (λ = 1), so k = 2π and the
// free-space impedance η ≈ 376.73 Ω sets the ohmic scale.
#pragma once
#include <vector>
#include <complex>
#include <cmath>

namespace phys { namespace em {

using cd = std::complex<double>;

struct Dipole {
    double L = 0.5;      // total length   (wavelengths)
    double a = 0.005;    // wire radius    (wavelengths)
    int    N = 41;       // segments (odd → a segment centred on the feed)
    // --- results (filled by solve(), for a 1 V delta-gap excitation) ---
    std::vector<double> z;   // segment-centre positions
    std::vector<cd>     I;   // current per segment  (A, phasor)
    cd Zin;                  // input impedance V/I(feed)  (Ω)

    static constexpr double PI  = 3.14159265358979323846;
    static constexpr double ETA = 376.730313668;             // √(μ0/ε0)

    // segment-integrated thin-wire Green's function ∫_seg e^{-jkR}/(4πR) dz',
    // R = √((zf−z')²+a²). We split off the 1/(4πR) singularity (integrated
    // analytically via asinh) and integrate the smooth, bounded remainder
    // (e^{-jkR}−1)/(4πR) numerically — accurate for the self cell and all others,
    // which the radiation resistance is sensitive to.
    static cd psi(double zf, double zn, double d, double a, double k, int Q = 16) {
        double off = zf - zn;
        double re = (std::asinh((0.5 * d - off) / a) + std::asinh((0.5 * d + off) / a)) / (4.0 * PI);
        cd sm(0, 0);
        for (int q = 0; q < Q; q++) {
            double u = -0.5 * d + (q + 0.5) * d / Q;         // segment-local coordinate (midpoint rule)
            double R = std::sqrt((off - u) * (off - u) + a * a);
            sm += (std::exp(cd(0.0, -k * R)) - 1.0) / (4.0 * PI * R) * (d / Q);
        }
        return cd(re, 0.0) + sm;
    }

    // solve Z·I = V by Gaussian elimination with partial pivoting (complex)
    static std::vector<cd> solveLin(std::vector<std::vector<cd>> A, std::vector<cd> b) {
        int n = (int)b.size();
        for (int col = 0; col < n; col++) {
            int piv = col; double best = std::abs(A[col][col]);
            for (int r = col + 1; r < n; r++) { double m = std::abs(A[r][col]); if (m > best) { best = m; piv = r; } }
            std::swap(A[col], A[piv]); std::swap(b[col], b[piv]);
            cd d = A[col][col];
            for (int r = col + 1; r < n; r++) {
                cd f = A[r][col] / d;
                for (int c = col; c < n; c++) A[r][c] -= f * A[col][c];
                b[r] -= f * b[col];
            }
        }
        std::vector<cd> x(n);
        for (int r = n - 1; r >= 0; r--) {
            cd s = b[r];
            for (int c = r + 1; c < n; c++) s -= A[r][c] * x[c];
            x[r] = s / A[r][r];
        }
        return x;
    }

    void solve(double V = 1.0) {
        const double k = 2.0 * PI, d = L / N;
        z.resize(N); for (int i = 0; i < N; i++) z[i] = -0.5 * L + (i + 0.5) * d;
        // field points span [-1 .. N] so the central difference is defined at every match point
        auto zf = [&](int m) { return -0.5 * L + (m + 0.5) * d; };
        std::vector<std::vector<cd>> Z(N, std::vector<cd>(N));
        for (int m = 0; m < N; m++)
            for (int n = 0; n < N; n++) {
                cd g0 = psi(zf(m),     z[n], d, a, k);
                cd gp = psi(zf(m + 1), z[n], d, a, k);
                cd gm = psi(zf(m - 1), z[n], d, a, k);
                Z[m][n] = (gp - 2.0 * g0 + gm) / (d * d) + k * k * g0;   // (∂²/∂z² + k²) G
            }
        // delta-gap: E_z^inc = V/d over the feed segment only;  RHS = -jωε·E_inc, ωε = k/η
        std::vector<cd> b(N, cd(0, 0));
        int feed = N / 2;
        b[feed] = cd(0.0, -k / ETA) * (V / d);
        I = solveLin(Z, b);
        Zin = V / I[feed];
    }

    // far-field pattern factor F(θ):  E_θ ∝ sinθ · ∫ I(z') e^{jk z' cosθ} dz'
    cd farField(double theta) const {
        const double k = 2.0 * PI, d = L / N;
        cd s(0, 0);
        for (int n = 0; n < N; n++) s += I[n] * std::exp(cd(0.0, k * z[n] * std::cos(theta)));
        return std::sin(theta) * s * d;
    }
    double pattern(double theta) const { return std::abs(farField(theta)); }

    // peak |I| and the total radiated-power proxy — small helpers for tests/plots
    double maxCurrent() const { double m = 0; for (auto& c : I) m = std::max(m, std::abs(c)); return m; }
};

}} // namespace phys::em
