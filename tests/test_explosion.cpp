// Explosion / fracture: grid fracture conserves the block, and detonate() throws
// every fragment outward with an upward plume, scales its energy with the charge,
// and leaves fragments outside the blast radius untouched.
#include "phys/fracture.h"
#include "phys/voronoi.h"
#include "phys/body.h"
#include "check.h"
#include <vector>
#include <memory>
#include <cmath>
using namespace phys;

// build fragment bodies for a block and return the Destructible + kept bodies
static void makeBlock(const Vector3& c, const Vector3& half, int n, real jit, unsigned seed,
                      std::vector<std::unique_ptr<RigidBody>>& own, Destructible& dest) {
    auto descs = fractureBoxGrid(c, half, n, n, n, jit, seed);
    for (auto& d : descs) {
        auto b = std::make_unique<RigidBody>();
        real vol = 8 * d.halfSize.x * d.halfSize.y * d.halfSize.z;
        b->setMass(vol * 2.4);                              // concrete density
        Matrix3 I; I.setBlockInertiaTensor(d.halfSize, b->getMass()); b->setInertiaTensor(I);
        b->setPosition(d.centre); b->calculateDerivedData(); b->setAwake(false);
        dest.fragments.push_back(b.get()); own.push_back(std::move(b));
    }
}

int main() {
    const Vector3 C(0, 0, 0), HALF(1.0, 1.0, 1.0);

    // A) grid fracture conserves total volume and fragment count
    {
        auto descs = fractureBoxGrid(C, HALF, 4, 4, 4, 0.4, 7u);
        CHECK(descs.size() == 64);
        real vol = 0; for (auto& d : descs) vol += 8 * d.halfSize.x * d.halfSize.y * d.halfSize.z;
        CHECK_NEAR(vol, 8 * HALF.x * HALF.y * HALF.z, 1e-9);          // = block volume, exactly
    }

    // B) detonate throws fragments outward and adds a net upward plume momentum
    {
        std::vector<std::unique_ptr<RigidBody>> own; Destructible dest;
        makeBlock(C, HALF, 5, 0.0, 7u, own, dest);                   // symmetric grid (jitter 0)
        dest.detonate(C, 6.0, 0.4, 10.0, 5.0, 11u);
        CHECK(dest.shattered);
        // mean outward radial speed is clearly positive
        real outward = 0; int cnt = 0; Vector3 P(0, 0, 0); real M = 0;
        for (auto* f : dest.fragments) {
            Vector3 d = f->getPosition() - C; real dist = d.magnitude();
            if (dist > 1e-6) outward += (f->getVelocity() * (d * (1.0 / dist)));
            cnt++; P += f->getVelocity() * f->getMass(); M += f->getMass();
        }
        CHECK(outward / cnt > 1.0);                                  // fragments fly apart
        // radial x/z momentum cancels by symmetry; the plume leaves a net +y
        CHECK(P.y > 0);
        CHECK(std::fabs(P.x) < 0.5 * P.y && std::fabs(P.z) < 0.5 * P.y);
        int awake = 0; for (auto* f : dest.fragments) awake += f->getAwake() ? 1 : 0;
        CHECK(awake == (int)dest.fragments.size());                  // all woken
    }

    // C) blast energy scales with the square of the charge strength
    {
        auto ke = [&](real strength) {
            std::vector<std::unique_ptr<RigidBody>> own; Destructible dest;
            makeBlock(C, HALF, 4, 0.3, 7u, own, dest);
            dest.detonate(C, strength, 0.4, 10.0, 5.0, 11u);
            real e = 0; for (auto* f : dest.fragments) e += f->getKineticEnergy(); return e;
        };
        real e1 = ke(3.0), e2 = ke(6.0);
        CHECK(e1 > 0);
        CHECK_NEAR(e2 / e1, 4.0, 0.4);                               // 2× strength → ~4× energy
    }

    // D) fragments outside the blast radius are left undisturbed
    {
        std::vector<std::unique_ptr<RigidBody>> own; Destructible dest;
        makeBlock(Vector3(0, 0, 0), Vector3(2.0, 2.0, 2.0), 5, 0.0, 7u, own, dest);
        dest.detonate(Vector3(0, 0, 0), 6.0, 0.4, 1.0, 5.0, 11u);    // small radius
        int thrown = 0, still = 0, stillAsleep = 0;
        for (auto* f : dest.fragments) {
            if (f->getVelocity().magnitude() > 1e-9) thrown++;
            else { still++; if (!f->getAwake()) stillAsleep++; }
        }
        CHECK(thrown > 0 && still > 0);                              // near ones fly, far ones stay
        CHECK(stillAsleep == still);                                 // untouched fragments remain asleep
    }

    // E) Voronoi fracture: cells partition the block (Σ volume ≈ box), are convex,
    //    and every cell has valid (positive-definite-ish) mass properties
    {
        Vector3 half(1.0, 1.2, 0.9); real boxV = 8 * half.x * half.y * half.z;
        unsigned st = 999u; auto rr = [&]() { st = st * 1103515245u + 12345u; return (real)(((st >> 16) & 0x7fff) / 32767.0); };
        std::vector<Vector3> seeds;
        for (int i = 0; i < 50; i++) seeds.push_back(Vector3((rr() * 2 - 1) * half.x * 0.9, (rr() * 2 - 1) * half.y * 0.9, (rr() * 2 - 1) * half.z * 0.9));
        auto cells = voronoiFracture(half, seeds, 20);
        CHECK(cells.size() == seeds.size());                            // every seed yields a cell
        real vol = 0; int badI = 0, notConvex = 0;
        for (auto& c : cells) {
            vol += c.volume;
            CHECK(c.volume > 0);
            if (!(c.inertiaUnit.data[0] > 0 && c.inertiaUnit.data[4] > 0 && c.inertiaUnit.data[8] > 0)) badI++;
            // convex: every vertex is inside (or on) every face's supporting plane
            for (size_t f = 0; f < c.tris.size(); f++) { Vector3 n = c.triN[f]; real d = n * c.verts[c.tris[f][0]];
                for (auto& v : c.verts) if (n * v > d + 1e-3) { notConvex++; break; } }
        }
        CHECK_NEAR(vol / boxV, 1.0, 0.02);                              // cells tile the block (≈ exact partition)
        CHECK(badI == 0);
        CHECK(notConvex == 0);                                          // all fragments are convex
    }

    return test::report("explosion");
}
