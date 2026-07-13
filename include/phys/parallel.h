// Island-parallel solver. A rigid-body world partitions into *islands* — groups of
// bodies/contacts that can only interact within the group — so each island's solve
// is independent and can run on its own thread (ODE/PhysX/Bullet all do this). We
// farm N islands across a small fixed thread pool and join. Because every worker
// writes ONLY into its own island's slot, the output is bit-for-bit identical to
// solving the islands serially in order: threading changes the schedule, never the
// arithmetic. Compile with -pthread.
#pragma once
#include "precision.h"
#include <thread>
#include <vector>
#include <atomic>

namespace phys {

// Solve every island in order on the calling thread — the deterministic reference.
template<class Island, class Solve>
inline void solveIslandsSerial(std::vector<Island>& islands, Solve solve) {
    for (unsigned i = 0; i < islands.size(); i++) solve(i, islands[i]);
}

// Run `solve(i, islands[i])` for every island across a pool of worker threads and
// join. Workers pull islands off a shared atomic counter (dynamic load balance);
// since each writes only island i's own slot, completion order is irrelevant and
// the result matches solveIslandsSerial byte-for-byte. `solve` MUST be independent
// per island (no shared writes) for this guarantee to hold. `maxThreads == 0`
// means "use hardware_concurrency". Returns the number of worker threads used.
template<class Island, class Solve>
inline unsigned solveIslands(std::vector<Island>& islands, Solve solve,
                             unsigned maxThreads = 0) {
    const unsigned n = (unsigned)islands.size();
    if (n == 0) return 0;

    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) hw = 1;                                  // unknown → assume single core
    unsigned nthreads = maxThreads ? maxThreads : hw;
    if (nthreads > n) nthreads = n;                       // no more threads than work
    if (nthreads < 1) nthreads = 1;

    if (nthreads == 1) {                                  // serial fast path (also deterministic)
        solveIslandsSerial(islands, solve);
        return 1;
    }

    std::atomic<unsigned> next{0};
    auto worker = [&]() {
        for (;;) {
            unsigned i = next.fetch_add(1, std::memory_order_relaxed);
            if (i >= n) break;
            solve(i, islands[i]);                         // writes island i only
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(nthreads - 1);
    for (unsigned t = 0; t + 1 < nthreads; t++) pool.emplace_back(worker);
    worker();                                             // the calling thread is a worker too
    for (auto& th : pool) th.join();
    return nthreads;
}

} // namespace phys
