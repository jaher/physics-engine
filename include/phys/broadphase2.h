// Alternative broad phases to the core bounding-sphere BVH (collide_coarse.h):
//   1. Sweep-and-Prune — sort AABB endpoints on one axis, sweep for overlaps
//      (the classic I-Collide / Bullet btAxisSweep approach).
//   2. Dynamic AABB tree — an incremental, self-balancing tree of fat AABBs with
//      insert / remove / rebalance, queryAABB and queryPairs, modelled on Bullet's
//      btDbvt and Box2D's b2DynamicTree.
// Both return the set of candidate overlapping pairs; the narrow phase confirms.
#pragma once
#include "query.h"          // Aabb
#include <vector>
#include <algorithm>
#include <numeric>
#include <utility>

namespace phys {

using BroadPair = std::pair<int, int>;   // (idA, idB) with idA < idB

inline BroadPair makePair(int a, int b) { return a < b ? BroadPair(a, b) : BroadPair(b, a); }

// ============================================================ Sweep-and-Prune
// Insert AABBs (each with a caller id, default = insertion index); overlappingPairs
// sorts the box lower-x endpoints and sweeps a single axis, doing the full 3-D test
// only on the x-overlapping neighbourhood.
class SweepAndPrune {
public:
    struct Proxy { Aabb box; int id; };
    std::vector<Proxy> proxies;

    int add(const Aabb& box, int id = -1) {
        int idx = (int)proxies.size();
        proxies.push_back({box, id < 0 ? idx : id});
        return idx;
    }
    void update(int idx, const Aabb& box) { proxies[(size_t)idx].box = box; }
    void clear() { proxies.clear(); }
    size_t size() const { return proxies.size(); }

    std::vector<BroadPair> overlappingPairs() const {
        std::vector<BroadPair> pairs;
        const size_t n = proxies.size();
        if (n < 2) return pairs;
        std::vector<int> order(n);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return proxies[(size_t)a].box.min.x < proxies[(size_t)b].box.min.x; });
        for (size_t a = 0; a < n; a++) {
            const Proxy& pa = proxies[(size_t)order[a]];
            for (size_t b = a + 1; b < n; b++) {
                const Proxy& pb = proxies[(size_t)order[b]];
                if (pb.box.min.x > pa.box.max.x) break;      // sorted: no later box can overlap in x
                if (pa.box.overlaps(pb.box)) pairs.push_back(makePair(pa.id, pb.id));
            }
        }
        return pairs;
    }
};

// ============================================================ Dynamic AABB tree
// Node pool with a free list (index-stable, header-only). Internal nodes hold the
// union of their children's fat AABBs; leaves also keep the tight AABB so queryPairs
// reproduces the exact broad-phase overlap set regardless of the fattening margin.
class DynamicAabbTree {
public:
    explicit DynamicAabbTree(real margin = 0) : margin_(margin) {}

    // Insert a leaf for `id` with tight bounds `box`; returns the proxy handle.
    int insert(int id, const Aabb& box) {
        int leaf = allocateNode();
        nodes_[leaf].tight  = box;
        nodes_[leaf].aabb   = box.expanded(margin_);
        nodes_[leaf].id     = id;
        nodes_[leaf].height = 0;
        nodes_[leaf].child1 = nodes_[leaf].child2 = NIL;
        insertLeaf(leaf);
        return leaf;
    }
    void remove(int proxy) { removeLeaf(proxy); freeNode(proxy); }

    // Reinsert with new bounds (remove+insert keeps the tree balanced) — the
    // incremental "rebalance" a moving body triggers each frame.
    int update(int proxy, const Aabb& box) {
        int id = nodes_[proxy].id;
        remove(proxy);
        return insert(id, box);
    }

    // All leaf ids whose fat AABB overlaps `q`.
    std::vector<int> queryAABB(const Aabb& q) const {
        std::vector<int> out;
        if (root_ == NIL) return out;
        std::vector<int> stack; stack.push_back(root_);
        while (!stack.empty()) {
            int n = stack.back(); stack.pop_back();
            if (n == NIL || !nodes_[n].aabb.overlaps(q)) continue;
            if (nodes_[n].isLeaf()) out.push_back(nodes_[n].id);
            else { stack.push_back(nodes_[n].child1); stack.push_back(nodes_[n].child2); }
        }
        return out;
    }

    // All overlapping leaf-id pairs (tight-AABB test → exact broad-phase set).
    std::vector<BroadPair> queryPairs() const {
        std::vector<BroadPair> pairs;
        if (root_ == NIL || nodes_[root_].isLeaf()) return pairs;
        // gather leaves, then descend from each with its fat AABB
        std::vector<int> leaves; collectLeaves(root_, leaves);
        for (int leaf : leaves) {
            std::vector<int> stack; stack.push_back(root_);
            const Aabb& fat = nodes_[leaf].aabb;
            const Aabb& tight = nodes_[leaf].tight;
            int selfId = nodes_[leaf].id;
            while (!stack.empty()) {
                int n = stack.back(); stack.pop_back();
                if (n == NIL || !nodes_[n].aabb.overlaps(fat)) continue;
                if (nodes_[n].isLeaf()) {
                    if (nodes_[n].id > selfId && nodes_[n].tight.overlaps(tight))
                        pairs.push_back(makePair(selfId, nodes_[n].id));
                } else { stack.push_back(nodes_[n].child1); stack.push_back(nodes_[n].child2); }
            }
        }
        return pairs;
    }

    int  height()   const { return root_ == NIL ? 0 : nodes_[root_].height; }
    bool empty()    const { return root_ == NIL; }
    int  leafCount() const { std::vector<int> l; if (root_ != NIL) collectLeaves(root_, l); return (int)l.size(); }

private:
    static const int NIL = -1;
    struct Node {
        Aabb aabb;             // fat bounds (this node & subtree)
        Aabb tight;            // tight bounds (leaves)
        int  id     = -1;      // user id (leaves)
        int  parent = NIL;
        int  child1 = NIL, child2 = NIL;
        int  height = -1;      // -1 ⇒ free
        int  next   = NIL;     // free-list link
        bool isLeaf() const { return child1 == NIL; }
    };
    std::vector<Node> nodes_;
    int   root_     = NIL;
    int   freeList_ = NIL;
    real  margin_;

    int allocateNode() {
        if (freeList_ == NIL) {
            nodes_.push_back(Node());
            return (int)nodes_.size() - 1;
        }
        int n = freeList_;
        freeList_ = nodes_[n].next;
        nodes_[n] = Node();
        return n;
    }
    void freeNode(int n) { nodes_[n].height = -1; nodes_[n].next = freeList_; freeList_ = n; }

    void collectLeaves(int n, std::vector<int>& out) const {
        if (n == NIL) return;
        if (nodes_[n].isLeaf()) out.push_back(n);
        else { collectLeaves(nodes_[n].child1, out); collectLeaves(nodes_[n].child2, out); }
    }

    void insertLeaf(int leaf) {
        if (root_ == NIL) { root_ = leaf; nodes_[leaf].parent = NIL; return; }
        // 1. find the best sibling by surface-area heuristic
        Aabb leafBox = nodes_[leaf].aabb;
        int index = root_;
        while (!nodes_[index].isLeaf()) {
            int c1 = nodes_[index].child1, c2 = nodes_[index].child2;
            real area = nodes_[index].aabb.area();
            Aabb combined = nodes_[index].aabb.merged(leafBox);
            real combinedArea = combined.area();
            real cost = 2 * combinedArea;                       // cost of a new parent here
            real inheritance = 2 * (combinedArea - area);       // extra cost pushed to descendants
            auto descendCost = [&](int child) {
                Aabb merged = leafBox.merged(nodes_[child].aabb);
                if (nodes_[child].isLeaf()) return merged.area() + inheritance;
                return (merged.area() - nodes_[child].aabb.area()) + inheritance;
            };
            real cost1 = descendCost(c1), cost2 = descendCost(c2);
            if (cost < cost1 && cost < cost2) break;            // best to split here
            index = (cost1 < cost2) ? c1 : c2;
        }
        int sibling = index;
        // 2. create a new parent for sibling and leaf
        int oldParent = nodes_[sibling].parent;
        int newParent = allocateNode();
        nodes_[newParent].parent = oldParent;
        nodes_[newParent].aabb   = leafBox.merged(nodes_[sibling].aabb);
        nodes_[newParent].height = nodes_[sibling].height + 1;
        nodes_[newParent].child1 = sibling;
        nodes_[newParent].child2 = leaf;
        nodes_[sibling].parent = newParent;
        nodes_[leaf].parent    = newParent;
        if (oldParent != NIL) {
            if (nodes_[oldParent].child1 == sibling) nodes_[oldParent].child1 = newParent;
            else                                     nodes_[oldParent].child2 = newParent;
        } else root_ = newParent;
        // 3. walk back up refitting bounds and rebalancing
        index = nodes_[leaf].parent;
        while (index != NIL) {
            index = balance(index);
            int c1 = nodes_[index].child1, c2 = nodes_[index].child2;
            nodes_[index].height = 1 + std::max(nodes_[c1].height, nodes_[c2].height);
            nodes_[index].aabb   = nodes_[c1].aabb.merged(nodes_[c2].aabb);
            index = nodes_[index].parent;
        }
    }

    void removeLeaf(int leaf) {
        if (leaf == root_) { root_ = NIL; return; }
        int parent = nodes_[leaf].parent;
        int grand  = nodes_[parent].parent;
        int sibling = (nodes_[parent].child1 == leaf) ? nodes_[parent].child2 : nodes_[parent].child1;
        if (grand != NIL) {
            if (nodes_[grand].child1 == parent) nodes_[grand].child1 = sibling;
            else                                nodes_[grand].child2 = sibling;
            nodes_[sibling].parent = grand;
            freeNode(parent);
            int index = grand;                                  // refit + rebalance up
            while (index != NIL) {
                index = balance(index);
                int c1 = nodes_[index].child1, c2 = nodes_[index].child2;
                nodes_[index].aabb   = nodes_[c1].aabb.merged(nodes_[c2].aabb);
                nodes_[index].height = 1 + std::max(nodes_[c1].height, nodes_[c2].height);
                index = nodes_[index].parent;
            }
        } else {
            root_ = sibling;
            nodes_[sibling].parent = NIL;
            freeNode(parent);
        }
    }

    // AVL-style rotation (faithful port of Box2D b2DynamicTree::Balance): if node A
    // leans by more than one level, promote its taller child. Returns the new root
    // of the subtree that was rooted at A. `refit` helpers keep it readable.
    void refit(int i) {
        int c1 = nodes_[i].child1, c2 = nodes_[i].child2;
        nodes_[i].aabb   = nodes_[c1].aabb.merged(nodes_[c2].aabb);
        nodes_[i].height = 1 + std::max(nodes_[c1].height, nodes_[c2].height);
    }
    int balance(int iA) {
        if (nodes_[iA].isLeaf() || nodes_[iA].height < 2) return iA;
        int iB = nodes_[iA].child1, iC = nodes_[iA].child2;
        int bal = nodes_[iC].height - nodes_[iB].height;

        if (bal > 1) {                         // right child C is taller → rotate C up
            int iF = nodes_[iC].child1, iG = nodes_[iC].child2;
            nodes_[iC].child1 = iA;
            nodes_[iC].parent = nodes_[iA].parent;
            nodes_[iA].parent = iC;
            reparentOldParent(iA, iC);
            if (nodes_[iF].height > nodes_[iG].height) {
                nodes_[iC].child2 = iF; nodes_[iA].child2 = iG; nodes_[iG].parent = iA;
            } else {
                nodes_[iC].child2 = iG; nodes_[iA].child2 = iF; nodes_[iF].parent = iA;
            }
            refit(iA); refit(iC);
            return iC;
        }
        if (bal < -1) {                        // left child B is taller → rotate B up
            int iD = nodes_[iB].child1, iE = nodes_[iB].child2;
            nodes_[iB].child1 = iA;
            nodes_[iB].parent = nodes_[iA].parent;
            nodes_[iA].parent = iB;
            reparentOldParent(iA, iB);
            if (nodes_[iD].height > nodes_[iE].height) {
                nodes_[iB].child2 = iD; nodes_[iA].child1 = iE; nodes_[iE].parent = iA;
            } else {
                nodes_[iB].child2 = iE; nodes_[iA].child1 = iD; nodes_[iD].parent = iA;
            }
            refit(iA); refit(iB);
            return iB;
        }
        return iA;
    }
    // After swapping A under its promoted child `up`, point A's former parent at up.
    void reparentOldParent(int iA, int up) {
        int p = nodes_[up].parent;
        if (p != NIL) {
            if (nodes_[p].child1 == iA) nodes_[p].child1 = up;
            else                        nodes_[p].child2 = up;
        } else root_ = up;
    }
};

} // namespace phys
