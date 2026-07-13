// A tiny robot-description loader: parses a SMALL subset of URDF (ROS/Gazebo) or
// MJCF (MuJoCo) straight from a std::string into an in-memory Model, with NO
// external XML library — just a hand-rolled tag/attribute scanner. Enough to
// pull out the links (mass, diagonal inertia, one box/sphere/cylinder geom) and
// the joints (type, parent/child, axis, origin) that the robotics layer needs.
#pragma once
#include "core.h"
#include <string>
#include <vector>
#include <sstream>
#include <cctype>

namespace phys {

// ------------------------------------------------------------------ model types
enum class GeomType  { None, Box, Sphere, Cylinder };
enum class JointKind { Fixed, Revolute, Prismatic };

struct Link {
    std::string name;
    real mass = 0;
    Vector3 inertiaDiag;                 // (ixx, iyy, izz) about the COM
    GeomType geom = GeomType::None;
    Vector3 size;                        // box: extents; sphere: radius→.x; cylinder: (radius, length)
    Vector3 origin;                      // frame offset from parent (used by MJCF `pos`)
};

struct ModelJoint {
    std::string name;
    JointKind type = JointKind::Fixed;
    std::string parent, child;           // link names
    Vector3 axis = Vector3(0, 0, 1);
    Vector3 origin;                      // xyz translation of the joint frame
};

struct Model {
    std::vector<Link>  links;
    std::vector<ModelJoint> joints;
    const Link*  findLink (const std::string& n) const { for (auto& l : links)  if (l.name == n) return &l; return nullptr; }
    const ModelJoint* findJoint(const std::string& n) const { for (auto& j : joints) if (j.name == n) return &j; return nullptr; }
};

// ---------------------------------------------------------------- tiny XML scan
// One parsed tag: name, attribute list, and whether it's a close (</x>) or a
// self-closing (<x/>) element. Text nodes are ignored — the subset is attribute
// driven, as both URDF and MJCF are.
struct XmlTag {
    std::string name;
    std::vector<std::pair<std::string, std::string>> attrs;
    bool close = false, selfClose = false;
    std::string attr(const std::string& k, const std::string& d = "") const {
        for (auto& a : attrs) if (a.first == k) return a.second;
        return d;
    }
    bool has(const std::string& k) const { for (auto& a : attrs) if (a.first == k) return true; return false; }
};

inline std::string loader_stripComments(const std::string& s) {
    std::string out; out.reserve(s.size());
    for (size_t i = 0; i < s.size(); ) {
        if (s.compare(i, 4, "<!--") == 0) {                 // drop <!-- ... -->
            size_t e = s.find("-->", i + 4);
            if (e == std::string::npos) break;
            i = e + 3;
        } else out.push_back(s[i++]);
    }
    return out;
}

inline std::vector<XmlTag> loader_scan(const std::string& raw) {
    std::string xml = loader_stripComments(raw);
    std::vector<XmlTag> tags;
    auto sp = [](char c) { return std::isspace((unsigned char)c) != 0; };
    for (size_t i = 0, n = xml.size(); i < n; ) {
        if (xml[i] != '<') { i++; continue; }
        size_t end = xml.find('>', i);
        if (end == std::string::npos) break;
        std::string b = xml.substr(i + 1, end - i - 1);
        i = end + 1;
        if (b.empty() || b[0] == '?' || b[0] == '!') continue;   // <?xml?>, <!doctype>
        XmlTag t;
        size_t p = 0;
        if (b[p] == '/') { t.close = true; p++; }
        while (p < b.size() && sp(b[p])) p++;
        size_t ns = p;
        while (p < b.size() && !sp(b[p]) && b[p] != '/') p++;
        t.name = b.substr(ns, p - ns);
        while (p < b.size()) {                              // attributes
            while (p < b.size() && sp(b[p])) p++;
            if (p >= b.size()) break;
            if (b[p] == '/') { t.selfClose = true; break; }
            size_t ks = p;
            while (p < b.size() && b[p] != '=' && !sp(b[p]) && b[p] != '/') p++;
            std::string key = b.substr(ks, p - ks);
            while (p < b.size() && sp(b[p])) p++;
            if (p < b.size() && b[p] == '=') {
                p++;
                while (p < b.size() && sp(b[p])) p++;
                char q = (p < b.size()) ? b[p] : 0;
                if (q == '"' || q == '\'') {
                    p++; size_t vs = p;
                    while (p < b.size() && b[p] != q) p++;
                    t.attrs.push_back({key, b.substr(vs, p - vs)});
                    if (p < b.size()) p++;
                } else {
                    size_t vs = p;
                    while (p < b.size() && !sp(b[p]) && b[p] != '/') p++;
                    t.attrs.push_back({key, b.substr(vs, p - vs)});
                }
            } else if (!key.empty()) t.attrs.push_back({key, ""});
        }
        tags.push_back(t);
    }
    return tags;
}

// --------------------------------------------------------------- attr → numbers
inline real loader_real(const std::string& s, real d = 0) {
    if (s.empty()) return d;
    try { return (real)std::stod(s); } catch (...) { return d; }
}
// Parse "x y z" filling as many components as present (radius-only sizes → .x).
inline Vector3 loader_vec3(const std::string& s, const Vector3& d = Vector3()) {
    std::istringstream is(s);
    real x;
    if (!(is >> x)) return d;
    Vector3 v(x, 0, 0);
    real y, z;
    if (is >> y) { v.y = y; if (is >> z) v.z = z; }
    return v;
}
inline JointKind loader_jointType(const std::string& s) {
    if (s == "revolute" || s == "continuous" || s == "hinge") return JointKind::Revolute;
    if (s == "prismatic" || s == "slide")                     return JointKind::Prismatic;
    return JointKind::Fixed;                                  // fixed / free / ball / floating → fixed
}
inline void loader_setGeom(Link& L, const XmlTag& t, const std::string& type) {
    if (type == "box")      { L.geom = GeomType::Box;      L.size = loader_vec3(t.attr("size")); }
    else if (type == "sphere") { L.geom = GeomType::Sphere; L.size = Vector3(loader_real(t.attr("radius", t.attr("size"))), 0, 0); }
    else if (type == "cylinder" || type == "capsule") {
        L.geom = GeomType::Cylinder;
        if (t.has("size")) L.size = loader_vec3(t.attr("size"));   // MJCF: radius [half-length]
        else L.size = Vector3(loader_real(t.attr("radius")), loader_real(t.attr("length")), 0);
    }
}

// ------------------------------------------------------------------- loadURDF
// URDF: <robot> holds <link name>{<inertial><mass value/><inertia ixx iyy izz/>
// <geometry><box/sphere/cylinder/></geometry>} and <joint name type>{<parent
// link/><child link/><axis xyz/><origin xyz/>}. Indices (not pointers) keep the
// "current" link/joint valid across vector reallocations.
inline Model loadURDF(const std::string& xml) {
    Model m;
    int li = -1, ji = -1;                                   // current link / joint index
    for (const XmlTag& t : loader_scan(xml)) {
        if (t.close) {
            if (t.name == "link")  li = -1;
            else if (t.name == "joint") ji = -1;
            continue;
        }
        if (t.name == "link") {
            m.links.push_back(Link{}); li = (int)m.links.size() - 1;
            m.links[li].name = t.attr("name");
            if (t.selfClose) li = -1;
        } else if (t.name == "joint") {
            m.joints.push_back(ModelJoint{}); ji = (int)m.joints.size() - 1;
            m.joints[ji].name = t.attr("name");
            m.joints[ji].type = loader_jointType(t.attr("type"));
            if (t.selfClose) ji = -1;
        } else if (li >= 0 && t.name == "mass") {
            m.links[li].mass = loader_real(t.attr("value"));
        } else if (li >= 0 && t.name == "inertia") {
            m.links[li].inertiaDiag = Vector3(loader_real(t.attr("ixx")),
                                              loader_real(t.attr("iyy")),
                                              loader_real(t.attr("izz")));
        } else if (li >= 0 && (t.name == "box" || t.name == "sphere" || t.name == "cylinder")) {
            loader_setGeom(m.links[li], t, t.name);
        } else if (ji >= 0 && t.name == "parent") {
            m.joints[ji].parent = t.attr("link");
        } else if (ji >= 0 && t.name == "child") {
            m.joints[ji].child = t.attr("link");
        } else if (ji >= 0 && t.name == "axis") {
            m.joints[ji].axis = loader_vec3(t.attr("xyz"), Vector3(1, 0, 0));
        } else if (ji >= 0 && t.name == "origin") {
            m.joints[ji].origin = loader_vec3(t.attr("xyz"));
        }
    }
    return m;
}

// ------------------------------------------------------------------- loadMJCF
// MJCF: nested <body name pos> elements give the kinematic tree; each body's
// connecting <joint name type axis> ties it to the enclosing body, and a <geom
// type size mass> / <inertial mass diaginertia> carries its shape and inertia.
// The body-name stack recovers parent/child from the nesting.
inline Model loadMJCF(const std::string& xml) {
    Model m;
    std::vector<int> stack;                                 // enclosing body indices
    for (const XmlTag& t : loader_scan(xml)) {
        if (t.close) { if (t.name == "body" && !stack.empty()) stack.pop_back(); continue; }
        if (t.name == "body") {
            m.links.push_back(Link{}); int li = (int)m.links.size() - 1;
            m.links[li].name   = t.attr("name");
            m.links[li].origin = loader_vec3(t.attr("pos"));
            if (t.has("mass")) m.links[li].mass = loader_real(t.attr("mass"));
            stack.push_back(li);
            if (t.selfClose) stack.pop_back();
        } else if (!stack.empty() && t.name == "geom") {
            Link& L = m.links[stack.back()];
            loader_setGeom(L, t, t.attr("type", "sphere"));
            if (t.has("mass")) L.mass = loader_real(t.attr("mass"));
        } else if (!stack.empty() && t.name == "inertial") {
            Link& L = m.links[stack.back()];
            if (t.has("mass")) L.mass = loader_real(t.attr("mass"));
            L.inertiaDiag = loader_vec3(t.attr("diaginertia"), L.inertiaDiag);
        } else if (!stack.empty() && t.name == "joint") {
            m.joints.push_back(ModelJoint{}); ModelJoint& j = m.joints.back();
            j.name   = t.attr("name");
            j.type   = loader_jointType(t.attr("type", "hinge"));   // MJCF default joint is a hinge
            j.child  = m.links[stack.back()].name;
            j.parent = stack.size() >= 2 ? m.links[stack[stack.size() - 2]].name : "world";
            j.axis   = loader_vec3(t.attr("axis"), Vector3(0, 0, 1));
            j.origin = m.links[stack.back()].origin;
        }
    }
    return m;
}

} // namespace phys
