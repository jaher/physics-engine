// Minimal PNG writer (RGB8, via zlib) so the demos can save screenshots without
// an image library. Rows are flipped from OpenGL's bottom-up convention.
#pragma once
#include <zlib.h>
#include <vector>
#include <cstdio>
#include <cstdint>

namespace gfx {
inline void put32(std::vector<uint8_t>& v, uint32_t x) { v.push_back(x >> 24); v.push_back(x >> 16); v.push_back(x >> 8); v.push_back(x); }
inline void chunk(std::vector<uint8_t>& out, const char* type, const std::vector<uint8_t>& data) {
    put32(out, (uint32_t)data.size());
    size_t s = out.size();
    out.insert(out.end(), type, type + 4);
    out.insert(out.end(), data.begin(), data.end());
    uint32_t crc = (uint32_t)crc32(0, out.data() + s, (uInt)(4 + data.size()));
    put32(out, crc);
}
inline bool writePNG(const char* fn, const uint8_t* rgb, int w, int h, bool flip = true) {
    std::vector<uint8_t> raw; raw.reserve((size_t)(w * 3 + 1) * h);
    for (int y = 0; y < h; y++) { int sy = flip ? (h - 1 - y) : y; raw.push_back(0);
        const uint8_t* row = rgb + (size_t)sy * w * 3; raw.insert(raw.end(), row, row + w * 3); }
    uLongf clen = compressBound((uLong)raw.size()); std::vector<uint8_t> comp(clen);
    if (compress2(comp.data(), &clen, raw.data(), (uLong)raw.size(), 9) != Z_OK) return false;
    comp.resize(clen);
    std::vector<uint8_t> out; const uint8_t sig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
    out.insert(out.end(), sig, sig + 8);
    std::vector<uint8_t> ihdr; put32(ihdr, w); put32(ihdr, h);
    ihdr.push_back(8); ihdr.push_back(2); ihdr.push_back(0); ihdr.push_back(0); ihdr.push_back(0);
    chunk(out, "IHDR", ihdr); chunk(out, "IDAT", comp); chunk(out, "IEND", {});
    FILE* f = fopen(fn, "wb"); if (!f) return false;
    fwrite(out.data(), 1, out.size(), f); fclose(f); return true;
}
} // namespace gfx
