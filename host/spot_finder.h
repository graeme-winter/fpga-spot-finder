#pragma once

#include <cstdint>
#include <hls_stream.h>
#include <hls_streamofblocks.h>

#define NY 16384
#define NX 1024 // must be multiple of PACK

constexpr int PACK = 4; // must be 2^n
constexpr int KNL = 3;
constexpr int KNL2 = 2 * KNL + 2;

struct pack8_t {
  uint8_t data[PACK];
  uint8_t &operator[](int i) { return data[i]; }
};

struct pack16_t {
  uint16_t data[PACK];
  uint16_t &operator[](int i) { return data[i]; }
};

struct pack32_t {
  uint32_t data[PACK];
  uint32_t &operator[](int i) { return data[i]; }
};

typedef pack8_t block8_t[NX / PACK];
typedef pack16_t block16_t[NX / PACK];
typedef pack32_t block32_t[NX / PACK];

typedef hls::stream<pack8_t> pack8_stream;
typedef hls::stream<pack16_t> pack16_stream;
typedef hls::stream<pack32_t> pack32_stream;

typedef hls::stream_of_blocks<block8_t> block8_stream;
typedef hls::stream_of_blocks<block16_t> block16_stream;
typedef hls::stream_of_blocks<block32_t> block32_stream;
