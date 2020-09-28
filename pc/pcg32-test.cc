/*
MIT License

Copyright (c) 2019-2020 kyawakyawa

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdint.h>

#include <iostream>
#include <typeinfo>
#ifdef USE_STACK_TRACE_LOGGER
#include <glog/logging.h>
#endif

#include "pcg_random.hpp"

namespace kyawakyawa {
using RNG = pcg32;

// https://en.wikipedia.org/wiki/Permuted_congruential_generator
static uint64_t state = 0x4d595df4d0f33173;  // Or something seed-dependent
static uint64_t const multiplier = 6364136223846793005u;
static uint64_t const increment =
    1442695040888963407u;  // Or an arbitrary odd constant

static uint32_t Rotr32(uint32_t x, unsigned r) {
  return x >> r | x << (-r & 31);
}

static uint32_t Pcg32(void) {
  uint64_t x     = state;
  unsigned count = unsigned(x >> 59);  // 59 = 64 - 5

  state = x * multiplier + increment;
  x ^= x >> 18;                             // 18 = (64 - 27)/2
  return Rotr32(uint32_t(x >> 27), count);  // 27 = 32 - 5
}

static void Pcg32Init(uint64_t seed) {
  state = seed + increment;
  (void)pcg32();
}

static void ShowRandom(void) {
  // RNG rng(42u, 54u);
  Pcg32Init(114514u);
  for (size_t i = 0; i < 50; i++) {
    // uint32_t a = rng();
    const uint32_t a = Pcg32();
    std::cout << float(a) / float(UINT32_MAX) << std::endl;
  }
}
}  // namespace kyawakyawa

int main(int argc, char** argv) {
#ifdef USE_STACK_TRACE_LOGGER
  (void)argc;
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
#else
  (void)argc;
  (void)argv;
#endif
  kyawakyawa::ShowRandom();
  return EXIT_SUCCESS;
}
