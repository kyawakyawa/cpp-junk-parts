/*
MIT License

Copyright (c) 2020 kyawakyawa

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

#include <algorithm>
#include <bit>
#include <bitset>
#include <chrono>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

static uint64_t Xor128Sub(void) {
  static uint64_t x = 123456789, y = 362436069, z = 521288629, w = 88675123;
  uint64_t t;
  t = (x ^ (x << 11));
  x = y;
  y = z;
  z = w;
  return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
}
static uint64_t Xor128(void) {
  const uint64_t x = Xor128Sub();
  const uint64_t y = Xor128Sub();
  return (x << uint64_t(32)) + y;
}

struct alignas(256) Int256 {
  uint64_t v[4];
  inline Int256 operator^(const Int256& o) const {
    Int256 ret;
    ret.v[0] = (v[0] ^ o.v[0]);
    ret.v[1] = (v[1] ^ o.v[1]);
    ret.v[2] = (v[2] ^ o.v[2]);
    ret.v[3] = (v[3] ^ o.v[3]);

    return ret;
  }
};

static int SimplePopCount(uint64_t bt) {
  int count = 0;
  for (; bt; bt &= bt - 1) ++count;
  return count;
}

static int PopCount(uint64_t bt) {
  const uint64_t a =
      0b0101'0101'0101'0101'0101'0101'0101'0101'0101'0101'0101'0101'0101'0101'0101'0101;
  const uint64_t b =
      0b1010'1010'1010'1010'1010'1010'1010'1010'1010'1010'1010'1010'1010'1010'1010'1010;

  const uint64_t c =
      0b0011'0011'0011'0011'0011'0011'0011'0011'0011'0011'0011'0011'0011'0011'0011'0011;
  const uint64_t d =
      0b1100'1100'1100'1100'1100'1100'1100'1100'1100'1100'1100'1100'1100'1100'1100'1100;

  const uint64_t e =
      0b0000'1111'0000'1111'0000'1111'0000'1111'0000'1111'0000'1111'0000'1111'0000'1111;
  const uint64_t f =
      0b1111'0000'1111'0000'1111'0000'1111'0000'1111'0000'1111'0000'1111'0000'1111'0000;

  const uint64_t g =
      0b0000'0000'1111'1111'0000'0000'1111'1111'0000'0000'1111'1111'0000'0000'1111'1111;
  const uint64_t h =
      0b1111'1111'0000'0000'1111'1111'0000'0000'1111'1111'0000'0000'1111'1111'0000'0000;

  const uint64_t i =
      0b0000'0000'0000'0000'1111'1111'1111'1111'0000'0000'0000'0000'1111'1111'1111'1111;
  const uint64_t j =
      0b1111'1111'1111'1111'0000'0000'0000'0000'1111'1111'1111'1111'0000'0000'0000'0000;

  const uint64_t k =
      0b0000'0000'0000'0000'0000'0000'0000'0000'1111'1111'1111'1111'1111'1111'1111'1111;
  const uint64_t l =
      0b1111'1111'1111'1111'1111'1111'1111'1111'0000'0000'0000'0000'0000'0000'0000'0000;

  bt = (bt & a) + ((bt & b) >> 1);
  bt = (bt & c) + ((bt & d) >> 2);
  bt = (bt & e) + ((bt & f) >> 4);
  bt = (bt & g) + ((bt & h) >> 8);
  bt = (bt & i) + ((bt & j) >> 16);
  return int((bt & k) + ((bt & l) >> 32));
}

static int PopCount(const Int256& bt) {
  return PopCount(bt.v[0]) + PopCount(bt.v[1]) + PopCount(bt.v[2]) +
         PopCount(bt.v[3]);
}
static int PopCountCpp20(const Int256& bt) {
  return std::popcount(bt.v[0]) + std::popcount(bt.v[1]) +
         std::popcount(bt.v[2]) + std::popcount(bt.v[3]);
}

static std::bitset<256> Int256ToBitset(const Int256& bt) {
  std::bitset<256> ret(0);

  for (int i = 0; i < 64; ++i) {
    ret[size_t(i)] = ((bt.v[0] >> i) & 1);
  }

  for (int i = 0; i < 64; ++i) {
    ret[size_t(i + 64)] = ((bt.v[1] >> i) & 1);
  }

  for (int i = 0; i < 64; ++i) {
    ret[size_t(i + 128)] = ((bt.v[2] >> i) & 1);
  }

  for (int i = 0; i < 64; ++i) {
    ret[size_t(i + 192)] = ((bt.v[3] >> i) & 1);
  }

  return ret;
}

static bool TestPopCount(void) {
  const size_t num_test = 100000;

  for (size_t i = 0; i < num_test; ++i) {
    Int256 bt = {};

    bt.v[0] = Xor128();
    bt.v[1] = Xor128();
    bt.v[2] = Xor128();
    bt.v[3] = Xor128();

    const std::bitset<256> bs = Int256ToBitset(bt);

    const int correct = std::accumulate(bt.v, bt.v + 4, 0,
                                        [](const int& acc, const uint64_t& b) {
                                          return acc + SimplePopCount(b);
                                        });

    if (correct != PopCount(bt) || correct != PopCountCpp20(bt) ||
        correct != int(bs.count())) {
      std::cerr << "[Error]!" << std::endl;
      return false;
    }
  }
  std::cerr << "[Pass]!" << std::endl;
  return true;
}

static void BenchPopCount(void) {
  const size_t num_exe = 1000000;
  std::vector<std::pair<Int256, Int256>> bts(num_exe);
  std::for_each(bts.begin(), bts.end(), [](std::pair<Int256, Int256>& pbt) {
    pbt.first.v[0] = Xor128();
    pbt.first.v[1] = Xor128();
    pbt.first.v[2] = Xor128();
    pbt.first.v[3] = Xor128();

    pbt.second.v[0] = Xor128();
    pbt.second.v[1] = Xor128();
    pbt.second.v[2] = Xor128();
    pbt.second.v[3] = Xor128();
  });
  std::vector<std::pair<std::bitset<256>, std::bitset<256>>> bss(num_exe);
  std::transform(bts.begin(), bts.end(), bss.begin(),
                 [](const std::pair<Int256, Int256>& pbt) {
                   return std::make_pair(Int256ToBitset(pbt.first),
                                         Int256ToBitset(pbt.second));
                 });

  std::vector<int> result;

  printf("[Info] start hamming dist bench (%lu iteration)\n", num_exe);
  result.reserve(num_exe * 3);
  const std::chrono::system_clock::time_point start_popcount =
      std::chrono::system_clock::now();
  std::for_each(bts.begin(), bts.end(),
                [&result](const std::pair<Int256, Int256>& pbt) {
                  const Int256 tmp = pbt.first ^ pbt.second;
                  result.emplace_back(PopCount(tmp));
                });
  const std::chrono::system_clock::time_point end_popcount =
      std::chrono::system_clock::now();

  const std::chrono::system_clock::time_point start_popcount_cpp20 =
      std::chrono::system_clock::now();
  std::for_each(bts.begin(), bts.end(),
                [&result](const std::pair<Int256, Int256>& pbt) {
                  const Int256 tmp = pbt.first ^ pbt.second;
                  result.emplace_back(PopCountCpp20(tmp));
                });
  const std::chrono::system_clock::time_point end_popcount_cpp20 =
      std::chrono::system_clock::now();

  const std::chrono::system_clock::time_point start_popcount_bitset =
      std::chrono::system_clock::now();
  std::for_each(
      bss.begin(), bss.end(),
      [&result](const std::pair<std::bitset<256>, std::bitset<256>>& pbs) {
        const std::bitset<256> tmp = pbs.first ^ pbs.second;
        result.emplace_back(int(tmp.count()));
      });
  const std::chrono::system_clock::time_point end_popcount_bitset =
      std::chrono::system_clock::now();
  printf("[Info] end bench\n");

  printf("[Info] start check\n");
  for (size_t i = 0; i < num_exe; ++i) {
    if (result[i] != result[i + num_exe] ||
        result[i + num_exe] != result[i + 2 * num_exe] ||
        result[i + 2 * num_exe] != result[i]) {
      std::cerr << "[Error]" << std::endl;
    }
  }
  printf("[Info] end check\n");

  const auto time_popcount =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end_popcount -
                                                           start_popcount)
          .count();
  const auto time_popcount_cpp20 =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end_popcount_cpp20 -
                                                           start_popcount_cpp20)
          .count();
  const auto time_popcount_bitset =
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          end_popcount_bitset - start_popcount_bitset)
          .count();

  printf("      my popcount: %3.6f ms\n", double(time_popcount) / 1000000);
  printf("   c++20 popcount: %3.6f ms\n",
         double(time_popcount_cpp20) / 1000000);
  printf("std::bitset count: %3.6f ms\n",
         double(time_popcount_bitset) / 1000000);
}

int main(void) {
  TestPopCount();
  BenchPopCount();

  return 0;
}
