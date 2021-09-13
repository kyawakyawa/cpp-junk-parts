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

#include <assert.h>
#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#ifdef USE_STACK_TRACE_LOGGER
#include <glog/logging.h>
#endif

static std::vector<int> GenerateRandomVector(const size_t n) {
  std::random_device seed_gen;
  std::default_random_engine engine(seed_gen());

  std::uniform_int_distribution<> dist(1, 100);

  std::vector<int> ret(n);
  std::for_each(ret.begin(), ret.end(),
                [&](int &a) { return a = dist(engine); });

  return ret;
}

// get floor((n - 1) / 2)th element
//  TODO 長さ0の時
static int MedianWithSort(const std::vector<int> &input) {
  const size_t n = input.size();
  auto v = input;
  std::sort(v.begin(), v.end());
  return v[(n - 1) / 2];
}

enum SelectionMethod { MedianOfMedian, Randomized };

//  TODO 長さ0の時
static int MedianWithQuickSelect(
    const std::vector<int> &input,
    SelectionMethod sm = SelectionMethod::Randomized) {
  const size_t n = input.size();
  auto v = input;

  const size_t m_id = (n - 1) / 2;

  int l = 0;
  int r = int(n);
  // [l, r)

  int ret = 0;
  while (r - l > 0) {
    size_t pivot_id = 0;
    if (sm == Randomized) {
      // NOLINTNEXTLINE
      pivot_id = size_t(l + (rand() % (r - l)));
    }
    // TODO 中央値の中央値を実装する // http://www.flint.jp/blog/?entry=109

    // ピボットを先頭に退避
    const int pivot = v[pivot_id];
    std::swap(v[size_t(l)], v[pivot_id]);

    const long c =
        (std::partition(v.begin() + l + 1, v.begin() + r,
                        [&pivot](const int x) { return x < pivot; }) -
         v.begin()) -
        1;

    // ピボットが後半グループの先頭になるようにする
    std::swap(v[size_t(l)], v[size_t(c)]);

    if (size_t(c) == m_id) {
      ret = v[size_t(c)];
      break;
    }
    if (m_id < size_t(c)) {
      r = int(c);
    } else {
      l = int(c + 1);
    }
  }
  return ret;
}

static int NthElement(const std::vector<int> &input) {
  auto v = input;
  const size_t n = v.size();
  const size_t m_id = (n - 1) / 2;
  std::nth_element(v.begin(), v.begin() + long(m_id), v.end());
  return v[m_id];
}

static void Test(const size_t n, const bool output = true) {
  const auto a = GenerateRandomVector(n);
  if (output) {
    std::cout << "input:\n\t";
    for (const auto v : a) {
      std::cout << v << " ";
    }
    std::cout << std::endl;

    std::cout << "quicksort  : " << MedianWithSort(a) << std::endl;
    std::cout << "quickselct : " << MedianWithQuickSelect(a) << std::endl;
    std::cout << "nth_element : " << NthElement(a) << std::endl;
  }

  int s = MedianWithSort(a);
  int q = MedianWithQuickSelect(a);
  int nth = MedianWithQuickSelect(a);

  if (s != q || q != nth || s != nth) {
    std::runtime_error("error is occured");
  }
}

int main(int argc, char **argv) {
#ifdef USE_STACK_TRACE_LOGGER
  (void)argc;
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
#else
  (void)argc;
  (void)argv;
#endif
  size_t n = 10;
  size_t it = 10;

  for (size_t i = 0; i < it; ++i) {
    Test(n);
    std::cout << std::endl;
  }

  std::random_device seed_gen;
  std::default_random_engine engine(seed_gen());

  std::uniform_int_distribution<> dist(1, 5000);
  it = 100000;

  for (size_t i = 0; i < it; ++i) {
    if (i % std::max<size_t>(1, it / 100) == 0) {
      std::cout << 100. * double(i) / double(it) << "%" << std::endl;
    }
    n = size_t(dist(engine));
    Test(n, false);
  }
  std::cout << "pass all tests" << std::endl;
  return 0;
}
