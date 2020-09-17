#include <algorithm>
#include <assert.h>
#include <iostream>
#include <random>
#include <stdexcept>
#include <stdlib.h>
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
//  TODO 0の時
static int MedianWithSort(const std::vector<int> &input) {
  const size_t n = input.size();
  auto v = input;
  std::sort(v.begin(), v.end());
  return v[(n - 1) / 2];
}

enum SelectionMethod { MedianOfMedian, Randomized };

//  TODO 0の時
static int
MedianWithQuickSelect(const std::vector<int> &input,
                      SelectionMethod sm = SelectionMethod::Randomized) {
  const size_t n = input.size();
  auto v = input;

  const size_t m_id = (n - 1) / 2;

  int l = 0;
  int r = int(n);
  // [l, r)

  int ret = 0;
  while (true) {
    size_t pivot_id = 0;
    if (sm == Randomized) {
      pivot_id = size_t(l + (rand() % (r - l)));
    }

    const int pivot = v[pivot_id];
    std::swap(v[size_t(l)], v[pivot_id]);

    const long c =
        (std::partition(v.begin() + l + 1, v.begin() + r,
                        [&pivot](const int x) { return x < pivot; }) -
         v.begin()) -
        1;

    std::swap(v[size_t(l)], v[size_t(c)]);

    if (size_t(c) == m_id) {
      ret = v[size_t(c)];
      break;
    } else if (m_id < size_t(c)) {
      r = int(c);
    } else {
      l = int(c + 1);
    }
  }
  return ret;
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
  }

  if (MedianWithSort(a) != MedianWithQuickSelect(a)) {
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
  it = 10000;

  for (size_t i = 0; i < it; ++i) {
    Test(n, false);
  }
  std::cout << "pass all tests" << std::endl;
  return 0;
}
