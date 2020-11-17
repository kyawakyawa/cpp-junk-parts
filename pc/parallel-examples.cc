#include <stdint.h>

#include <algorithm>
#include <chrono>
#include <numeric>
#include <random>
#include <thread>
#include <vector>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#endif

#include "spdlog/spdlog.h"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

template <class InputIterator, class T>
T ParallelReduce(size_t num_threads, InputIterator first, InputIterator last,
                 T init) {
  num_threads            = std::max<size_t>(1, num_threads);
  const size_t num_terms = size_t(last - first);

  if (num_terms < num_threads) {
    return std::accumulate(first, last, init);
  }
  const long m = long((num_terms + num_threads - 1) / num_threads);

  std::vector<std::thread> workers;
  std::vector<T> sums(num_threads);
  for (size_t thread_id = 0; thread_id < num_threads; ++thread_id) {
    workers.emplace_back([&, thread_id](void) {
      const long offset = long(thread_id) * m;
      if (thread_id == num_threads - 1) {
        sums[thread_id] = std::accumulate(first + offset, last, T(0));
      } else {
        sums[thread_id] =
            std::accumulate(first + offset, first + (offset + m), T(0));
      }
    });
  }
  for (auto& v : workers) {
    v.join();
  }
  return std::accumulate(sums.begin(), sums.end(), init);
}

int main(void) {
  const size_t array_size    = (1 << 26);
  const size_t num_iteration = 10;

  const size_t num_threads = std::thread::hardware_concurrency();

  std::vector<int64_t> vec(array_size);
  std::iota(vec.begin(), vec.end(), 0);
  const int64_t mod = 1000000007;
  for (size_t i = 0; i < num_iteration; ++i) {
    std::random_device rnd;
    std::shuffle(vec.begin(), vec.end(), std::mt19937(rnd()));
    std::for_each(vec.begin(), vec.end(), [](int64_t& u) { u %= mod; });

    const std::chrono::system_clock::time_point start_parallel_reduce =
        std::chrono::system_clock::now();
    const int64_t sum_pr =
        ParallelReduce(num_threads, vec.begin(), vec.end(), int64_t(0));
    const std::chrono::system_clock::time_point end_parallel_reduce =
        std::chrono::system_clock::now();

    const std::chrono::system_clock::time_point start_accumulate =
        std::chrono::system_clock::now();
    const int64_t sum_ac = std::accumulate(vec.begin(), vec.end(), int64_t(0));
    const std::chrono::system_clock::time_point end_accumulate =
        std::chrono::system_clock::now();

    if (sum_pr != sum_ac) {
      spdlog::error("ParallelReduce -> {}, std::accumulate -> {}", sum_pr,
                    sum_ac);
    }
    const auto time_parallel_reduce =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            end_parallel_reduce - start_parallel_reduce)
            .count();
    const auto time_accumulate =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end_accumulate -
                                                             start_accumulate)
            .count();

    spdlog::info("parallel reduce: {} ms\n",
                 double(time_parallel_reduce) / 1000000);
    spdlog::info("std::accumulate: {} ms\n", double(time_accumulate) / 1000000);
  }

  return 0;
}
