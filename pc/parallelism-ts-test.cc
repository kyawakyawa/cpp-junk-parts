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

#include <execution>
#include <numeric>
#include <random>
#include <vector>

#ifdef USE_STACK_TRACE_LOGGER
#include <glog/logging.h>
#endif

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
#ifdef USE_STACK_TRACE_LOGGER
  (void)argc;
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
#endif

  const size_t array_size = (1 << 20);
  const size_t num_iteration = 1000;

  std::vector<int> vec(array_size);
  std::iota(vec.begin(), vec.end(), 0);
  std::random_device rnd;
  std::shuffle(vec.begin(), vec.end(), std::mt19937(rnd()));
  for (size_t i = 0; i < num_iteration; ++i) {
    std::vector<int> _vec(vec.size());
    std::copy(std::execution::par_unseq, vec.begin(), vec.end(), _vec.begin());
    std::sort(std::execution::par_unseq, _vec.begin(), _vec.end());
  }
  return 0;
}
