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

#include <math.h>

#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

// NOLINTNEXTLINE
using namespace std;

static bool cmp_function(double dummy, double x, double y) {
  return x + dummy < y + dummy;
}

int main() {
  struct Comp {
    [[maybe_unused]] bool operator()(double dummy, double x, double y) {
      return x + dummy < y + dummy;
    }
  } cmp_functor;

  const size_t N = 10000000;
  vector<double> v(N);
  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());

  const double max_value = 1000000000.0;
  std::uniform_real_distribution<> dist(0.0, max_value);
  for (auto& _v : v) {
    _v = dist(engine);
  }

  std::chrono::high_resolution_clock::time_point start_function;
  std::chrono::high_resolution_clock::time_point end_function;
  double elapsed_function = 0.0;
  {
    std::vector<double> _v = v;
    start_function = std::chrono::high_resolution_clock::now();
    sort(_v.begin(), _v.end(),
         std::bind(cmp_function, 100.0, std::placeholders::_1,
                   std::placeholders::_2));
    end_function = std::chrono::high_resolution_clock::now();
    elapsed_function =
        double(std::chrono::duration_cast<std::chrono::microseconds>(
                   end_function - start_function)
                   .count());
  }

  std::chrono::high_resolution_clock::time_point start_functor;
  std::chrono::high_resolution_clock::time_point end_functor;
  double elapsed_functor = 0.0;
  {
    std::vector<double> _v = v;
    start_functor = std::chrono::high_resolution_clock::now();
    sort(_v.begin(), _v.end(),
         std::bind(cmp_functor, 100.0, std::placeholders::_1,
                   std::placeholders::_2));
    end_functor = std::chrono::high_resolution_clock::now();
    elapsed_functor =
        double(std::chrono::duration_cast<std::chrono::microseconds>(
                   end_functor - start_functor)
                   .count());
  }

  std::chrono::high_resolution_clock::time_point start_lambda;
  std::chrono::high_resolution_clock::time_point end_lambda;
  double elapsed_lambda = 0.0;
  {
    std::vector<double> _v = v;
    const auto cmp_lambda = [](double dummy, double x, double y) {
      return x + dummy < y + dummy;
    };
    start_lambda = std::chrono::high_resolution_clock::now();
    sort(_v.begin(), _v.end(),
         std::bind(cmp_lambda, 100.0, std::placeholders::_1,
                   std::placeholders::_2));
    end_lambda = std::chrono::high_resolution_clock::now();
    elapsed_lambda =
        double(std::chrono::duration_cast<std::chrono::microseconds>(
                   end_lambda - start_lambda)
                   .count());
  }

  const double micro_seconds_to_seconds = 10000000.0;
  // NOLINTNEXTLINE
  printf("function: %.6f\nfunctor : %.6f\nlambda  : %.6f\n",
         elapsed_function / micro_seconds_to_seconds,
         elapsed_functor / micro_seconds_to_seconds,
         elapsed_lambda / micro_seconds_to_seconds);

  return 0;
}
