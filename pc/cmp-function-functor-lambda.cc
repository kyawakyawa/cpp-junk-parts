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

#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

using namespace std;

static bool cmp_function(double dummy, double x, double y) {
  return x + dummy < y + dummy;
}

int main() {
  struct Comp {
    bool operator()(double dummy, double x, double y) {
      return x + dummy < y + dummy;
    }
  } cmp_functor;

  const size_t N = 10000000;
  vector<double> v(N);
  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());

  std::uniform_real_distribution<> dist(0.0, 1000000000.0);
  for (auto& _v : v) {
    _v = dist(engine);
  }

  std::chrono::system_clock::time_point start_function, end_function;
  double elapsed_function;
  {
    std::vector<double> _v = v;
    start_function         = std::chrono::system_clock::now();
    sort(_v.begin(), _v.end(),
         std::bind(cmp_function, 100.0, std::placeholders::_1,
                   std::placeholders::_2));
    end_function = std::chrono::system_clock::now();
    elapsed_function =
        double(std::chrono::duration_cast<std::chrono::microseconds>(
                   end_function - start_function)
                   .count());
  }

  std::chrono::system_clock::time_point start_functor, end_functor;
  double elapsed_functor;
  {
    std::vector<double> _v = v;
    start_functor          = std::chrono::system_clock::now();
    sort(_v.begin(), _v.end(),
         std::bind(cmp_functor, 100.0, std::placeholders::_1,
                   std::placeholders::_2));
    end_functor = std::chrono::system_clock::now();
    elapsed_functor =
        double(std::chrono::duration_cast<std::chrono::microseconds>(
                   end_functor - start_functor)
                   .count());
  }

  std::chrono::system_clock::time_point start_lambda, end_lambda;
  double elapsed_lambda;
  {
    std::vector<double> _v = v;
    const auto cmp_lambda  = [](double dummy, double x, double y) {
      return x + dummy < y + dummy;
    };
    start_lambda = std::chrono::system_clock::now();
    sort(_v.begin(), _v.end(),
         std::bind(cmp_lambda, 100.0, std::placeholders::_1,
                   std::placeholders::_2));
    end_lambda = std::chrono::system_clock::now();
    elapsed_lambda =
        double(std::chrono::duration_cast<std::chrono::microseconds>(
                   end_lambda - start_lambda)
                   .count());
  }

  printf("function: %.6f\nfunctor : %.6f\nlambda  : %.6f\n",
         elapsed_function / 1000000.0, elapsed_functor / 1000000.0,
         elapsed_lambda / 1000000.0);

  return 0;
}
