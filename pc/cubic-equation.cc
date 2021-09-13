/*
MIT License

Copyright (c) 2019-2021 kyawakyawa

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
#include <math.h>

#include <algorithm>
#include <chrono>
#include <random>
#include <vector>

template <typename T = double>
T EvaluateCubicFunction(const T b, const T c, const T d, const T x) {
  return ((x + b) * x + c) * x + d;
}

template <typename T = double>
T NewtonMethodForCubicEquation(const T b, const T c, const T d, T x,
                               const T thr = 1e-12,
                               const size_t max_itr = 10000) {
  const T a2 = T(2.0) * b;
  for (size_t itr = 0; itr < max_itr; ++itr) {
    const T x2 = x * x;

    const T f = x2 * x + b * x2 + c * x + d;
    const T f_dash = T(3) * x2 + a2 * x + c;

    x -= f / f_dash;
    if (abs(EvaluateCubicFunction(b, c, d, x)) < thr) {
      break;
    }
  }
  return x;
}

template <typename T = double>
inline void ComputeRealSolutionOfCubicEquation(
    T b, const T c, const T d, T* solutions, size_t* num_solutions,
    T D_thr = 1e-14, T multiple_solution_thr = 1e-14) {
  // x^3 + b x^2  + c x + d =  0;
  //             |
  //             v
  // x^3 + 3b x^2  + c x + d =  0;
#if 0
  const T original_b = b;
#endif
  b *= (T(1) / T(3));

  const T b2 = b * b;

  const T p = b2 - c * (T(1) / T(3));
  const T q = (b * c - T(2.0) * b * b2 - d) * T(0.5);

  // y^3 - 3 p y - 2 q = 0
  //       where y = x + b

  const T p3 = p * p * p;
  const T D = q * q - p3;

  if (D < -D_thr) {
    assert(p > T(0.0));
    const T sq_p_2 = sqrt(p /* FIXME */) * T(2);
#if 0
    solutions[0] = NewtonMethodForCubicEquation(original_b, c, d, -sq_p_2 - b);
    solutions[1] = NewtonMethodForCubicEquation(original_b, c, d, -b);
    solutions[2] = NewtonMethodForCubicEquation(original_b, c, d, sq_p_2 - b);
#else
    const T theta = acos(q * std::pow(p, -T(3) / T(2)));
    constexpr T kPi2 = T(2.0 * 3.141592653589793238462643383279);
    constexpr T kPi4 = T(4.0 * 3.141592653589793238462643383279);
    constexpr T three_div = T(1) / T(3);
    solutions[0] = sq_p_2 * cos(theta * three_div) - b;
    solutions[1] = sq_p_2 * cos((theta + kPi2) * three_div) - b;
    solutions[2] = sq_p_2 * cos((theta + kPi4) * three_div) - b;

#endif

    *num_solutions = 3;
  } else {
    const T sqD = (D < D_thr ? T(0.0) : sqrt(D));
    const T alpha_3 =
        (!std::signbit(q)
             ? cbrt(q + sqD)
             : cbrt(q - sqD));  // alpha3 の絶対値が大きくなるように
    const T beta_3 =
        (D < D_thr ? alpha_3
                   : p / alpha_3);  // 重解のときはゼロ割に近くなるので場合分け

    const T y1 = alpha_3 + beta_3;

    solutions[0] = y1 - b;
    // ret.emplace_back(NewtonMethodForCubicEquation(original_b, c, d, y1 - b));

    *num_solutions = 1;

    if (D < D_thr) {
      const T y2_cad = y1 * T(-0.5);
      if (abs(y2_cad) > multiple_solution_thr) {
        solutions[1] = y2_cad - b;
        // ret.emplace_back(
        //     NewtonMethodForCubicEquation(original_b, c, d, y2_cad - b));

        ++(*num_solutions);
      }
    }
  }
}

static bool Test(const double b, const double c, const double d) {
  double solutions[3];
  size_t num_solutions;

  ComputeRealSolutionOfCubicEquation(b, c, d, solutions, &num_solutions);

  if (!(0 < num_solutions && num_solutions < 4)) {
    return false;
  }

  for (size_t i = 0; i < num_solutions; ++i) {
    if (abs(EvaluateCubicFunction(b, c, d, solutions[i])) > 1e-12) {
      printf("error %e\n", abs(EvaluateCubicFunction(b, c, d, solutions[i])));
      return false;
    }
  }
  return true;
}

static long Benchmark(const double b, const double c, const double d) {
  std::chrono::high_resolution_clock::time_point start;
  std::chrono::high_resolution_clock::time_point end;

  double solutions[3];
  size_t num_solutions;

  start = std::chrono::high_resolution_clock::now();
  ComputeRealSolutionOfCubicEquation(b, c, d, solutions, &num_solutions);
  end = std::chrono::high_resolution_clock::now();

  const long com_time =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

  return com_time;
}

static void TryAndOutput(const double b, const double c, const double d) {
  double solutions[3];
  size_t num_solutions;

  ComputeRealSolutionOfCubicEquation(b, c, d, solutions, &num_solutions);

  if (solutions[0] > solutions[1]) {
    std::swap(solutions[0], solutions[1]);
  }
  if (num_solutions == 3) {
    if (solutions[1] > solutions[2]) {
      std::swap(solutions[1], solutions[2]);
    }
    if (solutions[0] > solutions[1]) {
      std::swap(solutions[0], solutions[1]);
    }
  }

  printf("number of soluton: %lu\n", num_solutions);
  for (size_t i = 0; i < num_solutions; ++i) {
    printf("x_%lu = %f, f(x) = %e\n", i + 1, solutions[i],
           EvaluateCubicFunction(b, c, d, solutions[i]));
  }
  printf("\n");
}

int main(void) {
  size_t num_test_itr = 10000000;
  size_t num_bench_itr = 1000000;

  std::random_device seed_gen;
  std::default_random_engine engine(seed_gen());

  for (size_t itr_cnt = 0; itr_cnt < num_test_itr; ++itr_cnt) {
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    const double b = dist(engine);
    const double c = dist(engine);
    const double d = dist(engine);
    if (!Test(b, c, d)) {
      printf("Test failed\n");
      return 1;
    }
  }

  long sum_com_time = 0.0;
  for (size_t itr_cnt = 0; itr_cnt < num_bench_itr; ++itr_cnt) {
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    const double b = dist(engine);
    const double c = dist(engine);
    const double d = dist(engine);

    sum_com_time += Benchmark(b, c, d);
  }
  printf("\nbenchmark mean time: %f[ns]\n\n",
         double(sum_com_time) / double(num_bench_itr));
  {
    // (x - 1)(x - 2)(x - 3)
    const double b = -6;
    const double c = 11;
    const double d = -6;

    TryAndOutput(b, c, d);
  }
  {
    // (x - 2)(x - 4)^2
    const double b = -10;
    const double c = 32;
    const double d = -32;

    TryAndOutput(b, c, d);
  }
  {
    // (x - 4)^3
    const double b = -12;
    const double c = 48;
    const double d = -64;

    TryAndOutput(b, c, d);
  }
  {
    // (x - 1)(x^2 + 1)
    const double b = -1;
    const double c = 1;
    const double d = -1;

    TryAndOutput(b, c, d);
  }
  return 0;
}
