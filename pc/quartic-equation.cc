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
#include <optional>
#include <random>
#include <vector>

template <typename T = double>
T EvaluateQuadraticFunction(const T b, const T c, const T x) {
  return (x + b) * x + c;
}

template <typename T = double>
T EvaluateCubicFunction(const T b, const T c, const T d, const T x) {
  return ((x + b) * x + c) * x + d;
}

template <typename T = double>
T EvaluateQuaricFunction(const T b, const T c, const T d, const T e,
                         const T x) {
  return (((x + b) * x + c) * x + d) * x + e;
}

template <typename T = double>
inline void ComputeRealSolutionOfQuadraticEquation(const T b, const T c,
                                                   T* solutions,
                                                   size_t* num_solutions,
                                                   const T D_thr = 1e-14) {
  const T D = b * b - T(4.0) * c;

  if (D > D_thr) {
    *num_solutions = 2;
    solutions[0]   = (-b - sqrt(D)) * T(0.5);
    solutions[1]   = (-b + sqrt(D)) * T(0.5);
  } else if (D > -D_thr) {
    *num_solutions = 1;
    solutions[0]   = -b * T(0.5);
  } else {
    *num_solutions = 0;
  }
}

template <typename T = double>
inline std::optional<T> ComputeOneOfRealNonNegativeSolutionOfCubicEquation(
    T b, const T c, const T d, T D_thr = 1e-14) {
  // x^3 + b x^2  + c x + d =  0;
  //             |
  //             v
  // x^3 + 3b x^2  + c x + d =  0;

  b *= (T(1) / T(3));

  const T b2 = b * b;

  const T p = b2 - c * (T(1) / T(3));
  const T q = (b * c - T(2.0) * b * b2 - d) * T(0.5);

  // y^3 - 3 p y - 2 q = 0
  //       where y = x + b

  const T p3 = p * p * p;
  const T D  = q * q - p3;

  T x_candidate = T(-1.0);
  if (D < -D_thr) {
    assert(p > T(0.0));
    const T sq_p_2        = sqrt(p /* FIXME */) * T(2);
    const T theta         = acos(q * std::pow(p, -T(3) / T(2)));
    constexpr T kPi2      = T(2.0 * 3.141592653589793238462643383279);
    constexpr T kPi4      = T(4.0 * 3.141592653589793238462643383279);
    constexpr T three_div = T(1) / T(3);

    // TODO SIMD
    const T x0 = sq_p_2 * cos(theta * three_div) - b;
    const T x1 = sq_p_2 * cos((theta + kPi2) * three_div) - b;
    const T x2 = sq_p_2 * cos((theta + kPi4) * three_div) - b;

    x_candidate = std::max({x0, x1, x2});

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

    x_candidate = y1 - b;

    if (std::signbit(x_candidate) && D < D_thr) {
      x_candidate = y1 * T(-0.5) - b;
    }
  }
  return std::signbit(x_candidate) ? std::nullopt
                                   : std::optional<T>(x_candidate);
}

template <typename T = double>
inline void ComputeRealSolutionOfQuarticEquation(
    T b, const T c, const T d, const T e, T solutions[4], size_t* num_solutions,
    const T D_thr = 1e-14, const T biquadratic_equation_thr = 1e-6,
    const T multiple_solution_thr = 1e-14) {
  b *= T(0.25);

  const T b2 = b * b;

  const T p = c - T(6) * b2;
  const T q = d - T(2) * c * b + T(8) * b2 * b;
  const T r = e - d * b + c * b2 - T(3) * b2 * b2;

  if (abs(q) < biquadratic_equation_thr) {
    // printf("aaa\n");
    T quadratic_equation_solutions[2];
    size_t num_quadratic_equation_solutions;
    ComputeRealSolutionOfQuadraticEquation(p, r, quadratic_equation_solutions,
                                           &num_quadratic_equation_solutions,
                                           D_thr);
    *num_solutions = 0;
    if (num_quadratic_equation_solutions >= 1) {
      bool zero_flag = false;  // To deal with floating point arithmetic errors.
      if (quadratic_equation_solutions[0] > multiple_solution_thr) {
        const T sq                    = sqrt(quadratic_equation_solutions[0]);
        solutions[(*num_solutions)++] = -sq;
        solutions[(*num_solutions)++] = sq;
      } else if (quadratic_equation_solutions[0] > -multiple_solution_thr) {
        solutions[(*num_solutions)++] = T(0);
        zero_flag                     = true;
      }

      if (num_quadratic_equation_solutions == 2) {
        if (quadratic_equation_solutions[1] > multiple_solution_thr) {
          const T sq                    = sqrt(quadratic_equation_solutions[1]);
          solutions[(*num_solutions)++] = -sq;
          solutions[(*num_solutions)++] = sq;
        } else if (!zero_flag &&
                   quadratic_equation_solutions[1] > -multiple_solution_thr) {
          solutions[(*num_solutions)++] = T(0);
        }
      }
    }
  } else {
    // printf("bbb\n");
    const std::optional<T> u_op =
        ComputeOneOfRealNonNegativeSolutionOfCubicEquation(
            T(2) * p, (p * p - T(4) * r), -q * q, D_thr);

    if (!u_op) {
      *num_solutions = 0;
      return;
    }
    const T u = u_op.value();
    // printf("u %le\n", u);
    // printf("cubic equation error %e\n",
    //       EvaluateCubicFunction(T(2) * p, (p * p - T(4) * r), -q * q, u));

    const T sq_u      = sqrt(u);
    const T alpha     = (p + u) * T(0.5);
    const T beta      = (q * T(0.5)) / u;
    const T sq_u_beta = sq_u * beta;
    // printf("sq_u_beta %le\n", sq_u_beta);

    size_t num_quadratic_equation_solutions0;

    ComputeRealSolutionOfQuadraticEquation(sq_u, alpha - sq_u_beta, solutions,
                                           &num_quadratic_equation_solutions0,
                                           D_thr);

    T quadratic_equation_solutions1[2];
    size_t num_quadratic_equation_solutions1;

    ComputeRealSolutionOfQuadraticEquation(
        -sq_u, alpha + sq_u_beta, quadratic_equation_solutions1,
        &num_quadratic_equation_solutions1, D_thr);

    char duplicate_with_solution0[2] = {false, false};
    for (size_t i = 0; i < num_quadratic_equation_solutions0; ++i) {
      for (size_t j = 0; j < num_quadratic_equation_solutions1; ++j) {
        if (abs(solutions[i] - quadratic_equation_solutions1[j]) <
            multiple_solution_thr) {
          duplicate_with_solution0[j] = true;
        }
      }
    }

    *num_solutions = num_quadratic_equation_solutions0;

    if (num_quadratic_equation_solutions1 >= 1 &&
        !duplicate_with_solution0[0]) {
      solutions[(*num_solutions)++] = quadratic_equation_solutions1[0];
    }
    if (num_quadratic_equation_solutions1 == 2 &&
        !duplicate_with_solution0[1]) {
      solutions[(*num_solutions)++] = quadratic_equation_solutions1[1];
    }
  }

  for (size_t i = 0; i < *num_solutions; ++i) {
    // TODO SIMD
    solutions[i] -= b;
  }
}

static bool Test(const double b, const double c, const double d,
                 const double e) {
  double solutions[4];
  size_t num_solutions;

  ComputeRealSolutionOfQuarticEquation(b, c, d, e, solutions, &num_solutions);

  if (num_solutions != 0 && num_solutions != 2 && num_solutions != 4) {
    return false;
  }

  // 現状 qが小さいと精度が悪くなる
  const double b4 = b * 0.25;
  const double q  = d - 2.0 * c * b4 + 8.0 * b4 * b4 * b4;

  // qが小さいときは甘くする
  const double thr = (abs(q) < 1e-2 ? 5e-2 : 1e-7);

  // TODO newton法でエラーチェック

  for (size_t i = 0; i < num_solutions; ++i) {
    if (std::isfinite(solutions[i]) &&
        abs(EvaluateQuaricFunction(b, c, d, e, solutions[i])) > thr) {
      printf("error %e\n",
             abs(EvaluateQuaricFunction(b, c, d, e, solutions[i])));
      printf("x^4 + %f x^3 + %f x^2 + %f x + %f\n", b, c, d, e);
      printf("num solution %lu\n", num_solutions);
      for (size_t j = 0; j < num_solutions; ++j) {
        printf("x_%lu = %f\n", j, solutions[j]);
      }
      printf("q = %e\n", d - 2.0 * c * b4 + 8.0 * b4 * b4 * b4);
      return false;
    }
  }
  return true;
}

static long Benchmark(const double b, const double c, const double d,
                      const double e) {
  std::chrono::system_clock::time_point start, end;

  double solutions[4];
  size_t num_solutions;

  start = std::chrono::high_resolution_clock::now();
  ComputeRealSolutionOfQuarticEquation(b, c, d, e, solutions, &num_solutions);
  end = std::chrono::high_resolution_clock::now();

  const long com_time =
      std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

  return com_time;
}

static void TryAndOutput(const double b, const double c, const double d,
                         const double e) {
  double solutions[4];
  size_t num_solutions;

  ComputeRealSolutionOfQuarticEquation(b, c, d, e, solutions, &num_solutions);

  std::sort(solutions, solutions + num_solutions);

  printf("number of soluton: %lu\n", num_solutions);
  for (size_t i = 0; i < num_solutions; ++i) {
    printf("x_%lu = %f, f(x) = %e\n", i + 1, solutions[i],
           EvaluateQuaricFunction(b, c, d, e, solutions[i]));
  }
  printf("\n");
}

int main(void) {
  size_t num_test_itr  = 10000000;
  size_t num_bench_itr = 1000000;

#if 0
   std::random_device seed_gen;
  std::default_random_engine engine(seed_gen());
#else
  std::default_random_engine engine(12345);
#endif

  for (size_t itr_cnt = 0; itr_cnt < num_test_itr; ++itr_cnt) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    const double b = dist(engine);
    const double c = dist(engine);
    const double d = dist(engine);
    const double e = dist(engine);
    if (!Test(b, c, d, e)) {
      printf("Test failed\n");
      return 1;
    }
  }

  long sum_com_time = 0.0;
  for (size_t itr_cnt = 0; itr_cnt < num_bench_itr; ++itr_cnt) {
    std::uniform_real_distribution<double> dist(-5.0, 5.0);
    const double b = dist(engine);
    const double c = dist(engine);
    const double d = dist(engine);
    const double e = dist(engine);

    sum_com_time += Benchmark(b, c, d, e);
  }
  printf("\nbenchmark mean time: %f[ns]\n\n",
         double(sum_com_time) / double(num_bench_itr));
  {
    // biquadratic equation cases
    const double b = 0;
    const double c = -5;
    const double d = 0;
    const double e = 4;

    TryAndOutput(b, c, d, e);
  }
  {
    // biquadratic equation cases
    const double b = 0;
    const double c = 1;
    const double d = 0;
    const double e = -6;

    TryAndOutput(b, c, d, e);
  }
  {
    // biquadratic equation cases
    const double b = 0;
    const double c = -1;
    const double d = 0;
    const double e = -6;

    TryAndOutput(b, c, d, e);
  }
  {
    // biquadratic equation cases
    const double b = 0;
    const double c = 5;
    const double d = 0;
    const double e = 4;

    TryAndOutput(b, c, d, e);
  }
  {
    // biquadratic equation cases
    // (x - 1)(x - 2)(x - 3)(x - 4)
    const double b = -10;
    const double c = 35;
    const double d = -50;
    const double e = 24;

    TryAndOutput(b, c, d, e);
  }

  {
    // (x - 1)(x + 2)(x + 3)(x - 4)
    const double b = 0;
    const double c = -15;
    const double d = -10;
    const double e = +24;

    TryAndOutput(b, c, d, e);
  }

  {
    // (x - 1)(x + 2)(x + 3)(x + 3)
    const double b = 7;
    const double c = 13;
    const double d = -3;
    const double e = -18;

    TryAndOutput(b, c, d, e);
  }

  {
    // (x + 1)^3(x - 3)
    const double b = 0;
    const double c = -6;
    const double d = -8;
    const double e = -3;

    TryAndOutput(b, c, d, e);
  }

  {
    // (x + 1)^2(x - 3)^2
    const double b = -4;
    const double c = -2;
    const double d = 12;
    const double e = 9;

    TryAndOutput(b, c, d, e);
  }

  {
    // (x - 5)^4
    const double b = -20;
    const double c = 150;
    const double d = -500;
    const double e = 625;

    TryAndOutput(b, c, d, e);
  }

  {
    // (x - 2)^2(x^2 + 4x + 5)
    const double b = 0;
    const double c = -7;
    const double d = -4;
    const double e = 20;

    TryAndOutput(b, c, d, e);
  }

  {
    // (x^2 -2x + 10)(x^2 + 4x + 5)
    const double b = 2;
    const double c = 7;
    const double d = 30;
    const double e = 50;

    TryAndOutput(b, c, d, e);
  }

  return 0;
}
