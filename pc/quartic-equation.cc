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
#include <time.h>

#include <algorithm>
#include <limits>
#include <random>

#define KAHAN_ADD(sum, add)    \
  diff   = add - remain;       \
  tmp    = sum + diff;         \
  remain = (tmp - sum) - diff; \
  sum    = tmp

#define KAHAN_ADD_LAST(sum, add) \
  diff = add - remain;           \
  tmp  = sum + diff;             \
  sum  = tmp

template <typename T = double>
void RefineCubicEquationSolution(const T b, const T c, const T d, T* x) {
  const T x2 = *x * *x;

  T numerator = T(0), denominator = T(0);
  T diff, tmp, remain;

  remain = T(0);
  KAHAN_ADD(numerator, d);
  KAHAN_ADD(numerator, c * *x);
  KAHAN_ADD(numerator, b * x2);
  KAHAN_ADD_LAST(numerator, *x * x2);

  remain = T(0);
  KAHAN_ADD(denominator, c);
  KAHAN_ADD(denominator, T(2) * b * *x);
  KAHAN_ADD_LAST(denominator, T(3) * x2);

  if (fabs(denominator) >
      std::numeric_limits<T>::epsilon()) {  // TODO Use the better way with
                                            // Ofast option
    *x -= numerator / denominator;
  }
}

template <typename T = double>
void ComputeRealSolutionOfQuadraticEquation(T b, const T c, T* solutions,
                                            size_t* num_solutions,
                                            const T D_thr = 1e-14) {
  // ref http://www.math.twcu.ac.jp/ogita/lec/na_basic.pdf
  // x^2  + b x + c =  0;
  //     |
  //     v
  // x^2 + 2b x + c =  0;
  //
  b *= T(0.5);

  const T D = b * b - c;

  if (D > D_thr) {
    *num_solutions        = 2;
    const T sqD           = sqrt(D);
    const int signbit_b   = std::signbit(b);  // if b >= 0 => 0 else b < 0 => 1
    const T sign_b        = signbit_b ? T(-1.0) : T(1.0);
    solutions[signbit_b]  = (-b - sign_b * sqD);
    solutions[!signbit_b] = c / solutions[signbit_b];
  } else if (D > -D_thr) {
    *num_solutions = 1;
    solutions[0]   = -b;
  } else {
    *num_solutions = 0;
  }
}

template <typename T = double>
void ComputeOneOfRealNonNegativeSolutionOfCubicEquation(T b, const T c,
                                                        const T d, T D_thr,
                                                        T* solutoin,
                                                        bool* find_root) {
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
    assert(p > 0.0);
    // 他の2つの実数解よりこの解のほうが大きい
    x_candidate = T(2) * sqrt(p) * cos(acos(q * sqrt(T(1) / p3)) / T(3)) - b;

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
  if (std::signbit(x_candidate)) {
    *find_root = false;
  } else {
    *solutoin  = x_candidate;
    *find_root = true;
  }
}

template <typename T = double>
void ComputeRealSolutionOfQuarticEquation(
    T b, const T c, const T d, const T e, T solutions[4], size_t* num_solutions,
    const T D_thr = 1e-14, const T biquadratic_equation_thr = 1e-6,
    const T multiple_solution_thr = 1e-14) {
  b *= T(0.25);

  const T b2 = b * b;

  const T p = c - T(6) * b2;
  const T q = d - T(2) * c * b + T(8) * b2 * b;
  const T r = e - d * b + c * b2 - T(3) * b2 * b2;

  if (abs(q) < biquadratic_equation_thr) {
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
    T u;
    bool find_root = 0;
    const T c1     = T(2) * p;
    const T c2     = p * p - T(4) * r;
    const T c3     = -q * q;

    ComputeOneOfRealNonNegativeSolutionOfCubicEquation(c1, c2, c3, D_thr, &u,
                                                       &find_root);
    if (!find_root) {
      *num_solutions = 0;
      return;
    }

    RefineCubicEquationSolution(c1, c2, c3, &u);

    const T sq_u      = sqrt(u);
    const T alpha     = (p + u) * T(0.5);
    const T beta      = (q * T(0.5)) / u;
    const T sq_u_beta = sq_u * beta;

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

template <typename T = double>
T EvaluateQuarticFunction(const T b, const T c, const T d, const T e,
                          const T x) {
  const T x_sq = x * x;

  T sum = T(0), remain = T(0);
  T diff, tmp;

  KAHAN_ADD(sum, x_sq * x_sq);
  KAHAN_ADD(sum, b * x * x_sq);
  KAHAN_ADD(sum, c * x_sq);
  KAHAN_ADD(sum, d * x);
  KAHAN_ADD_LAST(sum, e);

  return sum;
}

static bool Test(const double b, const double c, const double d,
                 const double e) {
  double solutions[4];
  size_t num_solutions;

  ComputeRealSolutionOfQuarticEquation(b, c, d, e, solutions, &num_solutions);

  if (num_solutions != 0 && num_solutions != 2 && num_solutions != 4) {
    return false;
  }

  const double b4 = b * 0.25;
  const double q  = d - 2.0 * c * b4 + 8.0 * b4 * b4 * b4;

  const double thr = 1e-12;

  // TODO newton法でエラーチェック

  for (size_t i = 0; i < num_solutions; ++i) {
    if (std::isfinite(solutions[i]) &&
        abs(EvaluateQuarticFunction(b, c, d, e, solutions[i])) > thr) {
      printf("error %e\n",
             abs(EvaluateQuarticFunction(b, c, d, e, solutions[i])));
      printf("x^4 + %f x^3 + %f x^2 + %f x + %f\n", b, c, d, e);
      printf("num solution %lu\n", num_solutions);
      for (size_t j = 0; j < num_solutions; ++j) {
        printf("x_%lu = %f\n", j, solutions[j]);
      }
      printf("q = %e\n", q);
      return false;
    }
  }
  return true;
}

static long Benchmark(const double b, const double c, const double d,
                      const double e) {
  double solutions[4];
  size_t num_solutions;

  struct timespec start_ts, end_ts;

  timespec_get(&start_ts, TIME_UTC);

  ComputeRealSolutionOfQuarticEquation(b, c, d, e, solutions, &num_solutions);

  timespec_get(&end_ts, TIME_UTC);
  return (long)end_ts.tv_sec * 1000000000 + end_ts.tv_nsec -
         start_ts.tv_sec * 1000000000 - start_ts.tv_nsec;
}

static bool TestAndOutput(const double b, const double c, const double d,
                          const double e, const size_t gt_num_solutions,
                          double gt_solutions[]) {
  double solutions[4];
  size_t num_solutions;

  ComputeRealSolutionOfQuarticEquation(
      b, c, d, e, solutions, &num_solutions, /*D_thr*/ 1e-14,
      /*biquadratic_equation_thr*/ 1e-14, /*multiple_solution_thr*/ 1e-14);

  std::sort(solutions, solutions + num_solutions);

  printf("number of soluton: %lu\n", num_solutions);
  for (size_t i = 0; i < num_solutions; ++i) {
    double y = EvaluateQuarticFunction(b, c, d, e, solutions[i]);
    printf("x_%lu = %f, f(x) = %.15e\n", i + 1, solutions[i], y);
  }
  if (num_solutions != gt_num_solutions) {
    printf("the number of solutions is wrong. (%lu (result) vs %lu (gt)\n",
           num_solutions, gt_num_solutions);
    return 0;
  }

  std::sort(gt_solutions, gt_solutions + gt_num_solutions);

  const double thr = 1e-15;
  for (size_t i = 0; i < num_solutions; ++i) {
    if (fabs(solutions[i] - gt_solutions[i]) > thr) {
      printf("the solution is wrong. (%.15f (result) vs %.15f (gt))\n",
             solutions[i], gt_solutions[i]);
      return 0;
    }
  }
  printf("\n");

  return 1;
}

int main(void) {
  size_t num_test_itr  = 10000000;
  size_t num_bench_itr = 10000000;

#if 0
   std::random_device seed_gen;
  std::default_random_engine engine(seed_gen());
#else
  std::default_random_engine engine(12345);
#endif

  for (size_t itr_cnt = 0; itr_cnt < num_test_itr; ++itr_cnt) {
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    const double b = dist(engine);
    const double c = dist(engine);
    const double d = dist(engine);
    const double e = dist(engine);
    if (!Test(b, c, d, e)) {
      printf("Test failed\n");
      return 1;
    }
  }
  printf("Pass Test\n");

  long sum_com_time = 0.0;
  for (size_t itr_cnt = 0; itr_cnt < num_bench_itr; ++itr_cnt) {
    std::uniform_real_distribution<double> dist(-10.0,10.0);
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

    const size_t gt_num_solutions = 4;
    double gt_solutions[]         = {-2, -1, 1, 2};

    if (!TestAndOutput(b, c, d, e, gt_num_solutions, gt_solutions)) {
      return 1;
    }
  }
  {
    // biquadratic equation cases
    const double b = 0;
    const double c = 1;
    const double d = 0;
    const double e = -6;

    const size_t gt_num_solutions = 2;
    double gt_solutions[]         = {-sqrt(2.0), sqrt(2.0)};

    if (!TestAndOutput(b, c, d, e, gt_num_solutions, gt_solutions)) {
      return 1;
    }
  }
  {
    // biquadratic equation cases
    const double b = 0;
    const double c = -1;
    const double d = 0;
    const double e = -6;

    const size_t gt_num_solutions = 2;
    double gt_solutions[]         = {-sqrt(3.0), sqrt(3.0)};

    if (!TestAndOutput(b, c, d, e, gt_num_solutions, gt_solutions)) {
      return 1;
    }
  }
  {
    // biquadratic equation cases
    const double b = 0;
    const double c = 5;
    const double d = 0;
    const double e = 4;

    if (!TestAndOutput(b, c, d, e, 0, nullptr)) {
      return 1;
    }
  }
  {
    // biquadratic equation cases
    // (x - 1)(x - 2)(x - 3)(x - 4)
    const double b = -10;
    const double c = 35;
    const double d = -50;
    const double e = 24;

    const size_t gt_num_solutions = 4;
    double gt_solutions[]         = {1, 2, 3, 4};

    if (!TestAndOutput(b, c, d, e, gt_num_solutions, gt_solutions)) {
      return 1;
    }
  }

  {
    // (x - 1)(x + 2)(x + 3)(x - 4)
    const double b = 0;
    const double c = -15;
    const double d = -10;
    const double e = +24;

    const size_t gt_num_solutions = 4;
    double gt_solutions[]         = {-3, -2, 1, 4};

    if (!TestAndOutput(b, c, d, e, gt_num_solutions, gt_solutions)) {
      return 1;
    }
  }

  {
    // (x - 1)(x + 2)(x + 3)^2
    const double b = 7;
    const double c = 13;
    const double d = -3;
    const double e = -18;

    const size_t gt_num_solutions = 3;
    double gt_solutions[]         = {-3, -2, 1};

    if (!TestAndOutput(b, c, d, e, gt_num_solutions, gt_solutions)) {
      return 1;
    }
  }

  {
    // (x + 1)^3(x - 3)
    const double b = 0;
    const double c = -6;
    const double d = -8;
    const double e = -3;

    const size_t gt_num_solutions = 2;
    double gt_solutions[]         = {-1, 3};

    if (!TestAndOutput(b, c, d, e, gt_num_solutions, gt_solutions)) {
      return 1;
    }
  }

  {
    // (x + 1)^2(x - 3)^2
    const double b = -4;
    const double c = -2;
    const double d = 12;
    const double e = 9;

    const size_t gt_num_solutions = 2;
    double gt_solutions[]         = {-1, 3};

    if (!TestAndOutput(b, c, d, e, gt_num_solutions, gt_solutions)) {
      return 1;
    }
  }

  {
    // (x - 5)^4
    const double b = -20;
    const double c = 150;
    const double d = -500;
    const double e = 625;

    const size_t gt_num_solutions = 1;
    double gt_solutions[]         = {5};

    if (!TestAndOutput(b, c, d, e, gt_num_solutions, gt_solutions)) {
      return 1;
    }
  }

  {
    // (x - 2)^2(x^2 + 4x + 5)
    const double b = 0;
    const double c = -7;
    const double d = -4;
    const double e = 20;

    const size_t gt_num_solutions = 1;
    double gt_solutions[]         = {2};

    if (!TestAndOutput(b, c, d, e, gt_num_solutions, gt_solutions)) {
      return 1;
    }
  }

  {
    // (x^2 -2x + 10)(x^2 + 4x + 5)
    const double b = 2;
    const double c = 7;
    const double d = 30;
    const double e = 50;

    if (!TestAndOutput(b, c, d, e, 0, nullptr)) {
      return 1;
    }
  }

  return 0;
}
