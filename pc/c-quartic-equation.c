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
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define KAHAN_ADD(sum, add)    \
  diff = add - remain;         \
  tmp = sum + diff;            \
  remain = (tmp - sum) - diff; \
  sum = tmp

#define KAHAN_ADD_LAST(sum, add) \
  diff = add - remain;           \
  tmp = sum + diff;              \
  sum = tmp

// NOLINTNEXTLINE
static double evaluate_quadratic_function(const double b, const double c,
                                          const double x) {
  double sum = 0.0;
  double remain = 0.0;
  double diff;
  double tmp;

  KAHAN_ADD(sum, x * x);
  KAHAN_ADD(sum, b * x);
  KAHAN_ADD_LAST(sum, c);

  return sum;
}

// NOLINTNEXTLINE
static double evaluate_cubic_function(const double b, const double c,
                                      const double d, const double x) {
  const double x_sq = x * x;

  double sum = 0.0;
  double remain = 0.0;
  double diff;
  double tmp;

  KAHAN_ADD(sum, x * x_sq);
  KAHAN_ADD(sum, b * x_sq);
  KAHAN_ADD(sum, c * x);
  KAHAN_ADD_LAST(sum, d);

  return sum;
}

static double evaluate_quartic_function(const double b, const double c,
                                        const double d, const double e,
                                        const double x) {
  const double x_sq = x * x;

  double sum = 0.0;
  double remain = 0.0;
  double diff;
  double tmp;

  KAHAN_ADD(sum, x_sq * x_sq);
  KAHAN_ADD(sum, b * x * x_sq);
  KAHAN_ADD(sum, c * x_sq);
  KAHAN_ADD(sum, d * x);
  KAHAN_ADD_LAST(sum, e);

  return sum;
}

// NOLINTNEXTLINE
static double evaluate_derivative_function_of_quadratic_function(
    const double b, const double x) {
  return 2.0 * x + b;
}

// NOLINTNEXTLINE
static double evaluate_derivative_function_of_cubic_function(const double b,
                                                             const double c,
                                                             const double x) {
  double sum = 0.0;
  double remain = 0.0;
  double diff;
  double tmp;

  KAHAN_ADD(sum, 3.0 * x * x);
  KAHAN_ADD(sum, 2.0 * b * x);
  KAHAN_ADD_LAST(sum, c);

  return sum;
}

// NOLINTNEXTLINE
static double evaluate_derivative_function_of_quartic_function(const double b,
                                                               const double c,
                                                               const double d,
                                                               const double x) {
  const double x_sq = x * x;

  double sum = 0.0;
  double remain = 0.0;
  double diff;
  double tmp;

  KAHAN_ADD(sum, 4.0 * x * x_sq);
  KAHAN_ADD(sum, 3.0 * b * x_sq);
  KAHAN_ADD(sum, 2.0 * c * x);
  KAHAN_ADD_LAST(sum, d);

  return sum;
}

static void refine_cubic_equation_solution(const double b, const double c,
                                           const double d, double* x) {
  const double x2 = *x * *x;

  double numerator = 0;
  double denominator = 0;
  double diff;
  double tmp;
  double remain;

  remain = 0;
  KAHAN_ADD(numerator, d);
  KAHAN_ADD(numerator, c * *x);
  KAHAN_ADD(numerator, b * x2);
  KAHAN_ADD_LAST(numerator, *x * x2);

  remain = 0;
  KAHAN_ADD(denominator, c);
  KAHAN_ADD(denominator, 2.0 * b * *x);
  KAHAN_ADD_LAST(denominator, 3.0 * x2);

  if (fabs(denominator) >
      DBL_EPSILON) {  // TODO Use the better way with Ofast option
    *x -= numerator / denominator;
  }
}

static void compute_real_solution_of_quadratic_equation(double b,
                                                        const double c,
                                                        double* solutions,
                                                        size_t* num_solutions,
                                                        const double D_thr) {
  // ref http://www.math.twcu.ac.jp/ogita/lec/na_basic.pdf
  // x^2  + b x + c =  0;
  //     |
  //     v
  // x^2 + 2b x + c =  0;
  //
  b *= 0.5;

  const double D = b * b - c;
  if (D > D_thr) {
    *num_solutions = 2;
    const double sqD = sqrt(D);
    const int signbit_b =
        (signbit(b) == 0 ? 0 : 1);  // if b >= 0 => 0 else b < 0 => 1
    const double sign_b = signbit_b ? -1.0 : 1.0;
    solutions[signbit_b] = (-b - sign_b * sqD);
    solutions[!signbit_b] = c / solutions[signbit_b];
  } else if (D > -D_thr) {
    *num_solutions = 1;
    solutions[0] = -b;
  } else {
    *num_solutions = 0;
  }
}

static void compute_one_of_real_non_negative_solution_of_cubic_equation(
    double b, const double c, const double d, double D_thr, double* solution,
    char* find_root) {
  // x^3 + b x^2  + c x + d =  0;
  //             |
  //             v
  // x^3 + 3b x^2  + c x + d =  0;

  b *= 1.0 / 3.0;

  const double b2 = b * b;

  const double p = b2 - c * (1.0 / 3.0);
  const double q = (b * c - 2.0 * b * b2 - d) * 0.5;

  // y^3 - 3 p y - 2 q = 0
  //       where y = x + b

  const double p3 = p * p * p;
  const double D = q * q - p3;

  double x_candidate = -1.0;
  if (D < -D_thr) {
    assert(p > 0.0);
    // 他の2つの実数解よりこの解のほうが大きい
    x_candidate = 2.0 * sqrt(p) * cos(acos(q * sqrt(1.0 / p3)) / 3.0) - b;

  } else {
    const double sqD = (D < D_thr ? 0.0 : sqrt(D));
    const double alpha_3 =
        (!signbit(q) ? cbrt(q + sqD)
                     : cbrt(q - sqD));  // alpha3 の絶対値が大きくなるように
    const double beta_3 =
        (D < D_thr ? alpha_3
                   : p / alpha_3);  // 重解のときはゼロ割に近くなるので場合分け

    const double y1 = alpha_3 + beta_3;

    x_candidate = y1 - b;

    if (signbit(x_candidate) && D < D_thr) {
      x_candidate = y1 * -0.5 - b;
    }
  }
  if (signbit(x_candidate)) {
    *find_root = 0;
  } else {
    *solution = x_candidate;
    *find_root = 1;
  }
}

// NOLINTNEXTLINE
static void compute_real_solution_of_quartic_equation(
    double b, const double c, const double d, const double e,
    double solutions[4], size_t* num_solutions, const double D_thr,
    const double biquadratic_equation_thr, const double multiple_solution_thr) {
  b *= 0.25;

  const double b2 = b * b;

  const double p = c - 6.0 * b2;
  const double q = d - 2.0 * c * b + 8.0 * b2 * b;
  const double r = e - d * b + c * b2 - 3.0 * b2 * b2;

  if (fabs(q) < biquadratic_equation_thr) {
    double quadratic_equation_solutions[2];
    size_t num_quadratic_equation_solutions;
    compute_real_solution_of_quadratic_equation(
        p, r, quadratic_equation_solutions, &num_quadratic_equation_solutions,
        D_thr);
    *num_solutions = 0;
    if (num_quadratic_equation_solutions >= 1) {
      char zero_flag = 0;  // To deal with floating point arithmetic errors.
      if (quadratic_equation_solutions[0] > multiple_solution_thr) {
        const double sq = sqrt(quadratic_equation_solutions[0]);
        solutions[(*num_solutions)++] = -sq;
        solutions[(*num_solutions)++] = sq;
      } else if (quadratic_equation_solutions[0] > -multiple_solution_thr) {
        solutions[(*num_solutions)++] = 0.0;
        zero_flag = 1;
      }

      if (num_quadratic_equation_solutions == 2) {
        if (quadratic_equation_solutions[1] > multiple_solution_thr) {
          const double sq = sqrt(quadratic_equation_solutions[1]);
          solutions[(*num_solutions)++] = -sq;
          solutions[(*num_solutions)++] = sq;
        } else if (!zero_flag &&
                   quadratic_equation_solutions[1] > -multiple_solution_thr) {
          solutions[(*num_solutions)++] = 0.0;
        }
      }
    }
  } else {
    double u;
    char find_root = 0;
    const double c1 = 2.0 * p;
    const double c2 = p * p - 4.0 * r;
    const double c3 = -q * q;

    compute_one_of_real_non_negative_solution_of_cubic_equation(
        c1, c2, c3, D_thr, &u, &find_root);

    if (!find_root) {
      *num_solutions = 0;
      return;
    }

    refine_cubic_equation_solution(c1, c2, c3, &u);

    const double sq_u = sqrt(u);
    const double alpha = (p + u) * 0.5;
    const double beta = (q * 0.5) / u;
    const double sq_u_beta = sq_u * beta;

    size_t num_quadratic_equation_solutions0;

    compute_real_solution_of_quadratic_equation(
        sq_u, alpha - sq_u_beta, solutions, &num_quadratic_equation_solutions0,
        D_thr);

    double quadratic_equation_solutions1[2];
    size_t num_quadratic_equation_solutions1;

    compute_real_solution_of_quadratic_equation(
        -sq_u, alpha + sq_u_beta, quadratic_equation_solutions1,
        &num_quadratic_equation_solutions1, D_thr);

    char duplicate_with_solution0[2] = {0, 0};
    for (size_t i = 0; i < num_quadratic_equation_solutions0; ++i) {
      for (size_t j = 0; j < num_quadratic_equation_solutions1; ++j) {
        if (fabs(solutions[i] - quadratic_equation_solutions1[j]) <
            multiple_solution_thr) {
          duplicate_with_solution0[j] = 1;
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

static _Bool test(const double b, const double c, const double d,
                  const double e) {
  double solutions[4];
  size_t num_solutions;

  compute_real_solution_of_quartic_equation(
      b, c, d, e, solutions, &num_solutions, 1e-14, 1e-14, 1e-14);

  if (num_solutions != 0 && num_solutions != 2 && num_solutions != 4) {
    return 0;
  }

  const double b4 = b * 0.25;
  const double q = d - 2.0 * c * b4 + 8.0 * b4 * b4 * b4;

  const double thr = 1.0e-12;

  for (size_t i = 0; i < num_solutions; ++i) {
    if (isfinite(solutions[i]) &&
        fabs(evaluate_quartic_function(b, c, d, e, solutions[i])) > thr) {
      printf("error %.15e, when x_%lu is %.15f\n",
             fabs(evaluate_quartic_function(b, c, d, e, solutions[i])), i,
             solutions[i]);
      printf("x^4 + %f x^3 + %f x^2 + %f x + %f\n", b, c, d, e);
      printf("num solution %lu\n", num_solutions);
      for (size_t j = 0; j < num_solutions; ++j) {
        printf("x_%lu = %f\n", j, solutions[j]);
      }
      printf("q = %e\n", q);
      return 0;
    }
  }
  return 1;
}

static int compare_double(const void* a, const void* b) {
  const double da = *(const double*)a;
  const double db = *(const double*)b;
  if (da - db < -1e15) {
    return -1;
  }
  if (da - db > 1e-15) {
    return 1;
  }
  return 0;
}

static _Bool test_and_output(const double b, const double c, const double d,
                             const double e, const size_t gt_num_solutions,
                             double gt_solutions[]) {
  double solutions[4];
  size_t num_solutions;

  compute_real_solution_of_quartic_equation(
      b, c, d, e, solutions, &num_solutions, /*D_thr*/ 1e-14,
      /*biquadratic_equation_thr*/ 1e-14, /*multiple_solution_thr*/ 1e-14);

  qsort(solutions, num_solutions, sizeof(double), compare_double);

  printf("number of soluton: %lu\n", num_solutions);
  for (size_t i = 0; i < num_solutions; ++i) {
    double y = evaluate_quartic_function(b, c, d, e, solutions[i]);
    printf("x_%lu = %f, f(x) = %.15e\n", i + 1, solutions[i], y);
  }
  if (num_solutions != gt_num_solutions) {
    printf("the number of solutions is wrong. (%lu (result) vs %lu (gt)\n",
           num_solutions, gt_num_solutions);
    return 0;
  }

  qsort(gt_solutions, gt_num_solutions, sizeof(double), compare_double);

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

static uint64_t benchmark(const double b, const double c, const double d,
                          const double e) {
  double solutions[4];
  size_t num_solutions;

  struct timespec start_ts;
  struct timespec end_ts;

  timespec_get(&start_ts, TIME_UTC);

  compute_real_solution_of_quartic_equation(
      b, c, d, e, solutions, &num_solutions, 1e-14, 1e-14, 1e-14);

  timespec_get(&end_ts, TIME_UTC);
  return (uint64_t)end_ts.tv_sec * 1000000000 + (uint64_t)end_ts.tv_nsec -
         (uint64_t)start_ts.tv_sec * 1000000000 - (uint64_t)start_ts.tv_nsec;
}

static uint64_t xor128(void) {
  static uint64_t x = 123456789;
  static uint64_t y = 362436069;
  static uint64_t z = 521288629;
  static uint64_t w = 88675123;
  uint64_t t;
  t = (x ^ (x << 11));
  x = y;
  y = z;
  z = w;
  return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
}

int main(void) {
  size_t num_test_itr = 10000000;
  size_t num_bench_itr = 10000000;

  for (size_t itr_cnt = 0; itr_cnt < num_test_itr; ++itr_cnt) {
    const double b = ((double)(xor128()) / (double)(UINT64_MAX - 1) - 0.5) * 10;
    const double c = ((double)(xor128()) / (double)(UINT64_MAX - 1) - 0.5) * 10;
    const double d = ((double)(xor128()) / (double)(UINT64_MAX - 1) - 0.5) * 10;
    const double e = ((double)(xor128()) / (double)(UINT64_MAX - 1) - 0.5) * 10;

    if (!test(b, c, d, e)) {
      printf("Test failed\n");
      return 1;
    }
  }
  printf("Pass Test\n");

  uint64_t sum_com_time = 0.0;
  for (size_t itr_cnt = 0; itr_cnt < num_bench_itr; ++itr_cnt) {
    const double b = (double)(xor128()) / (double)(UINT64_MAX - 1) * 10;
    const double c = (double)(xor128()) / (double)(UINT64_MAX - 1) * 10;
    const double d = (double)(xor128()) / (double)(UINT64_MAX - 1) * 10;
    const double e = (double)(xor128()) / (double)(UINT64_MAX - 1) * 10;

    sum_com_time += benchmark(b, c, d, e);
  }
  printf("bench %f ns\n", (double)sum_com_time / (double)num_bench_itr);

  {
    // biquadratic equation cases
    const double b = 0;
    const double c = -5;
    const double d = 0;
    const double e = 4;

    const size_t gt_num_solutions = 4;
    double gt_solutions[] = {-2, -1, 1, 2};

    if (!test_and_output(b, c, d, e, gt_num_solutions, gt_solutions)) {
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
    double gt_solutions[] = {-sqrt(2.0), sqrt(2.0)};

    if (!test_and_output(b, c, d, e, gt_num_solutions, gt_solutions)) {
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
    double gt_solutions[] = {-sqrt(3.0), sqrt(3.0)};

    if (!test_and_output(b, c, d, e, gt_num_solutions, gt_solutions)) {
      return 1;
    }
  }
  {
    // biquadratic equation cases
    const double b = 0;
    const double c = 5;
    const double d = 0;
    const double e = 4;

    if (!test_and_output(b, c, d, e, 0, NULL)) {
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
    double gt_solutions[] = {1, 2, 3, 4};

    if (!test_and_output(b, c, d, e, gt_num_solutions, gt_solutions)) {
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
    double gt_solutions[] = {-3, -2, 1, 4};

    if (!test_and_output(b, c, d, e, gt_num_solutions, gt_solutions)) {
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
    double gt_solutions[] = {-3, -2, 1};

    if (!test_and_output(b, c, d, e, gt_num_solutions, gt_solutions)) {
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
    double gt_solutions[] = {-1, 3};

    if (!test_and_output(b, c, d, e, gt_num_solutions, gt_solutions)) {
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
    double gt_solutions[] = {-1, 3};

    if (!test_and_output(b, c, d, e, gt_num_solutions, gt_solutions)) {
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
    double gt_solutions[] = {5};

    if (!test_and_output(b, c, d, e, gt_num_solutions, gt_solutions)) {
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
    double gt_solutions[] = {2};

    if (!test_and_output(b, c, d, e, gt_num_solutions, gt_solutions)) {
      return 1;
    }
  }

  {
    // (x^2 -2x + 10)(x^2 + 4x + 5)
    const double b = 2;
    const double c = 7;
    const double d = 30;
    const double e = 50;

    if (!test_and_output(b, c, d, e, 0, NULL)) {
      return 1;
    }
  }

  return 0;
}
