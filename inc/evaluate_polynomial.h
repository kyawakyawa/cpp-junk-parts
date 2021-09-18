#include <math.h>

#include <algorithm>
#include <vector>

template <typename T = double>
T EvaluateQuadraticFunction(const T b, const T c, const T x) {
  // return (x + b) * x + c;
  std::vector<double> terms(3, 0.0);
  terms[0] = std::pow(x, 2);
  terms[1] = b * std::pow(x, 1);
  terms[2] = c;

  std::sort(terms.begin(), terms.end(),
            [&](const T l, const T r) { return abs(l) < abs(r); });

  double sum = T(0), remain = T(0);
  for (const auto t : terms) {
    const double diff = t + remain;
    const double tmp  = sum + diff;
    remain            = diff - (tmp - sum);
    sum               = tmp;
  }
  return sum;
}

template <typename T = double>
T EvaluateCubicFunction(const T b, const T c, const T d, const T x) {
  // return (((x + b) * x + c) * x + d);
  std::vector<double> terms(4, 0.0);
  terms[0] = std::pow(x, 3);
  terms[1] = b * std::pow(x, 2);
  terms[2] = c * std::pow(x, 1);
  terms[3] = d;

  std::sort(terms.begin(), terms.end(),
            [&](const T l, const T r) { return abs(l) < abs(r); });

  double sum = T(0), remain = T(0);
  for (const auto t : terms) {
    const double diff = t + remain;
    const double tmp  = sum + diff;
    remain            = diff - (tmp - sum);
    sum               = tmp;
  }
  return sum;
}

template <typename T = double>
T EvaluateQuaricFunction(const T b, const T c, const T d, const T e,
                         const T x) {
  // return (((x + b) * x + c) * x + d) * x + e;
  std::vector<double> terms(5, 0.0);
  terms[0] = std::pow(x, 4);
  terms[1] = b * std::pow(x, 3);
  terms[2] = c * std::pow(x, 2);
  terms[3] = d * std::pow(x, 1);
  terms[4] = e;

  std::sort(terms.begin(), terms.end(),
            [&](const T l, const T r) { return abs(l) < abs(r); });

  double sum = T(0), remain = T(0);
  for (const auto t : terms) {
    const double diff = t + remain;
    const double tmp  = sum + diff;
    remain            = diff - (tmp - sum);
    sum               = tmp;
  }
  return sum;
}

template <typename T = double>
T EvaluateDerivativeFunctionOfQuadraticFunction(const T b, const T x) {
  std::vector<double> terms(2, 0.0);
  terms[0] = T(2) * std::pow(x, 1);
  terms[1] = b;

  std::sort(terms.begin(), terms.end(),
            [&](const T l, const T r) { return abs(l) < abs(r); });

  double sum = T(0), remain = T(0);
  for (const auto t : terms) {
    const double diff = t + remain;
    const double tmp  = sum + diff;
    remain            = diff - (tmp - sum);
    sum               = tmp;
  }
  return sum;
}

template <typename T = double>
T EvaluateDerivativeFunctionOfCubicFunction(const T b, const T c, const T x) {
  std::vector<double> terms(3, 0.0);
  terms[0] = T(3) * std::pow(x, 2);
  terms[1] = T(2) * b * std::pow(x, 1);
  terms[2] = c;

  std::sort(terms.begin(), terms.end());

  double sum = T(0), remain = T(0);
  for (const auto t : terms) {
    const double diff = t + remain;
    const double tmp  = sum + diff;
    remain            = diff - (tmp - sum);
    sum               = tmp;
  }
  return sum;
}

template <typename T = double>
T EvaluateDerivativeFunctionOfQuaricFunction(const T b, const T c, const T d,
                                             const T x) {
  std::vector<double> terms(4, 0.0);
  terms[0] = T(4) * std::pow(x, 3);
  terms[1] = T(3) * b * std::pow(x, 2);
  terms[2] = T(2) * c * std::pow(x, 1);
  terms[3] = d;

  std::sort(terms.begin(), terms.end());

  double sum = T(0), remain = T(0);
  for (const auto t : terms) {
    const double diff = t + remain;
    const double tmp  = sum + diff;
    remain            = diff - (tmp - sum);
    sum               = tmp;
  }
  return sum;
}
