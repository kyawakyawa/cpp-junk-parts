/*
MIT License

Copyright (c) 2020 kyawakyawa

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

#include <iostream>
#include <utility>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#endif

#include <Eigen/Core>

#ifdef __clang__
#pragma clang diagnostic pop
#endif

// reference http://takashiijiri.com/study/miscs/QRfactorization.htm

template <typename Real, int N>
Eigen::Matrix<Real, N, N> HouseholderTransformation(
    const Eigen::Matrix<Real, N, 1>& x, const Eigen::Matrix<Real, N, 1>& y) {
  static_assert(N >= 1);
  const auto aux                    = x - y;
  const Eigen::Matrix<Real, N, 1> u = aux.normalized();

  return Eigen::Matrix<Real, N, N>::Identity() - 2 * u * u.transpose();
}

template <typename Real>
Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> HouseholderTransformation(
    const Eigen::Matrix<Real, Eigen::Dynamic, 1>& x,
    const Eigen::Matrix<Real, Eigen::Dynamic, 1>& y) {
  assert(x.rows() >= 1 && x.rows() == y.rows());
  const auto aux                                 = x - y;
  const Eigen::Matrix<Real, Eigen::Dynamic, 1> u = aux.normalized();

  return Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>::Identity(
             x.rows(), x.rows()) -
         2 * u * u.transpose();
}

template <typename Real, int N, int M>
std::pair<Eigen::Matrix<Real, N, N>, Eigen::Matrix<Real, N, M>>
QRDecompositionWithHouseholderTransformation(
    const Eigen::Matrix<Real, N, M>& input) {
  static_assert(M >= 1 && N >= M);
  Eigen::Matrix<Real, N, M> R = input;
  Eigen::Matrix<Real, N, N> Q = Eigen::Matrix<Real, N, N>::Identity();

  for (int k = 0; k < M; ++k) {
    const Real sign = (R(k, k) >= Real(0) ? Real(1) : Real(-1));
    const int Mmk   = M - k;

    const Eigen::Matrix<Real, Eigen::Dynamic, 1> x = R.block(k, k, Mmk, 1);
    // グラフィック用途(ローテーションとスケールの分離)の場合、
    // y(0,0)の符号は+で固定したほうがいい(この値がRの対角要素になる)
    Eigen::Matrix<Real, Eigen::Dynamic, 1> y =
        Eigen::Matrix<Real, Eigen::Dynamic, 1>::Zero(Mmk, 1);
    y(0, 0) = -sign * x.norm();

    const Eigen::Matrix<Real, Eigen::Dynamic, 1> u = (x - y).normalized();

    R.bottomRightCorner(Mmk, Mmk) -=
        Real(2) * u * (u.transpose() * R.bottomRightCorner(Mmk, Mmk));
    //////////////////// Rのみ欲しい場合はここまでの計算で良い
    ///////////////////////

    if (k > 0) {
      Q.topRightCorner(k, Mmk) -=
          Real(2) * (Q.topRightCorner(k, Mmk) * u) * u.transpose();
    }
    Q.bottomRightCorner(Mmk, Mmk) -=
        Real(2) * (Q.bottomRightCorner(Mmk, Mmk) * u) * u.transpose();
  }

  return std::make_pair(Q, R);
}

#ifdef USE_STACK_TRACE_LOGGER
#include <glog/logging.h>
#endif

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {
#ifdef USE_STACK_TRACE_LOGGER
  (void)argc;
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
#endif

  const Eigen::Vector3d x = {1, 2, 3};
  const Eigen::Vector3d y = {3, 2, 1};

  std::cout << "----------HouseholderTransformation----------" << std::endl;

  std::cout << " x = (" << x.transpose() << ")^T\n y = (" << y.transpose()
            << ")^T\n"
            << std::endl;

  const Eigen::Matrix3d H = HouseholderTransformation(x, y);
  std::cout << " H = \n" << H << "\n" << std::endl;

  const Eigen::Vector3d Hx = H * x;
  const Eigen::Vector3d Hy = H * y;

  std::cout << " (Hx)^T = (" << Hx.transpose() << ")\n" << std::endl;
  std::cout << " (Hy)^T = (" << Hy.transpose() << ")\n" << std::endl;

  std::cout << " H*H^T = \n" << H * H.transpose() << "\n" << std::endl;

  std::cout << "---------- ----------" << std::endl;

  std::cout
      << "----------QRDecompositionWithHouseholderTransformation----------"
      << std::endl;
  Eigen::Matrix3d A;

  A(0, 0) = 1.;
  A(0, 1) = 4.;
  A(0, 2) = 9.;

  A(1, 0) = 16.;
  A(1, 1) = 25.;
  A(1, 2) = 36.;

  A(2, 0) = 49.;
  A(2, 1) = 64.;
  A(2, 2) = 81.;

  const auto result = QRDecompositionWithHouseholderTransformation(A);

  std::cout << " Q = \n" << result.first << "\n" << std::endl;
  std::cout << " R = \n" << result.second << "\n" << std::endl;

  std::cout << " QR = \n" << result.first * result.second << "\n" << std::endl;

  std::cout << " QQ^T = \n"
            << result.first * result.first.transpose() << "\n"
            << std::endl;

  std::cout << "---------- ----------" << std::endl;
  return 0;
}
