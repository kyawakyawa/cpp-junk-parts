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

#include <iostream>
#include <variant>
#include <vector>

#include "macros.h"

template <class T>
concept Debug = requires(T& x) {
  x.debug();
};

template <class T>
concept Display = requires(T& x) {
  x.display();
};

template <class T>
concept Error = requires(T& x) {
  x.error();
};

template <class T>
concept DebugDisplay = Debug<T> && Display<T>;

class A {
 public:
  void debug() const { std::cout << "debug A" << std::endl; }

  void display() const { std::cout << "display A" << std::endl; }

  void error() const { std::cerr << "error A" << std::endl; }
};

class B {
 public:
  void debug() const { std::cout << "debug B" << std::endl; }

  void display() const { std::cout << "display B" << std::endl; }
};

// This is also easy to use, but may cause some inconvenience when making
// forward declarations.
//
// template <DebugDisplay... Args>
// using VariantOfAB = std::variant<Args...>;
// using Base = VariantOfAB<A, B>;

using Base = std::variant<A, B /*, int */ /* error */>;

CHECK_VARIANT_CONSTRAINTS(Base, DebugDisplay);

// error
// CHECK_VARIANT_CONSTRAINTS(Base, Error);

VARIANT_CONSTRAINTS_CHECKER(DebugDisplay)
VARIANT_CONSTRAINTS_CHECKER(Error)

[[maybe_unused]] static void CheckVariantSatisfiesConcepts(void) {
  CheckVariantSatisfiesDebugDisplay<Base>();

  // error
  // CheckVariantSatisfiesError<Base>();
}

int main(void) {
  std::vector<Base> hoge;

  hoge.emplace_back(A());
  hoge.emplace_back(B());

  for (const auto& v : hoge) {
    std::visit(
        []<DebugDisplay T>(T&& u) {
          u.debug();
          u.display();
          if constexpr (Error<decltype(u)>) {
            u.error();
          }
        },
        v);
  }
}
