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

class A {
public:
  A(void)  = default;
  ~A(void) = default;

  void Print(void) const {
    std::cout << msg << " in a member function" << std::endl;
  }

  friend struct PrintVisitor;

private:
  std::string msg = "class A";
};

class B {
public:
  B(void)  = default;
  ~B(void) = default;

  void Print(void) const {
    std::cout << msg << " in a member function" << std::endl;
  }

  friend struct PrintVisitor;

private:
  std::string msg = "class B";
};

struct PrintVisitor {
  void operator()(const A& a) {
    std::cout << a.msg << " in a functional object" << std::endl;
  }

  void operator()(const B& b) {
    std::cout << b.msg << " in a functional object" << std::endl;
  }
};

using C = std::variant<A, B>;

int main(void) {
  C c = A();

  std::visit([](const auto& v) { v.Print(); }, c);
  std::visit(PrintVisitor{}, c);

  c = B();

  std::visit([](const auto& v) { v.Print(); }, c);
  std::visit(PrintVisitor{}, c);

  return 0;
}
