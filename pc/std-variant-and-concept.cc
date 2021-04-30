#include <iostream>
#include <variant>
#include <vector>

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
concept DebugDisplay = Debug<T>&& Display<T>;

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

template <DebugDisplay... Args>
using Base = std::variant<Args...>;

int main(void) {
  std::vector<Base<A, B>> hoge;

  hoge.emplace_back(A());
  hoge.emplace_back(B());

  for (const auto& v : hoge) {
    std::visit(
        [](auto&& u) {
          u.debug();
          u.display();
          if constexpr (Error<decltype(u)>) {
            u.error();
          }
        },
        v);
  }
}
