#include <execution>
#include <numeric>
#include <random>
#include <vector>

#ifdef USE_STACK_TRACE_LOGGER
#include <glog/logging.h>
#endif

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
#ifdef USE_STACK_TRACE_LOGGER
  (void)argc;
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
#endif

  std::vector<int> vec(1000000);
  std::iota(vec.begin(), vec.end(), 1);
  std::random_device rnd;
  std::shuffle(vec.begin(), vec.end(), std::mt19937(rnd()));
  for (int i = 0; i < 1000; i++) {
    std::vector<int> _vec(vec.size());
    std::copy(std::execution::par_unseq, vec.begin(), vec.end(), _vec.begin());
    std::sort(std::execution::par_unseq, _vec.begin(), _vec.end());
  }
  return 0;
}
