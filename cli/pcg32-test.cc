#include <iostream>
#include <typeinfo>
#ifdef USE_STACK_TRACE_LOGGER
#include <glog/logging.h>
#endif

#include "pcg_random.hpp"

namespace kyawakyawa {
using RNG = pcg32;

static void ShowRandom(void) {
  RNG rng(42u, 54u);
  for (size_t i = 0; i < 50; i++) {
    uint32_t a = rng();
    std::cout << float(a) / UINT32_MAX << std::endl;
  }
}
}  // namespace kyawakyawa

int main(int argc, char** argv) {
#ifdef USE_STACK_TRACE_LOGGER
  (void)argc;
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
#else
  (void)argc;
  (void)argv;
#endif
  kyawakyawa::ShowRandom();
  return EXIT_SUCCESS;
}

