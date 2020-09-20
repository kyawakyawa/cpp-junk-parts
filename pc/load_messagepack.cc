#include <fstream>
#include <iostream>

#ifdef USE_STACK_TRACE_LOGGER
#include <glog/logging.h>
#endif

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#include <spdlog/spdlog.h>

#include <nlohmann/json.hpp>
#pragma clang diagnostic pop

namespace kyawakyawa {
static std::vector<uint8_t> load_msgpack(const std::string& filename) {
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  if (!ifs.is_open()) {
    spdlog::error("cannnot laod the file at {}", filename);
    throw std::runtime_error("cannnot laod the file at " + filename);
  }
  std::vector<uint8_t> msgpack;
  while (true) {
    uint8_t buffer;
    ifs.read(reinterpret_cast<char*>(&buffer), sizeof(uint8_t));
    if (ifs.eof()) {
      break;
    }
    msgpack.emplace_back(buffer);
  }
  return msgpack;
}
}  // namespace kyawakyawa

int main(int argc, char** argv) {
#ifdef USE_STACK_TRACE_LOGGER
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
#endif

  if (argc <= 1) {
    spdlog::error("Please specify msg filename");
    return EXIT_FAILURE;
  }

  const std::string filename = argv[1];

  std::vector<uint8_t> msgpack = kyawakyawa::load_msgpack(filename);

  const auto json = nlohmann::json::from_msgpack(msgpack);

  std::cout << json;

  return 0;
}
