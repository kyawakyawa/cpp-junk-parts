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

#include <fstream>
#include <iostream>

#include "macros.h"

#ifdef USE_STACK_TRACE_LOGGER
#include <glog/logging.h>
#endif

IGNORE_STRICT_WARNING_PUSH

#include <spdlog/spdlog.h>

#include <nlohmann/json.hpp>

IGNORE_STRICT_WARNING_POP

namespace kyawakyawa {
static std::vector<uint8_t> load_msgpack(const std::string& filename) {
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  if (!ifs.is_open()) {
    spdlog::error("cannnot laod the file at {}", filename);
    throw std::runtime_error("cannnot laod the file at " + filename);
  }
  std::vector<uint8_t> msgpack;
  while (true) {
    uint8_t buffer = 0;
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
