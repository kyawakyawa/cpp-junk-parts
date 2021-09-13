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

#include <stdlib.h>

#include <iostream>
#include <map>
#include <string>

#include "macros.h"

IGNORE_STRICT_WARNING_PUSH

#include <spdlog/spdlog.h>

#include "staticjson/staticjson.hpp"
#ifdef USE_STACK_TRACE_LOGGER
#include <glog/logging.h>
#endif

IGNORE_STRICT_WARNING_POP

struct A {
  int hoge = 0;
  std::string str;
};

namespace staticjson {

template <>
void init(A *t, staticjson::ObjectHandler *h) {
  h->add_property("hoge", &(t->hoge), staticjson::Flags::Optional);
  h->add_property("str", &(t->str), staticjson::Flags::Optional);
}

}  // namespace staticjson

/* json example

{
  "aaa": {
    "hoge" : 1,
    "str"  : "yes"
  },
  "bbb": {
    "hoge" : 0,
    "str" : "no"
  }
}

*/

int main([[maybe_unused]] int argc, [[maybe_unused]] char **argv) {
#ifdef USE_STACK_TRACE_LOGGER
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
#endif
  std::map<std::string, A> m;

  if (argc < 2) {
    std::string str =
        "{ \"aaa\": { \"hoge\" : 1,\"str\"  : \"yes\" }, \"bbb\": { \"hoge\" : "
        "0, \"str\" : \"no\" } }";

    staticjson::ParseStatus res;
    if (!staticjson::from_json_string(str.c_str(), &m, &res)) {
      spdlog::error("failed to load json file");
      return EXIT_FAILURE;
    }

    for (const auto &v : m) {
      spdlog::info("{}:\n\t{}\t{}", v.first, v.second.hoge, v.second.str);
    }

    return EXIT_SUCCESS;
  }

  printf("%s\n", argv[1]);

  staticjson::ParseStatus res;
  if (!staticjson::from_json_file(argv[1], &m, &res)) {
    spdlog::error("failed to load json file");
    return EXIT_FAILURE;
  }

  for (const auto &v : m) {
    spdlog::info("{}:\n\t{}\t{}", v.first, v.second.hoge, v.second.str);
  }

  return EXIT_SUCCESS;
}
