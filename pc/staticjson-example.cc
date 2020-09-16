#include <stdlib.h>

#include <iostream>
#include <map>
#include <string>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#endif

#include <spdlog/spdlog.h>

#include "staticjson/staticjson.hpp"
#ifdef USE_STACK_TRACE_LOGGER
#include <glog/logging.h>
#endif

#ifdef __clang__
#pragma clang diagnostic pop
#endif

struct A {
  int hoge        = 0;
  std::string str = "";
};

namespace staticjson {

template <>
void init(A *d, staticjson::ObjectHandler *h) {
  h->add_property("hoge", &(d->hoge), staticjson::Flags::Optional);
  h->add_property("str", &(d->str), staticjson::Flags::Optional);
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
