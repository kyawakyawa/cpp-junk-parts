cmake_minimum_required(VERSION 3.10)
project(cpp-junk-parts LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release")
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# CCACHE
option(USE_CCACHE "Use ccache to accelerate a build(if available)" ON)
find_program(CCACHE_EXE ccache)
if(USE_CCACHE AND CCACHE_EXE)
    message(STATUS "ccache: Enabled")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_EXE}")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_EXE}")
else()
    message(STATUS "ccache: Disabled")
    set(USE_CCACHE OFF)
endif()

# Sanitizer
option(USE_SANITIZER "Enable Address/Memory sanitizer to check memory leak etc..." OFF)
if(USE_SANITIZER)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
    message(STATUS "Address/Memory sanitizer: ENABLED")
else()
    message(STATUS "Address/Memory sanitizer: DISABLED")
endif()

# Overall compile options
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Og")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

option(BUILD_WITH_MARCH_NATIVE "Enable architecture-aware optimization (-march=native)" OFF)
if(BUILD_WITH_MARCH_NATIVE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=native -march=native")
    message(STATUS "Architecture-aware optimization (-march=native): ENABLED")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=native")
    message(STATUS "Architecture-aware optimization (-march=native): DISABLED")
endif()

# Check C++17 support
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++17" COMPILER_SUPPORTS_CXX17)
if (COMPILER_SUPPORTS_CXX17)
    message(STATUS "C++17 support: OK (-std=c++17)")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++17 support. Please use a different C++ compiler.")
endif()

# find package and check third party
# - json
# - spdlog
# - pcg-cpp
# - google glog (optional)

# json
if (NOT EXISTS "${PROJECT_SOURCE_DIR}/third_party/json/CMakeLists.txt")
    message(FATAL_ERROR "json was not found!"
        "Please execute the command \"git submodule update --init --recursive\"")
endif()
set(JSON_BuildTests OFF CACHE INTERNAL "")
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/json)

# spdlog
if (NOT EXISTS "${PROJECT_SOURCE_DIR}/third_party/spdlog/include")
    message(FATAL_ERROR "spdlog was not found!"
        "Please execute the command \"git submodule update --init --recursive\"")
endif()
set(SPDLOG_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third_party/spdlog/include)

set(USE_STACK_TRACE_LOGGER OFF CACHE BOOL "Enable automatic stack trace logger of google-glog")
if(USE_STACK_TRACE_LOGGER)
    message(STATUS "Stack trace logger: ENABLED")
    # Glog
    find_package(Glog REQUIRED)
    # include directories is ${GLOG_INCLUDE_DIR}
else()
    message(STATUS "Stack trace logger: DISABLED")
endif()

# pcg-cpp
if (NOT EXISTS "${PROJECT_SOURCE_DIR}/third_party/pcg-cpp/include")
  message(FATAL_ERROR "pcg-cpp was not found!"
      "Please execute the command \"git submodule update --init --recursive\"")
endif()
set(PCG_CPP_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third_party/pcg-cpp/include)

# build selection
add_subdirectory(pc)

#message(STATUS "*** dump start cmake variables ***")
#get_cmake_property(_variableNames VARIABLES)
#foreach(_variableName ${_variableNames})
#        message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#message(STATUS "*** dump end ***")
#

