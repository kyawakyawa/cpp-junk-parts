cmake_minimum_required(VERSION 3.10)
project(cpp-junk-parts LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

option(BUILD_WITH_MARCH_NATIVE "Enable architecture-aware optimization (-march=native)" OFF)

option(USE_CCACHE "Use ccache to accelerate a build(if available)" ON)
option(USE_CPP20 "use c++20" OFF)
option(USE_SANITIZER "Enable Address/Memory sanitizer to check memory leak etc..." OFF)
option(USE_STACK_TRACE_LOGGER "Enable automatic stack trace logger of google-glog" OFF)

# build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release")
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# CCACHE
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
if (MSVC)
else (MSVC)
  if(USE_SANITIZER)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
      message(STATUS "Address/Memory sanitizer: ENABLED")
  else(USE_SANITIZER)
      message(STATUS "Address/Memory sanitizer: DISABLED")
  endif(USE_SANITIZER)
endif(MSVC)

if (MSVC)
else (MSVC)
  # Check C++17 support
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-std=c++17" COMPILER_SUPPORTS_CXX17)
  if (COMPILER_SUPPORTS_CXX17)
      message(STATUS "C++17 support: OK (-std=c++17)")
  else()
      message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++17 support. Please use a different C++ compiler.")
  endif()
  
  # Check C++20 support if USE_CPP20 is ON
  if (USE_CPP20)
      check_cxx_compiler_flag("-std=c++20" COMPILER_SUPPORTS_CXX20)
      if (COMPILER_SUPPORTS_CXX20)
          message(STATUS "C++20 support: OK (-std=c++20)")
      else()
          message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++20 support. Please use a different C++ compiler.")
      endif()
  endif()
endif(MSVC)

# find package and check third party
# - thread
# - eigen
# - json
# - StaticJSON + rapidjson
# - spdlog
# - pcg-cpp
# - oneTBB
# - google glog (optional)

set(GIT_SUBMODULE_UPDATE_COMMAND "git submodule update --init --recursive --depth 1")

# Thread
find_package(Threads REQUIRED)

# Eigen
find_package(Eigen3 3.3 REQUIRED)

# GLFW3
find_package( glfw3 REQUIRED )

# GLEW
find_package( GLEW REQUIRED )

# json
if (NOT EXISTS "${PROJECT_SOURCE_DIR}/third_party/json/CMakeLists.txt")
    message(FATAL_ERROR "[error] json was not found!\n"
        "Please execute the command\n\"${GIT_SUBMODULE_UPDATE_COMMAND}\"")
endif()
set(JSON_BuildTests OFF CACHE INTERNAL "")
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/json)

# Static JSON + rapidjson

set (POSSIBLE_RAPIDJSON_DIR "${PROJECT_SOURCE_DIR}/third_party/rapidjson/include")
find_path(RAPIDJSON_INCLUDE_DIR rapidjson/rapidjson.h PATHS ${POSSIBLE_RAPIDJSON_DIR})
if (NOT RAPIDJSON_INCLUDE_DIR)
    message("Rapid JSON was not found")
endif()

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/third_party/StaticJSON/src/staticjson.cpp")
    message(FATAL_ERROR "[error] StaticJson was not found!\n"
        "Please execute the command\n\"${GIT_SUBMODULE_UPDATE_COMMAND}\"")
endif()

set(STATICJSON_SOURCE_FILES "${PROJECT_SOURCE_DIR}/third_party/StaticJSON/src/staticjson.cpp")
add_library(staticjson SHARED ${STATICJSON_SOURCE_FILES})

target_include_directories(staticjson PUBLIC ${RAPIDJSON_INCLUDE_DIR})
target_include_directories(staticjson PUBLIC ${PROJECT_SOURCE_DIR}/third_party/StaticJSON/include)

target_compile_options(staticjson PRIVATE
          $<$<CXX_COMPILER_ID:Clang>:  -Wall -Wextra -pedantic -g>
          $<$<CXX_COMPILER_ID:GCC>:    -Wall -Wextra -pedantic -g>
          $<$<CXX_COMPILER_ID:MSVC>:>)
 
target_compile_features(staticjson
                        PRIVATE cxx_std_11)   

# spdlog
set(SPDLOG_BUILD_SHARED ON CACHE BOOL "Build libspdlog shared library")
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/spdlog)

# google glog
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
  message(FATAL_ERROR "[error] pcg-cpp was not found!\n"
      "Please execute the command\n\"${GIT_SUBMODULE_UPDATE_COMMAND}\"")
endif()
set(PCG_CPP_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third_party/pcg-cpp/include)

# include(${PROJECT_SOURCE_DIR}/cmake/oneTBB/TBBGet.cmake)
# include(${PROJECT_SOURCE_DIR}/cmake/oneTBB/TBBBuild.cmake)
# tbb_get(TBB_ROOT tbb_root SOURCE_CODE)
# tbb_build(TBB_ROOT ${tbb_root} CONFIG_DIR TBB_DIR)
# find_package(TBB REQUIRED)

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/third_party/oneTBB/cmake/TBBBuild.cmake")
  message(FATAL_ERROR "\n[error] oneTBB was not found!\n"
      "Please execute the command\n\"${GIT_SUBMODULE_UPDATE_COMMAND}\"")
endif()
include(${PROJECT_SOURCE_DIR}/third_party/oneTBB/cmake/TBBBuild.cmake)
tbb_build(TBB_ROOT ${PROJECT_SOURCE_DIR}/third_party/oneTBB CONFIG_DIR TBB_DIR)
find_package(TBB REQUIRED)

# build selection
add_subdirectory(pc)

#message(STATUS "*** dump start cmake variables ***")
#get_cmake_property(_variableNames VARIABLES)
#foreach(_variableName ${_variableNames})
#        message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#message(STATUS "*** dump end ***")
#

