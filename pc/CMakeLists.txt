set(EXECUTABLE_TARGETS "")

# executables
# - glfw-shader
# - glfw-texture
# - load_messagepack
# - parallelism-ts-test
# - pcg32-test
# - popcount
# - quick-select
# - spdlog-test
# - staticjson-example
# - std-variant-test 

add_executable(glfw-shader ./glfw-shader.cc)
target_link_libraries(glfw-shader PRIVATE GL GLEW glfw)
list(APPEND EXECUTABLE_TARGETS glfw-shader)

add_executable(glfw-texture ./glfw-texture.cc)
target_link_libraries(glfw-texture PRIVATE GL GLEW glfw)
list(APPEND EXECUTABLE_TARGETS glfw-texture)

add_executable(load_messagepack load_messagepack.cc)
target_link_libraries(load_messagepack PRIVATE nlohmann_json::nlohmann_json)
list(APPEND EXECUTABLE_TARGETS load_messagepack)

add_executable(parallelism-ts-test ./parallelism-ts-test.cc)
target_link_libraries(parallelism-ts-test PRIVATE ${TBB_IMPORTED_TARGETS})
list(APPEND EXECUTABLE_TARGETS parallelism-ts-test)

add_executable(pcg32-test ./pcg32-test.cc)
target_include_directories(pcg32-test SYSTEM PRIVATE ${PCG_CPP_INCLUDE_DIR})
target_compile_options(pcg32-test PRIVATE -Wno-date-time)
list(APPEND EXECUTABLE_TARGETS pcg32-test)

if (USE_CPP20)
    add_executable(popcount ./popcount.cc)
    list(APPEND EXECUTABLE_TARGETS popcount)
endif()

add_executable(quick-select ./quick-select.cc)
target_compile_definitions(quick-select PRIVATE  _LIBCPP_DEBUG _GLIBCXX_DEBUG)
list(APPEND EXECUTABLE_TARGETS quick-select)

if (USE_CPP20)
  add_executable(spdlog-test ./spdlog-test.cc)
  list(APPEND EXECUTABLE_TARGETS spdlog-test)
endif()

add_executable(staticjson-example ./staticjson-example.cc)
target_link_libraries(staticjson-example PRIVATE staticjson)
list(APPEND EXECUTABLE_TARGETS staticjson-example)

add_executable(std-variant-test ./std-variant-test.cc)
list(APPEND EXECUTABLE_TARGETS std-variant-test)

foreach(EXECUTABLE_TARGET IN LISTS EXECUTABLE_TARGETS)
    # set output directory for executables
    set_target_properties(${EXECUTABLE_TARGET} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

    if (USE_CPP20)
        target_compile_features(${EXECUTABLE_TARGET} PRIVATE cxx_std_20)
    else()
        target_compile_features(${EXECUTABLE_TARGET} PRIVATE cxx_std_17)
    endif()

    target_compile_options(${EXECUTABLE_TARGET}
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>:-Wall -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded>
            $<$<CXX_COMPILER_ID:GCC>:-pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Werror -Wno-unused>
            $<$<CXX_COMPILER_ID:MSVC>:/W4>)
    target_link_libraries(${EXECUTABLE_TARGET} PRIVATE spdlog::spdlog)
    if(USE_STACK_TRACE_LOGGER)
        target_compile_definitions(${EXECUTABLE_TARGET}
                PRIVATE -DUSE_STACK_TRACE_LOGGER)
        target_link_libraries(${EXECUTABLE_TARGET} PRIVATE glog::glog)
    endif()
endforeach()
