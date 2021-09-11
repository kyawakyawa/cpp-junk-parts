set(EXECUTABLE_TARGETS "")

# executables
# - cmp-function-functor-lambda
# - cubic-equation
# - glfw-shader
# - glfw-texture
# - load_messagepack
# - parallelism-ts-test
# - pcg32-test
# - popcount
# - qr-decomposition
# - quaric-equation
# - quick-select
# - spdlog-test
# - staticjson-example
# - std-variant-and-concept 
# - std-variant-test 

add_executable(cmp-function-functor-lambda ./cmp-function-functor-lambda.cc)
list(APPEND EXECUTABLE_TARGETS cmp-function-functor-lambda)

add_executable(cubic-equation ./cubic-equation.cc)
list(APPEND EXECUTABLE_TARGETS cubic-equation)

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

add_executable(parallel-examples ./parallel-examples.cc)
target_link_libraries(parallel-examples PRIVATE Threads::Threads ${TBB_IMPORTED_TARGETS})
list(APPEND EXECUTABLE_TARGETS parallel-examples)

add_executable(pcg32-test ./pcg32-test.cc)
target_include_directories(pcg32-test PRIVATE ${PCG_CPP_INCLUDE_DIR})
target_compile_options(pcg32-test PRIVATE -Wno-date-time)
list(APPEND EXECUTABLE_TARGETS pcg32-test)

if (USE_CPP20)
    add_executable(popcount ./popcount.cc)
    list(APPEND EXECUTABLE_TARGETS popcount)
endif(USE_CPP20)

add_executable(quick-select ./quick-select.cc)
target_compile_definitions(quick-select PRIVATE  _LIBCPP_DEBUG _GLIBCXX_DEBUG)
list(APPEND EXECUTABLE_TARGETS quick-select)

add_executable(qr-decomposition ./qr-decomposition.cc)
target_link_libraries(qr-decomposition PRIVATE Eigen3::Eigen)
list(APPEND EXECUTABLE_TARGETS qr-decomposition)

add_executable(quartic-equation quartic-equation.cc)
list(APPEND EXECUTABLE_TARGETS quartic-equation)

if (USE_CPP20)
  add_executable(spdlog-test ./spdlog-test.cc)
  list(APPEND EXECUTABLE_TARGETS spdlog-test)
endif()

add_executable(staticjson-example ./staticjson-example.cc)
target_link_libraries(staticjson-example PRIVATE staticjson)
list(APPEND EXECUTABLE_TARGETS staticjson-example)

if (USE_CPP20)
    add_executable(std-variant-and-concept ./std-variant-and-concept.cc)
    list(APPEND EXECUTABLE_TARGETS std-variant-and-concept)
endif(USE_CPP20)

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
            $<$<CXX_COMPILER_ID:AppleClang>:-Wall -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded>
            $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wredundant-decls -Wcast-align -Wmissing-declarations -Wmissing-include-dirs -Wswitch-enum -Wswitch-default -Winvalid-pch -Wredundant-decls -Wformat=2 -Wmissing-format-attribute -Wformat-nonliteral>
            $<$<CXX_COMPILER_ID:MSVC>:/W4>)
    target_link_libraries(${EXECUTABLE_TARGET} PRIVATE spdlog::spdlog)

    target_include_directories(${EXECUTABLE_TARGET} PRIVATE ${PROJECT_SOURCE_DIR}/inc)

    if(USE_STACK_TRACE_LOGGER)
        target_compile_definitions(${EXECUTABLE_TARGET}
                PRIVATE -DUSE_STACK_TRACE_LOGGER)
        target_link_libraries(${EXECUTABLE_TARGET} PRIVATE glog::glog)
    endif()

endforeach()

set(C_PROGRAMS "")
add_executable(c-quartic-equation c-quartic-equation.c)
target_compile_options(c-quartic-equation
        PRIVATE
        $<$<CXX_COMPILER_ID:Clang>:-Ofast>
        $<$<CXX_COMPILER_ID:AppleClang>:-Ofast>
        $<$<CXX_COMPILER_ID:GNU>:-Ofast>
        $<$<CXX_COMPILER_ID:MSVC>:/fp:fast>)
 
list(APPEND C_PROGRAMS c-quartic-equation)

foreach(C_PROGRAM IN LISTS C_PROGRAMS)

    # set output directory for c programs
    set_target_properties(${C_PROGRAM} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

    target_compile_options(${C_PROGRAM}
            PRIVATE
            $<$<CXX_COMPILER_ID:Clang>:-Wall -Weverything -Wno-padded>
            $<$<CXX_COMPILER_ID:AppleClang>:-Wall -Weverything -Wno-padded>
            $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wredundant-decls -Wcast-align -Wmissing-declarations -Wmissing-include-dirs -Wswitch-enum -Wswitch-default -Winvalid-pch -Wredundant-decls -Wformat=2 -Wmissing-format-attribute -Wformat-nonliteral>
            $<$<CXX_COMPILER_ID:MSVC>:/W4>)

    target_include_directories(${C_PROGRAM} PRIVATE ${PROJECT_SOURCE_DIR}/inc)

    target_link_libraries(${C_PROGRAM} PRIVATE m)
    
endforeach()