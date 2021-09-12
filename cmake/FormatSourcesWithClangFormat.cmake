# reference https://qiita.com/tenmyo/items/f8548ee9bab78f18cd25

find_program(CLANG_FORMAT_EXE clang-format)

function(clang_format target sources)
  if(CLANG_FORMAT_EXE)
    message(STATUS "Enable Clang-Format ${target}")
    # message(STATUS "${CLANG_FORMAT_EXE} -i -style=file ${sources}")
    add_custom_target(
      "${target}_format-with-clang-format"
      COMMAND ${CLANG_FORMAT_EXE} -i -style=file ${sources}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      )
    add_dependencies(${target} "${target}_format-with-clang-format")
  endif(CLANG_FORMAT_EXE)
endfunction()
