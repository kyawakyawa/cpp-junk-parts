set(EXECUTABLE_TARGETS "")
add_executable(load_messagepack load_messagepack.cc)
list(APPEND EXECUTABLE_TARGETS load_messagepack)

foreach(EXECUTABLE_TARGET IN LISTS EXECUTABLE_TARGETS)
endforeach()
