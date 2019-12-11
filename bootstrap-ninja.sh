#!/bin/bash

rm -rf cmake-build-debug

CC=clang CXX=clang++ cmake -Bcmake-build-debug -H. \
                           -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
                           -DCMAKE_BUILD_TYPE=Debug \
                           -DUSE_CCACHE=ON \
                           -DUSE_SANITIZER=OFF \
                           -DUSE_STACK_TRACE_LOGGER=ON \
                           -DBUILD_WITH_MARCH_NATIVE=ON \
                           -GNinja \
                           ..

mv cmake-build-debug/compile_commands.json .
