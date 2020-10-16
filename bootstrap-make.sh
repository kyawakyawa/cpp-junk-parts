#!/bin/bash

rm -rf cmake-build-debug
rm -rf cmake-build-release

CC=clang CXX=clang++ cmake -Bcmake-build-debug -H. \
                           -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
                           -DCMAKE_BUILD_TYPE=Debug \
                           -DBUILD_WITH_MARCH_NATIVE=ON \
                           -DUSE_CCACHE=ON \
                           -DUSE_CPP20=ON \
                           -DUSE_SANITIZER=ON \
                           -DUSE_STACK_TRACE_LOGGER=OFF

mv cmake-build-debug/compile_commands.json .

CC=clang CXX=clang++ cmake -Bcmake-build-release -H. \
                           -DCMAKE_EXPORT_COMPILE_COMMANDS=0 \
                           -DCMAKE_BUILD_TYPE=Release \
                           -DBUILD_WITH_MARCH_NATIVE=ON \
                           -DUSE_CCACHE=ON \
                           -DUSE_CPP20=ON \
                           -DUSE_SANITIZER=OFF \
                           -DUSE_STACK_TRACE_LOGGER=OFF
