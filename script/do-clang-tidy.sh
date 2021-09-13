#!/bin/bash

clang-tidy pc/*.cc --extra-arg=-std=c++20 --fix --fix-errors --header-filter=$(pwd)/inc/*.h
clang-tidy pc/*.c --fix --fix-errors --header-filter=$(pwd)/inc/*.h