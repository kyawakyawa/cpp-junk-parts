#!/bin/bash

# https://qiita.com/rtakasuke/items/fd35ef60370e3d8c225d
set -e

clang-tidy pc/*.cc --extra-arg=-std=c++20 --fix --fix-errors --header-filter=$(pwd)/inc/*.h
clang-tidy pc/*.c --fix --fix-errors --header-filter=$(pwd)/inc/*.h