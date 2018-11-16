#!/usr/bin/env bash

DIRNAME=../build
if [ -d "$DIRNAME" ]; then
    echo "dir exists"
    cd "$DIRNAME"
    cmake ".."
    make
else
    mkdir "$DIRNAME"
    cd "$DIRNAME"
    cmake ".."
    make
fi