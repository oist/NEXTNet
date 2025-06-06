#!/bin/bash
set -e

ver=$1
binary=$2

# Activate Holy Build Box environment.
source /hbb_exe/activate

set -x

# Extract and enter source
tar xzf /io/archives/"NEXTNet-v$ver-source.tar.gz"
cd NEXTNet

# Compile
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBOOST_OVERRIDE=/ext ../
make -j16

# Copy result to host
cp nextnet "/io/$binary"
strip "/io/$binary"
