#!/bin/bash
set -e

ver=$1
binary=$2
if [ "$ver" == "" ] || [ "$binary" == "" ]; then
	echo "Usage: $0 version binary" >&2
	exit 1
fi

rm -rf .build-release-mac
mkdir .build-release-mac
cd .build-release-mac

tar xf ../archives/"NEXTNet-v$ver-source.tar.gz"
mkdir build
cd build

cmake \
	-DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" \
	-DCMAKE_BUILD_TYPE=Release \
	-DBOOST_OVERRIDE=~/Installs/boost/1.86-headeronly/include \
	../

make -j 16

cd ..
cd ..
cp .build-release-mac/build/nextnet "$binary"

rm -rf .build-release-mac
