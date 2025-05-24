#!/bin/bash
set -e

ver=$1
binary=$2
if [ "$ver" == "" ] || [ "$binary" == "" ]; then
	echo "Usage: $0 version binary" >&2
	exit 1
fi

echo "*** Building container holy-build-box-boost using compile-release-docker/Dockerfile"
docker build -t holy-build-box-boost compile-release-docker

echo "*** Building NEXTNet-v$ver"
docker run -t -i --rm \
  -v `pwd`:/io \
  holy-build-box-boost \
  bash /io/scripts/build-release-docker-hbb.sh $ver "$binary"
