#!/bin/bash
set -e
set -o pipefail 

ver=$1
if [ "$ver" == "" ]; then
	echo "Usage: $0 version" >&2
	exit 1
fi

if [ $(git ls-remote --tags origin v$ver | wc -l) != 0 ]; then
	echo "Re-using exiting tag for version $ver" >&2
else
	./scripts/tag-release.sh $ver
fi

echo "Building full archive for version $ver (including submodules)" >&2
mkdir -p archives
rm -f archives/NEXTNet-v$ver-full.tar
./scripts/git-archive-all.sh --tree-ish v$ver archives/NEXTNet-v$ver-full.tar
gzip --best archives/NEXTNet-v$ver-full.tar
