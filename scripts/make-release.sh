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
./scripts/git-archive-all.sh --prefix NEXTNet/ --tree-ish v$ver "archives/NEXTNet-v$ver-source.tar"
gzip --best "archives/NEXTNet-v$ver-source.tar"

echo "Building Linux x86_64 release binaries for version $ver" >&2
mkdir -p binaries/linux-x86_64
./scripts/build-release-linux.sh $ver binaries/linux-x86_64/nextnet-v$ver

echo "Building Mac x86_64+arm64 release binaries for version $ver" >&2
mkdir -p binaries/mac
./scripts/build-release-mac.sh $ver binaries/mac/nextnet-v$ver

for arch in linux-x86_64 mac; do
	echo "Creating NEXTNet-v$ver-$arch.tar.gz"
	rm -rf .build-binary-archive-$arch
	mkdir .build-binary-archive-$arch
	cp binaries/$arch/nextnet-v$ver .build-binary-archive-$arch/nextnet
	(cd .build-binary-archive-$arch; tar czf ../archives/NEXTNet-v$ver-$arch.tar.gz nextnet)
	rm -rf .build-binary-archive-$arch
done
