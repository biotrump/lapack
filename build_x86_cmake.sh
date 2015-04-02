#!/bin/bash

case `uname` in
"Darwin")
	# Should also work on other BSDs
	CORE_COUNT=`sysctl -n hw.ncpu`
	;;
"Linux")
	CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
	;;
CYGWIN*)
	CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
	;;
*)
	echo Unsupported platform: `uname`
	exit -1
esac

#There are more make.inc.x in INSTALL
if [ -f make.inc.example ]; then
cp -f make.inc.example make.inc
fi

if [ ! -d build-x86 ]; then
mkdir build-x86
else
rm -rf build-x86/*
fi
pushd build-x86
cmake ..
make -j${CORE_COUNT}
popd
