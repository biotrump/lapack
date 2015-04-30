#!/bin/bash
# Make sure you have NDK_ROOT defined in .bashrc or .bash_profile

#export CMAKE_BUILD_TYPE "Debug"
export CMAKE_BUILD_TYPE="Release"

#get cpu counts
case $(uname -s) in
  Darwin)
    CONFBUILD=i386-apple-darwin`uname -r`
    HOSTPLAT=darwin-x86
    CORE_COUNT=`sysctl -n hw.ncpu`
  ;;
  Linux)
    CONFBUILD=x86-unknown-linux
    HOSTPLAT=linux-`uname -m`
    CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
  ;;
CYGWIN*)
	CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
	;;
  *) echo $0: Unknown platform; exit
esac

export LAPACK_SRC=`pwd`

# Modify INSTALL_DIR to suit your situation
#Lollipop	5.0 - 5.1	API level 21, 22
#KitKat	4.4 - 4.4.4	API level 19
#Jelly Bean	4.3.x	API level 18
#Jelly Bean	4.2.x	API level 17
#Jelly Bean	4.1.x	API level 16
#Ice Cream Sandwich	4.0.3 - 4.0.4	API level 15, NDK 8
#Ice Cream Sandwich	4.0.1 - 4.0.2	API level 14, NDK 7
#Honeycomb	3.2.x	API level 13
#Honeycomb	3.1	API level 12, NDK 6
#Honeycomb	3.0	API level 11
#Gingerbread	2.3.3 - 2.3.7	API level 10
#Gingerbread	2.3 - 2.3.2	API level 9, NDK 5
#Froyo	2.2.x	API level 8, NDK 4

#NDK_ROOT=${NDK_ROOT:-/home/thomas/aosp/NDK/android-ndk-r10d}
#export NDK_ROOT
#gofortran is supported in r9
export NDK_ROOT=${HOME}/NDK/android-ndk-r9
export ANDROID_NDK=${NDK_ROOT}

if [[ ${NDK_ROOT} =~ .*"-r9".* ]]
then
#ANDROID_APIVER=android-8
#ANDROID_APIVER=android-9
#android 4.0.1 ICS and above
ANDROID_APIVER=android-14
#TOOL_VER="4.6"
#gfortran is in r9d V4.8.0
TOOL_VER="4.8.0"
else
#android 4.0.1 ICS and above
ANDROID_APIVER=android-14
TOOL_VER="4.9"
fi

if [ $# -ge 1 ]; then
	export ARCHI=$1
else
#default
	export ARCHI=arm
fi
echo ARCHI=$ARCHI

#default is arm
case $ARCHI in
  arm)
    TARGPLAT=arm-linux-androideabi
    CONFTARG=arm-eabi
	echo "Using: $NDK_ROOT/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin"
	#export PATH="$NDK_ROOT/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin/:\
	#$NDK_ROOT/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/${TARGPLAT}/bin/:$PATH"
	export PATH="${NDK_ROOT}/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin/:$PATH"
  ;;
  x86)
    TARGPLAT=i686-linux-android
    CONFTARG=x86
	echo "Using: $NDK_ROOT/toolchains/x86-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin"
	export PATH="${NDK_ROOT}/toolchains/x86-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin/:$PATH"
#specify assembler for x86 SSE3, but ffts's sse.s needs 64bit x86.
#intel atom z2xxx and the old atoms are 32bit, so 64bit x86 in android can't work in
#most atom devices.
#http://forum.cvapp.org/viewtopic.php?f=13&t=423&sid=4c47343b1de899f9e1b0d157d04d0af1
#	export  CCAS="${TARGPLAT}-as"
#	export  CCASFLAGS="--64 -march=i686+sse3"
#	export  CCASFLAGS="--64"

  ;;
  mips)
  ## probably wrong
    TARGPLAT=mipsel-linux-android
    CONFTARG=mips
  ;;
  *) echo $0: Unknown target; exit
esac
#: ${NDK_ROOT:?}
echo $PATH

export SYS_ROOT="${NDK_ROOT}/platforms/${ANDROID_APIVER}/arch-${ARCHI}/"
export CC="${TARGPLAT}-gcc --sysroot=$SYS_ROOT"
export LD="${TARGPLAT}-ld"
export AR="${TARGPLAT}-ar"
export ARCH=${AR}
export RANLIB="${TARGPLAT}-ranlib"
export STRIP="${TARGPLAT}-strip"
#export CFLAGS="-Os -fPIE"
export CFLAGS="-Os -fPIE -fPIC --sysroot=$SYS_ROOT"
export CXXFLAGS="-fPIE -fPIC --sysroot=$SYS_ROOT"
export FORTRAN="${TARGPLAT}-gfortran --sysroot=$SYS_ROOT"

#!!! quite importnat for cmake to define the NDK's fortran compiler.!!!
#Don't let cmake decide it.
export FC=${FORTRAN}

case $ARCHI in
  arm)
	if [ -f make.inc.armv7-a ]; then
	cp -f make.inc.armv7-a make.inc
	fi

	if [ ! -d build_NDK_$ARCHI ]; then
	mkdir build_NDK_$ARCHI
	else
	rm -rf build_NDK_$ARCHI/*
	fi
	pushd build_NDK_$ARCHI

	cmake -DCMAKE_TOOLCHAIN_FILE=${LAPACK_SRC}/android.toolchain.cmake -DANDROID_NATIVE_API_LEVEL=${ANDROID_APIVER} \
	-DANDROID_NDK=${NDK_ROOT} -DANDROID_TOOLCHAIN_NAME=${TARGPLAT}-${TOOL_VER} \
	-DCMAKE_BUILD_TYPE=Release -DANDROID_ABI="armeabi-v7a with VFPV3" ..
  ;;
  x86)
	if [ -f make.inc.NDK-x86 ]; then
	cp -f make.inc.NDK-x86 make.inc
	fi

	if [ ! -d build_NDK_$ARCHI ]; then
	mkdir build_NDK_$ARCHI
	else
	rm -rf build_NDK_$ARCHI/*
	fi
	pushd build_NDK_$ARCHI

	cmake -DCMAKE_TOOLCHAIN_FILE=${LAPACK_SRC}/android.toolchain.cmake -DANDROID_NATIVE_API_LEVEL=${ANDROID_APIVER} \
	-DANDROID_NDK=${NDK_ROOT} -DANDROID_TOOLCHAIN_NAME=x86-${TOOL_VER} \
	-DCMAKE_BUILD_TYPE=Release -DANDROID_ABI="x86" ..
  ;;
  mips)
  ;;
  *) echo $0: Unknown target; exit
esac

make -j${CORE_COUNT}
#make

