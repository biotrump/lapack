#!/bin/bash
# Compiles ffts for Android
# Make sure you have NDK_ROOT defined in .bashrc or .bash_profile
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
export NDK_ROOT=/home/thomas/aosp/NDK/android-ndk-r9

INSTALL_DIR="`pwd`/java/android/bin"

if [ -f make.inc.armv7-a ]; then
cp -f make.inc.armv7-a make.inc
fi

#ANDROID_APIVER=android-8
#ANDROID_APIVER=android-9
#android 4.0.1 ICS and above
ANDROID_APIVER=android-14
#TOOL="4.6"
#gfortran is in r9d V4.8.0
TOOL="4.8.0"

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

arm=${arm:-arm}
echo arm=$arm
case arm in
  arm)
    TARGPLAT=arm-linux-androideabi
    ARCHI=arm
    CONFTARG=arm-eabi
  ;;
  x86)
    TARGPLAT=x86
    ARCHI=x86
    CONFTARG=x86
  ;;
  mips)
  ## probably wrong
    TARGPLAT=mipsel-linux-android
    ARCHI=mips
    CONFTARG=mips
  ;;
  *) echo $0: Unknown target; exit
esac

: ${NDK_ROOT:?}

echo "Using: $NDK_ROOT/toolchains/${TARGPLAT}-${TOOL}/prebuilt/${HOSTPLAT}/bin"
export ARCHI
#export PATH="$NDK_ROOT/toolchains/${TARGPLAT}-${TOOL}/prebuilt/${HOSTPLAT}/bin/:\
#$NDK_ROOT/toolchains/${TARGPLAT}-${TOOL}/prebuilt/${HOSTPLAT}/${TARGPLAT}/bin/:$PATH"
export PATH="$NDK_ROOT/toolchains/${TARGPLAT}-${TOOL}/prebuilt/${HOSTPLAT}/bin/:$PATH"
export SYS_ROOT="$NDK_ROOT/platforms/${ANDROID_APIVER}/arch-${ARCHI}/"
export CC="${TARGPLAT}-gcc --sysroot=$SYS_ROOT"
export LD="${TARGPLAT}-ld"
export AR="${TARGPLAT}-ar"
export ARCH=${AR}
export RANLIB="${TARGPLAT}-ranlib"
export STRIP="${TARGPLAT}-strip"
#export CFLAGS="-Os -fPIE"
export CFLAGS="-Os -fPIE --sysroot=$SYS_ROOT"
export CXXFLAGS="-fPIE --sysroot=$SYS_ROOT"
export FORTRAN="${TARGPLAT}-gfortran --sysroot=$SYS_ROOT"
#include path :
#platforms/android-21/arch-arm/usr/include/

BLISLIB=${BLISLIB:-/home/thomas/build/libFLAME/blis/lib/armv7a/libblis.a}
export BLISLIB
echo $BLISLIB

make clean

make -j${CORE_COUNT}

#pushd BLAS/SRC
#make clean
#make -j${CORE_COUNT}
#popd

#make -j${CORE_COUNT} blaslib

#make -j${CORE_COUNT} cblas_example

make -j${CORE_COUNT} blas_testing_android
make -j${CORE_COUNT} eig_testing_android
make -j${CORE_COUNT} lin_testing_android

#pushd CBLAS
#make clean
#make -j${CORE_COUNT}
#popd

#make -j${CORE_COUNT} cblaslib
#make -j${CORE_COUNT} lapacklib

#make -j${CORE_COUNT} lapacke_example

