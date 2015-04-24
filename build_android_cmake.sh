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

#default is arm
arm=${arm:-arm}
echo arm=$arm
case arm in
  arm)
    TARGPLAT=arm-linux-androideabi
    CONFTARG=arm-eabi
    ARCHI=arm
  ;;
  x86)
    TARGPLAT=i686-linux-android
    CONFTARG=x86
    ARCHI=x86
  ;;
  mips)
  ## probably wrong
    TARGPLAT=mipsel-linux-android
    CONFTARG=mips
    ARCHI=mips
  ;;
  *) echo $0: Unknown target; exit
esac
echo ARCHI=$ARCHI

#: ${NDK_ROOT:?}

echo "Using: $NDK_ROOT/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin"
export ARCHI
#export PATH="$NDK_ROOT/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin/:\
#$NDK_ROOT/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/${TARGPLAT}/bin/:$PATH"
export PATH="${NDK_ROOT}/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin/:$PATH"
echo $PATH
export SYS_ROOT="${NDK_ROOT}/platforms/${ANDROID_APIVER}/arch-${ARCHI}/"
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

#!!! quite importnat for cmake to define the NDK's fortran compiler.!!!
#Don't let cmake decide it.
export FC=${FORTRAN}


if [ -f make.inc.armv7-a ]; then
cp -f make.inc.armv7-a make.inc
fi

if [ ! -d build_and ]; then
mkdir build_and
else
rm -rf build_and/*
fi
pushd build_and

cmake -DCMAKE_TOOLCHAIN_FILE=${LAPACK_SRC}/android.toolchain.cmake -DANDROID_NATIVE_API_LEVEL=${ANDROID_APIVER} \
-DANDROID_NDK=${NDK_ROOT} -DANDROID_TOOLCHAIN_NAME=${TARGPLAT}-${TOOL_VER} \
-DCMAKE_BUILD_TYPE=Release -DANDROID_ABI="armeabi-v7a with VFPV3" ..

make -j${CORE_COUNT}
#make

