#!/bin/bash


if [ "$#" != 2 ]; then
  echo "Wrong number of arguments"
  echo "Command line should be:"
  echo "$0 'path to installation directory' 'number of threads'"
  exit
fi

prefix=$1
threads=$2
current_dir=$(pwd)

NOCOLOR='\033[0m'
RED='\033[0;31m'
GREEN='\033[0;32m'
ORANGE='\033[0;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
LIGHTGRAY='\033[0;37m'
DARKGRAY='\033[1;30m'
LIGHTRED='\033[1;31m'
LIGHTGREEN='\033[1;32m'
YELLOW='\033[1;33m'
LIGHTBLUE='\033[1;34m'
LIGHTPURPLE='\033[1;35m'
LIGHTCYAN='\033[1;36m'
WHITE='\033[1;37m'

clang_bin=$prefix/bin/clang
approx_runtime_lib=$prefix/lib/libapprox.so

if [ ! -f $clang_bin ]; then
  mkdir -p build_compiler
  mkdir -p $prefix
  pushd build_compiler
  cmake -G Ninja \
    -DCMAKE_INSTALL_PREFIX=$prefix \
    -DLLVM_CCACHE_BUILD='Off'\
    -DCMAKE_EXPORT_COMPILE_COMMANDS='On'\
    -DCMAKE_BUILD_TYPE='RelWithDebInfo' \
    -DLLVM_FORCE_ENABLE_STATS='On' \
    -DLLVM_ENABLE_PROJECTS='clang;openmp' \
    -DLLVM_OPTIMIZED_TABLEGEN='On' \
    -DCLANG_BUILD_EXAMPLES='On' \
    -DBUILD_SHARED_LIBS='On' \
    -DLLVM_ENABLE_ASSERTIONS='On' \
    ../llvm

    ninja -j $threads
    ninja -j $threads install 
    popd
    echo "#!/bin/bash" > puppet_env.sh
    echo "export PATH=$prefix/bin/:\$PATH" >> puppet_env.sh
    echo "export LD_LIBRARY_PATH=$prefix/lib/:\$LD_LIBRARY_PATH" >> puppet_env.sh
    echo "export CC=clang" >> puppet_env.sh
    echo "export CPP=clang++" >> puppet_env.sh
fi

source puppet_env.sh

if [ ! -f $approx_runtime_lib ]; then
  mkdir build_hpac
  pushd build_hpac
  CC=clang CPP=clang++ cmake -G Ninja \
      -DCMAKE_INSTALL_PREFIX=$prefix \
      -DLLVM_EXTERNAL_CLANG_SOURCE_DIR=${current_dir}/clang/ \
      -DPACKAGE_VERSION=11.0.0git \
     ../approx
    ninja -j $threads
    ninja -j $threads install
    popd
fi
