#!/bin/bash

module load gcc/8.3.1
module load cuda/11.1.0
module load ninja

build_dir=$(realpath $1)
prefix=$(realpath $2)
root_dir=$(realpath $3)
hpac_build_dir=$build_dir/hpac_build

echo $build_dir
echo $prefix
echo $root_dir
echo $hpac_build_dir


if [ ! -d $build_dir ]; then 
  pushd $build_dir
  mkdir $build_dir
  CC=gcc CPP=g++ cmake -G Ninja   \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo\
  -DLLVM_TARGETS_TO_BUILD='Native;NVPTX' \
  -DLLVM_ENABLE_PROJECTS='clang' \
  -DLLVM_ENABLE_RUNTIMES='openmp' \
  -DCMAKE_INSTALL_PREFIX=${prefix}\
  -DLLVM_OPTIMIZED_TABLEGEN=on   \
  -DBUILD_SHARED_LIBS=on   \
  -DLLVM_ENABLE_ASSERTIONS=On   \
  -DOPENMP_ENABLE_LIBOMPTARGET_PROFILING=off   \
  -DLLVM_CCACHE_BUILD=off    \
  -DLIBOMPTARGET_ENABLE_DEBUG='debug' ${root_dir}/llvm
else
  pushd $build_dir
fi

ninja -j 40 
ninja -j 40 install

popd

export PATH=${prefix}/bin/:$PATH
export LD_LIBRARY_PATH=${prefix}/lib/:$LD_LIBRARY_PATH


if [ ! -d $hpac_build_dir ]; then 
  mkdir -p $hpac_build_dir
  pushd $hpac_build_dir
  CC=clang CPP=clang++ cmake -G Ninja \
    -DCMAKE_INSTALL_PREFIX=$prefix \
    -DLLVM_EXTERNAL_CLANG_SOURCE_DIR=${root_dir}/clang/ \
    -DWITH_TORCH=On -DTorch_DIR=$(spack location -i py-torch)/lib/python3.8/site-packages/torch/share/cmake/Torch \
    -DPACKAGE_VERSION=15.0.0git ${root_dir}/approx
else
  pushd $hpac_build_dir
fi

ninja -j 40 
ninja -j 40 install

popd

echo "module load gcc/8.3.1" > env.sh
echo "module load cuda/11.1.0" >> env.sh
echo "export PATH=\"${prefix}/bin/:\$PATH\"" >> env.sh
echo "export LD_LIBRARY_PATH=\"${prefix}/lib/:\$LD_LIBRARY_PATH\"" >> env.sh
module load ninja

