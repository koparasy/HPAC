# HPAC

HPAC is comprised of two components, namely the core and the harness. The core implements the approximate 
programming model. It extends the Clang/LLVM compiler and provides runtime support. The harness 
facilitates easy exploration of the approximate design space.

Currently, HPAC supports:

- Perforation:
   - ini: Drop the first iterations of a for loop.
   - fini: Drop the last iterations of a for loop.
   - random: Drop randomly iterations of a for loop.
   - small: Drop 1 iteration every N iterations.
   - large: Drop N-1 iterations every N iterations.
- Memoization:
   - iACT: Apply approximate memoization based on the inputs of a code region.
   - TAF: Apply approximate memoization based on the outputs of a code region.


## Build Requirements:

HPAC on a x86_64 system has been successfully built with the following dependencies: 
following free projects:
- Ninja 1.9.0
- CMAKE 3.14
- gcc 8.3.1
- python 3.7.6
- pandas 1.2.2
- numpy 1.20.1
- matplotlib 3.3.4
- seaborn-0.11.1
- yaml 5.4.1

## Build HPAC:

### Build HPAC Compiler Extensions

To build the HPAC compiler and all the infrastructure please execute the following commands:

```bash
git clone git@github.com:koparasy/HPAC.git
cd HPAC
./setup.sh 'PREFIX' 'NUM THREADS' 
```

The 'PREFIX' argument defines where to install all the HPAC related binaries and executables. The 'NUM THREADS' parameter
define how many threads should the installation use. The installation script performs the following actions:

1. Configures, builds and installs clang/LLVM/OpenMP including the approximation extensions.
2. Configures, builds and installs the approximation library. 
3. Creates a file, called 'hpac_env.sh', at the root of the project which should always be sourced before using HPAC.
4. Fetches the original applications of our evaluation.
5. Patches the applications with our extensions.
6. Configures and compiler the extended approximated applications.

The installation process can take a considerable amount of time. In our tested system, the complete 
installation took ~5 hours. 


## Replicate Analysis of SC-22

The software and the data are being reviewed for open source release at the moment. Once, the release is done we will update the instructions.

## Contributing
To contribute to this repo please send a [pull
request](https://help.github.com/articles/using-pull-requests/) on the
develop branch of this repo.

## Authors

This code was created by Konstantinos Parasyris (parasyris1@llnl.gov) and Giorgis Georgakoudis,
(georgakoudis1@llnl.gov), assisted with design input from Harshitha Menon (gopalakrishn1@llnl.gov).


## License

This repo is distributed under the terms of the Apache License (Version
2.0) with LLVM exceptions. Other software that is part of this
repository may be under a different license, documented by the file
LICENSE in its sub-directory.

All new contributions to this repo must be under the Apache License (Version 2.0) with LLVM exceptions.

See files [LICENSE](LICENSE) and [NOTICE](NOTICE) for more information.

SPDX License Identifier: "Apache-2.0 WITH LLVM-exception"

LLNL-CODE- 825539
