# Puppeteer 

Puppeteer is a set of extensions on the HPAC programming model to allow sensitivity analysis of applications in respect to approximation errors. Puppeteer is comprised of two components. The core, and the sensitivity analysis script. The core implements the programming model and runtime extensions. The sensitivity analysis scripts analyze annotated applications and extracts the sensitivities of the annotated kernels. 

## Requirements of Puppeteer:

Puppeteer on a x86_64 system has been successfully built with the following dependencies: 
- Ninja 1.9.0
- CMAKE 3.14
- gcc 8.3.1
- python 3.7.6
- cycler 0.11.0
- dill 0.3.5.1
- docopt 0.6.2
- importlib-metadata 4.12.0
- kaleido 0.2.1
- kiwisolver 1.4.3
- matplotlib 3.3.4
- multiprocess 0.70.13
- numpy 1.20.1
- pandas 1.2.2
- pathos 0.2.9
- Pillow 9.2.0
- pipreqs 0.4.11
- plotly 5.9.0
- pox 0.3.1
- ppft 1.7.6.5
- pyparsing 3.0.9
- python-dateutil 2.8.2
- pytz 2022.1
- PyYAML 6.0
- SALib 1.4.5
- scipy 1.7.3
- seaborn 0.11.1
- six 1.16.0
- tenacity 8.0.1
- yarg 0.1.9
- zipp 3.8.0


## Build Puppeteer Core components:

### Build Puppeteer Compiler Extensions

To build the Puppeteer compiler and runtime extensions please execute the following commands:

```bash
git clone git@github.com:koparasy/HPAC.git -b Puppeteer
cd HPAC
./setup.sh 'PREFIX' 'NUM THREADS' 
```

The 'PREFIX' argument defines where to install all the Puppeteer related binaries and executables. The 'NUM THREADS' parameter defines how many threads should the installation use. The installation script performs the following actions:

1. Configures, builds and installs clang/LLVM/OpenMP including the HPAC/Puppeteer extensions.
2. Configures, builds and installs the HPAC/Puppeteer runtime library. 
3. Creates a file, called 'puppet_env.sh', at the root of the project which should always be sourced before using HPAC.

The installation process can take a considerable amount of time. In our tested system, the complete 
installation took ~2 hours. 

## Replicate Analysis of Puppeteer (Approximate Computing Through the Lens of Uncertainty Quantification)

Please follow the instructions provided in [a relative link](approx/puppeteer/README.md)


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
