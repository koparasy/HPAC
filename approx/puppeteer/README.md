# Approximate Computing Through the Lens of Uncertainty Quantification

This directory contains the source code of the analyzed benchmarks, the approximate version of those benchmarks, the input files used, a script responsible to perform the sensitivity analysis and the raw data used to make the plots of the published work.

## Structure of directories

We are mainly interested in the following directories :
- data
- libraries 
- benchmarks
- errors

### The RAW data

The 'data' directory contains all the date we used in our published work. To visualize the data one can issue the following commands:

```bash
mkdir visualize
cd visualize
cmake -DWITH_VISUALIZE='On' ../ 
```

These commands will create all the plots of our work. "cmake" outputs the path to each plot.

### The 'libraries' directory

The directory contains source code files that are used to read/write data in a uniform manner and also creates a binary that can compute the global error function (Equation 4 in the manuscript). The library and the binary is used by the analysis script.

### The 'benchmarks' directory

The benchmarks directory contains all the source code of the analyzed benchmarks and also the approximate version of those benchmarks. The approximated versions filename contains the suffix '\_eval.cpp'.
Moreover it contains configuration files named with a suffix 'yaml.in'. The configuration files contain build commands and command line arguments. Finally, each directory contains a CMakeLists.txt file that builds all the required executables.

