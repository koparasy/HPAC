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

### The 'errors' directory

The directory contains a description of the error domain of each application. The error domain is used by the sensitivity analysis script to sample error points and perform perturbation on the requested kernels. 

## Build benchmarks and perform simple executions.

As a first step we need to pull the input files we used for our analysis. To do so, please issue the following command:

```bash
git lfs pull
```

Next we need to configure our project using cmake. The $PUPPET_ROOT environment variable is set by the
'puppet\_env.sh' file.
```bash
mkdir $PROJECT_BUILD_DIR 
cd $PROJECT_BUILD_DIR 
cmake $PUPPET_ROOT
```

These commands configure our project and also update the configuration files (\*.yaml.in). 
You can set any path in the variable "$PROJECT\_BUILD\_DIR".
To build a specific benchmark (for example blackscholes) perform the following steps:
```bash
cd $PROJECT_BUILD_DIR
cd benchmarks/blackscholes
make
```

The command will create three binaries (blck\_accurate blck\_analysis  blck\_approx).
blck\_accurate does not contain any Puppeteer extensions. Whereas the blck\_analysis does. Finally, the blck\_approx binary executes is the approximate version of blackscholes. The rest of the benchmarks create binaries with similar suffixes. 

To run the 'accurate' blackscholes executable issue the following command:
```bash
./blck_accurate 1 $PUPPET_ROOT/inputs/random_input.bin accurate.out
```

The command executes the binary and stores the output in the file accurate. Respectively, one can execute the approximate version by issuing the following command:

```bash
./blck_approx 1 $PUPPET_ROOT/inputs/random_input.bin approx.out
```

Both execution list the execution time. To compute the error you can issue the command:

```bash
$PROJECT_BUILD_DIR/bin/quality -m RE -a ./accurate.out -t approx.out
```

The option '-m' RE instructs our framework to compute the equation 4 between the two files.

### Perform Trace mode using Puppeteer

To perform a trace execution of the analysis binary you can issue the following command.

```bash
export PETRUBATE_FILE=tmp.json
export PETRUBATE_TYPE=RECORD

./blck_analysis 1 $PUPPET_ROOT/inputs/random_input.bin test.out
```

```json
{
    "CNDF_1": {
        "NofXd1": [
            1,
            0.0,
            "double",
            1,
            0
        ]
    },
    "CNDF_2": {
        "NofXd2": [
            1,
            0.0,
            "double",
            1,
            0
        ]
    },
    "EXP": {
        "FutureValueX": [
            1,
            0.0,
            "double",
            1,
            0
        ]
    },
    "LOG_1": {
        "logValues": [
            1,
            0.0,
            "double",
            1,
            0
        ]
    },
    "SQRT_1": {
        "xSqrtTime": [
            1,
            0.0,
            "double",
            1,
            0
        ]
    }
}
```

The json file simillar to the one presented in Listing 4 of the original publication.

To perform an execution in which we inject an error using Equation 2 please modify tmp.json to:

```json
{
    "CNDF_1": {
        "NofXd1": [
            1,
            0.001,
            "double",
            2,
            0
        ]
    },
    "CNDF_2": {
        "NofXd2": [
            1,
            0.001,
            "double",
            2,
            0
        ]
    },
    "EXP": {
        "FutureValueX": [
            1,
            0.001,
            "double",
            2,
            0
        ]
    },
    "LOG_1": {
        "logValues": [
            1,
            0.001,
            "double",
            2,
            0
        ]
    },
    "SQRT_1": {
        "xSqrtTime": [
            1,
            0.001,
            "double",
            2,
            0
        ]
    }
}
```

This configuration instructs Puppeteer to inject an error of 0.001 on all output kernels.
To perform the execution issue the following commands:

```bash
export PETRUBATE_FILE=tmp.json
export PETRUBATE_TYPE=PETRUBATE

./blck_analysis 1 $PUPPET_ROOT/inputs/random_input.bin test.out

```

The erroneous output will we stored in the file 'test.out'. You can compare this output with the 'accurate.out' using the quality utility.
