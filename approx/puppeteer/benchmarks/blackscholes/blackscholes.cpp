// Copyright (c) 2007 Intel Corp.

// Black-Scholes
// Analytical method for calculating European Options
//
// 
// Reference Source: Options, Futures, and Other Derivatives, 3rd Edition, Prentice 
// Hall, John C. Hull,

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <executionTime.h>
#include <datautils.hpp>
#include <approx.h>
#include <omp.h>

//Precision to use for calculations
#define fptype double 

#define DOUBLE 0
#define FLOAT 1
#define INT 2


fptype *prices;
size_t numOptions;

int    * otype;
fptype * sptprice;
fptype * strike;
fptype * rate;
fptype * volatility;
fptype * otime;
int numError = 0;
int nThreads;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Cumulative Normal Distribution Function
// See Hull, Section 11.8, P.243-244
//
//CONSTANTS
const fptype inv_sqrt_2xPI=0.39894228040143270286f;
const fptype zero = 0.0;
const fptype half = 0.5;
const fptype const1=0.2316419;
const fptype one=1.0;
const fptype const2=0.319381530;
const fptype const3=0.356563782;
const fptype const4=1.781477937;
const fptype const5=1.821255978;
const fptype const6=1.330274429;

fptype CNDF ( fptype InputX ) 
{
    int sign;

    fptype OutputX;
    fptype xInput;
    fptype xNPrimeofX;
    fptype expValues;
    fptype xK2;
    fptype xK2_2, xK2_3;
    fptype xK2_4, xK2_5;
    fptype xLocal, xLocal_1;
    fptype xLocal_2, xLocal_3;
    fptype temp;

    // Check for negative value of InputX
    if (InputX < zero) {
        InputX = -InputX;
        sign = 1;
    } else 
        sign = 0;

    xInput = InputX;

    // Compute NPrimeX term common to both four & six decimal accuracy calcs
    temp = -half * InputX * InputX;

    expValues = expf(temp);

    xNPrimeofX = expValues;
    xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;

    xK2 = const1* xInput;
    xK2 = one + xK2;
    xK2 = one / xK2;
    xK2_2 = xK2 * xK2;
    xK2_3 = xK2_2 * xK2;
    xK2_4 = xK2_3 * xK2;
    xK2_5 = xK2_4 * xK2;

    xLocal_1 = xK2 * const2;
    xLocal_2 = xK2_2 * (-const3);
    xLocal_3 = xK2_3 * const4;
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_4 * (-const5);
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_5 * const6;
    xLocal_2 = xLocal_2 + xLocal_3;

    xLocal_1 = xLocal_2 + xLocal_1;
    xLocal   = xLocal_1 * xNPrimeofX;
    xLocal   = one - xLocal;

    OutputX  = xLocal;

    if (sign) {
        OutputX = one - OutputX;
    }

    return OutputX;
}

fptype BlkSchlsEqEuroNoDiv( fptype sptprice,
        fptype strike, fptype rate, fptype volatility,
        fptype time, int otype, float timet )
{
    fptype OptionPrice;

    // local private working variables for the calculation
    fptype xStockPrice;
    fptype xStrikePrice;
    fptype xRiskFreeRate;
    fptype xVolatility;

    fptype xTime;
    fptype xSqrtTime;

    fptype logValues;
    fptype xLogTerm;
    fptype xD1;
    fptype xD2;
    fptype xPowerTerm;
    fptype xDen;
    fptype d1;
    fptype d2;
    fptype FutureValueX;
    fptype NofXd1;
    fptype NofXd2;
    fptype NegNofXd1;
    fptype NegNofXd2;    
    fptype temp;

    xStockPrice = sptprice;
    xStrikePrice = strike;
    xRiskFreeRate = rate;
    xVolatility = volatility;

    xTime = time;

#pragma approx label("SQRT_1") petrubate(out) out(xSqrtTime)
    xSqrtTime = sqrtf(xTime);


#pragma approx label("LOG_1") petrubate(out) out(logValues)
    logValues = logf( sptprice / strike );

    xLogTerm = logValues;

    xPowerTerm = xVolatility * xVolatility;
    xPowerTerm = xPowerTerm * half;

    xD1 = xRiskFreeRate + xPowerTerm;
    xD1 = xD1 * xTime;
    xD1 = xD1 + xLogTerm;

    xDen = xVolatility * xSqrtTime;
    xD1 = xD1 / xDen;
    xD2 = xD1 -  xDen;

    d1 = xD1;
    d2 = xD2;

#pragma approx label("CNDF_1") petrubate(out) out(NofXd1)
    NofXd1 = CNDF( d1 );

#pragma approx label("CNDF_2") petrubate(out) out(NofXd2)
    NofXd2 = CNDF( d2 );

    temp = -(rate*time);

#pragma approx label("EXP") petrubate(out) out(FutureValueX)
    FutureValueX =  ( expf( temp  ) );

    FutureValueX *=strike;

    if (otype == 0) {
        OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
    } else {
        NegNofXd1 = (one - NofXd1);
        NegNofXd2 = (one - NofXd2);
        OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
    }

    return OptionPrice;
}


int bs_thread(void *tid_ptr) {
    int i, j,k;
    fptype price;
    fptype priceDelta;
    int end = numOptions;
#pragma omp parallel for firstprivate(end)
    for (i=0; i<end; i++) {
        prices[i] = BlkSchlsEqEuroNoDiv( sptprice[i], strike[i],
                rate[i], volatility[i], otime[i],
                otype[i], 0);

    }

    return 0;
}

int main (int argc, char **argv)
{
    FILE *file;
    int i;
    int loopnum;
    int rv;
    struct timeval start, end; 

    // start timer. 


    printf("PARSEC Benchmark Suite\n");
    fflush(NULL);
    if (argc != 4)
    {
        printf("Usage:\n\t%s <nthreads> <inputFile> <outputFile>\n", argv[0]);
        exit(1);
    }
    nThreads = atoi(argv[1]);
    char *inputFile = argv[2];
    char *outputFile = argv[3];

    //Read input data from file
    file = fopen(inputFile, "rb");
    if(file == NULL) {
        printf("ERROR: Unable to open file `%s'.\n", inputFile);
        exit(1);
    }

    if(nThreads != 1) {
        printf("Error: <nthreads> must be 1 (serial version)\n");
        exit(1);
    }

#define PAD 256
#define LINESIZE 64
    readData(file,&otype, &numOptions);  
    readData(file,&sptprice, &numOptions);  
    readData(file,&strike, &numOptions);  
    readData(file,&rate, &numOptions);  
    readData(file,&volatility, &numOptions);  
    readData(file,&otime, &numOptions);  
    prices = (fptype*) malloc(sizeof(fptype)*numOptions);


    int tid=0;

    startMeasure();
    bs_thread(&tid);
    stopMeasure();

    //Write prices to output file
    writeData(prices, numOptions, DOUBLE, outputFile);
    free(sptprice);
    free(strike);
    free(rate);
    free(volatility);
    free(otime);
    free(otype);
    free(prices);

    return 0;
}

