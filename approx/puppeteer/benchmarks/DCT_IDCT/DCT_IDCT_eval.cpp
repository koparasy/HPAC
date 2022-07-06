#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <executionTime.h>
#include <datautils.hpp>
#include <approx.h>
#define N 512
#define DOUBLE 0

double COS[8][8], C[8];
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const char *coeff_name[8][8] = {
    "Coeff_0","Coeff_1",  "Coeff_2" ,"Coeff_3" , "Coeff_4","Coeff_5" ,"Coeff_6","Coeff_7", 
    "Coeff_8","Coeff_9",  "Coeff_10","Coeff_11","Coeff_12","Coeff_13","Coeff_14","Coeff_15", 
    "Coeff_16","Coeff_17","Coeff_18","Coeff_19","Coeff_20","Coeff_21","Coeff_22","Coeff_23", 
    "Coeff_24","Coeff_25","Coeff_26","Coeff_27","Coeff_28","Coeff_29","Coeff_30","Coeff_31", 
    "Coeff_32","Coeff_33","Coeff_34","Coeff_35","Coeff_36","Coeff_37","Coeff_38","Coeff_39", 
    "Coeff_40","Coeff_41","Coeff_42","Coeff_43","Coeff_44","Coeff_45","Coeff_46","Coeff_47", 
    "Coeff_48","Coeff_49","Coeff_50","Coeff_51","Coeff_52","Coeff_53","Coeff_54","Coeff_55", 
    "Coeff_56","Coeff_57","Coeff_58","Coeff_59","Coeff_60","Coeff_61","Coeff_62","Coeff_63" 
  };

#ifdef WITH_QUANT 
const int sensitivity[] = { 0,8,1,16,9,24,17,2,32,10,25,48,35,42,28,34,49,5,18,3,41,20,33,56,40,11,21,26,13,27,19,12,4,43,6,36,50,14,15,57,7,22,23,58,44,29,51,30,39,52,31,37,59,55,38,45,60,62,63,54,46,47,61,53 }; 
#else
const int sensitivity[] = {0,8,1,16,9,24,17,2,10,32,25,18,3,33,11,40,26,19,4,41,12,48,34,27,42,56,20,49,5,28,35,57,13,21,50,36,43,6,29,14,58,22,7,51,37,44,15,30,59,52,45,38,23,31,53,60,46,61,39,47,54,62,55,63};
#endif

unsigned char** malloc2D(int numRows, int numCols){
    long len = sizeof(void*) * numRows + sizeof(unsigned char) * numRows * numCols;
    unsigned char **ptr = (unsigned char **)malloc(len);
    unsigned char *start;
    start = (unsigned char*)(ptr + numRows);
    // for loop to point rows pointer to appropriate location in 2D array
    for(int i = 0; i < numRows; i++)
        ptr[i] = (start + numCols * i);
    return ptr;
}

double** malloc2DD(int numRows, int numCols){
    long len = sizeof(void*) * numRows + sizeof(double) * numRows * numCols;
    double **ptr = (double **)malloc(len);
    double *start;
    start = (double*)(ptr + numRows);
 
    // for loop to point rows pointer to appropriate location in 2D array
    for(int i = 0; i < numRows; i++)
        ptr[i] = (start + numCols * i);
    return ptr;
}


void init() {
  int i, j;
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++)
      COS[i][j] = cos((2 * i + 1) * j * acos(-1) / 16.0);
    if (i) C[i] = 1;
    else C[i] = 1 / sqrt(2);
  }
}

void DCT(unsigned char **img, double **dct, int ROWS, int COLS, int lsk) {
  int y_blocks = ROWS/8;
  int x_blocks = COLS/8;
  int r, c, i, j, x, y;
  char label[100];
  for (r = 0; r < y_blocks; r++)
    for (c = 0; c < x_blocks; c++)
      for (i = 0; i < 8; i++)
        for (j = 0; j < 8; j++) {
          if ( sensitivity[i*8+j] < lsk ){
              double sum = 0;
              int end_x = 8;
              for (x = 0; x < end_x; x++){
                int end_y = 8;
                for (y = 0; y < end_y; y++){
                  sum += (img[r * 8 + x][c * 8 + y] - 128) * COS[x][i] * COS[y][j];
                }
              }
              sum *= C[i] * C[j] * 0.25;
              dct[r * 8 + i][c * 8 + j] = sum;
          }
      }
}

void IDCT(double **dct, double **idct, int ROWS, int COLS) {
  int y_blocks = ROWS/8;
  int x_blocks = COLS/8;
  int r, c, i, j, x, y;
  for (r = 0; r < y_blocks; r++)
    for (c = 0; c < x_blocks; c++)
      for (i = 0; i < 8; i++)
        for (j = 0; j < 8; j++) {
          double sum = 0;
          for (x = 0; x < 8; x++)
            for (y = 0; y < 8; y++)
              sum += C[x] * C[y] * dct[r * 8 + x][c * 8 + y] * COS[i][x] * COS[j][y];
          sum *= 0.25;
          sum += 128;
          idct[r * 8 + i][c * 8 + j] = sum;
      }
}

void quantization(double **dct, int ROWS, int COLS) {
  int table[8][8] = {
    16, 11, 10, 16, 24, 40, 51, 61,
    12, 12, 14, 19, 26, 58, 60, 55,
    14, 13, 16, 24, 40, 57, 69, 56,
    14, 17, 22, 29, 51, 87, 80, 82,
    18, 22, 37, 56, 68, 109, 103, 77,
    24, 35, 55, 64, 81, 104, 113, 92,
    49, 64, 78, 87, 103, 121, 120, 101,
    72, 92, 95, 98, 112, 100, 103, 99
  };
  int r, c, i, j;
  int y_blocks = ROWS/8;
  int x_blocks = COLS/8;
  for (r = 0; r < y_blocks; r++)
    for (c = 0; c < x_blocks; c++)
      for (i = 0; i < 8; i++)
        for (j = 0; j < 8; j++) {
          dct[r * 8 + i][c * 8 + j] = round(dct[r * 8 + i][c * 8 + j] / table[i][j]);
          dct[r * 8 + i][c * 8 + j] = dct[r * 8 + i][c * 8 + j] * table[i][j];
        }
}

void MSE(unsigned char **pic, double **idct, int COLS, int ROWS) {
  double MSE = 0;
  int r, c;
  for (r = 0; r < ROWS; r++)
    for (c = 0; c < COLS; c++) {
      MSE += (pic[r][c] - idct[r][c]) * (pic[r][c] - idct[r][c]);
    }
  printf("MSE is %g\n", MSE);
  MSE /= (COLS*ROWS);
  double PSNR = 10 * log10(255.0 * 255.0 / MSE);
  printf("%.2lf\n", PSNR);
}

int main(int argc, char *argv[]) {
  int COLS = atoi(argv[1]); 
  int ROWS = atoi(argv[2]); 
  char *in_image = argv[3]; 
  char *out_image = argv[4];
  int lsk = atoi(argv[5]);
  double **dct;
  double **idct;
  int r, c;
  long numBytes = COLS*ROWS;
  unsigned char **in_img = malloc2D(ROWS, COLS); 
  unsigned char **out_img = malloc2D(ROWS, COLS); 
  dct = malloc2DD(ROWS, COLS); 
  idct = malloc2DD(ROWS, COLS); 
  FILE *fd = fopen(in_image,"r");
  fread(&in_img[0][0], sizeof(unsigned char), ROWS*COLS, fd);
  fclose(fd);
  double **dimg = malloc2DD(ROWS, COLS); 
  for(int i = 0; i < ROWS; i ++){
    for ( int j = 0; j < COLS; j++){
      dimg[i][j] = (double) in_img[i][j];
    }
  }
  writeData(&dimg[0][0], ROWS*COLS, DOUBLE, "d_in_img.raw");

  init();
  startMeasure();
  for (int i = 0; i < 20; i++){
    DCT(in_img, dct, ROWS, COLS, 64-lsk);
  }
  stopMeasure();
  quantization(dct, ROWS, COLS);
  IDCT(dct, idct, ROWS, COLS);
  writeData(&idct[0][0], ROWS*COLS, DOUBLE, out_image);
  return 0;
}
