#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cctype>
#include <executionTime.h>
#include <datautils.hpp>
#include <approx.h>
#include <omp.h>


#define DOUBLE 0
#define FLOAT 1
#define INT 2
#define LONG 3
#include <cmath>
using std::fabs;

const int max_external = 100000;
const int max_num_messages = 500;
const int max_num_neighbors = max_num_messages;

char   *HPC_title;
int HPC_start_row;
int HPC_stop_row;
int HPC_total_nrow;
long long HPC_total_nnz;
int HPC_local_nrow;
int HPC_local_ncol;  // Must be defined in make_local_matrix
int HPC_local_nnz;
int  * HPC_nnz_in_row;
double ** HPC_ptr_to_vals_in_row;
float**HPC_ptr_to_vals_in_row_single; 
int ** HPC_ptr_to_inds_in_row;
double ** HPC_ptr_to_diags;
double *HPC_list_of_vals;   //needed for cleaning up memory
float *HPC_list_of_vals_single;   //needed for cleaning up memory
int *HPC_list_of_inds;      //needed for cleaning up memory


template <typename T, typename U> 
void cast(const int n, T *__restrict x, U *__restrict__ y){
  int i;
#pragma omp parallel for firstprivate(n) 
    for ( i=0; i<n; i++) {
      y[i] = x[i];
    }
}

template <typename T, typename U, typename M> 
int ddot (const int n, const T* const x, const U* const y, 
    M* const result, double& time_allreduce)
{  
  U local_result = 0.0;
  int i;
  if ((void*) y== (void *) x){
#pragma omp parallel for reduction (+:local_result) firstprivate(n) 
    for ( i=0; i<n; i++) {
      local_result = local_result+ x[i]*x[i];
    }
  }
  else{
#pragma omp parallel for reduction (+:local_result)
    for (i=0; i<n; i++){
      local_result = local_result + x[i]*y[i];
    }
  }

  *result = local_result;
  return(0);
}

//waxpby.cpp
template <typename T, typename U, typename M> 
int waxpby (const int n, const T alpha, const T* const x, 
    const U beta, const U * const y, 
    M* const w)
{  
  int i;
  if (alpha==1.0) {
#pragma omp parallel for firstprivate(beta,n)
    for (i=0; i<n; i++) {
      w[i] = x[i] + beta * y[i];
    }
  }
  else if(beta==1.0) {
#pragma omp parallel for firstprivate(alpha,n)
    for (i=0; i<n; i++) {
      w[i] = alpha * x[i] + y[i];
    }
  }
  else {
#pragma omp parallel for firstprivate(alpha,beta,n)
    for (i=0; i<n; i++) {
      w[i] = alpha * x[i] + beta * y[i];
    }
  }
  return(0);
}

int compute_residual(const int n, const double * const v1, 
    const double * const v2, double * const residual)
{
  double local_residual = 0.0;
  int i;
  for (i=0; i<n; i++) {
    if ( *(long *) &v1[i] != *(long *) &v2[i]){
      double diff = fabs(v1[i] - v2[i]);
      if (diff > local_residual) local_residual = diff;
    }
  }
  *residual = local_residual;
  return(0);
}

template <typename T, typename U> 
int HPC_sparsemv_s(const T* const x, U* const y)
{
  const int nrow = HPC_local_nrow;
  int i;
#pragma omp parallel for firstprivate(nrow)
  for (i=0; i< nrow; i++)
  {
    int j;
    U sum = 0.0;
    const int    * const cur_inds = 
      (const int    * const) HPC_ptr_to_inds_in_row[i];

    const int cur_nnz = (const int) HPC_nnz_in_row[i];

    for (j=0; j< cur_nnz; j++)
      sum += HPC_ptr_to_vals_in_row_single[i][j]*x[cur_inds[j]];
    y[i] = sum;
  }
  return(0);
}



template <typename T, typename U> 
int HPC_sparsemv(const T* const x, U* const y)
{
  const int nrow = HPC_local_nrow;
  int i;
#pragma omp parallel for firstprivate(nrow)
  for (i=0; i< nrow; i++)
  {
    int j;
    U sum = 0.0;
    const int    * const cur_inds = 
      (const int    * const) HPC_ptr_to_inds_in_row[i];

    const int cur_nnz = (const int) HPC_nnz_in_row[i];

    for (j=0; j< cur_nnz; j++)
      sum += HPC_ptr_to_vals_in_row[i][j]*x[cur_inds[j]];
    y[i] = sum;
  }
  return(0);
}


//generate_matrix.cpp
//
void generate_matrix(int nx, int ny, int nz, double **x, double **b, double **xexact)
{
  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  HPC_title = 0;
  bool use_7pt_stencil = false;
  int local_nrow = nx*ny*nz; // This is the size of our subblock
  int local_nnz = 27*local_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)
  int total_nrow = local_nrow*size; // Total number of grid points in mesh
  long long total_nnz = 27* (long long) total_nrow; // Approximately 27 nonzeros per row (except for boundary nodes)
  int start_row = local_nrow*rank; // Each processor gets a section of a chimney stack domain
  int stop_row = start_row+local_nrow-1;
  // Allocate arrays that are of length local_nrow
  int tmp_size = local_nrow;
  HPC_nnz_in_row = new int[local_nrow];
  HPC_ptr_to_vals_in_row = new double*[local_nrow];
  HPC_ptr_to_vals_in_row_single = new float*[local_nrow];
  HPC_ptr_to_inds_in_row = new int*[local_nrow];
  HPC_ptr_to_diags       = new double*[local_nrow];

  *x =  new double[local_nrow];
  *b =  new double[local_nrow];
  *xexact =  new double[local_nrow];

  HPC_list_of_vals = new double[local_nnz];
  HPC_list_of_vals_single = new float[local_nnz];
  HPC_list_of_inds = new int   [local_nnz];

  double* curvalptr = HPC_list_of_vals;
  float * curvalptr_s = HPC_list_of_vals_single;
  int* curindptr = HPC_list_of_inds;

  long long nnzglobal = 0;
  int iz;
  int iy;
  int ix;
  for (iz=0; iz<nz; iz++) {
    for (iy=0; iy<ny; iy++) {
      for (ix=0; ix<nx; ix++) {
        int curlocalrow = iz*nx*ny+iy*nx+ix;
        int currow = start_row+iz*nx*ny+iy*nx+ix;
        int nnzrow = 0;
        HPC_ptr_to_vals_in_row[curlocalrow] = curvalptr;
        HPC_ptr_to_vals_in_row_single[curlocalrow] = curvalptr_s;
        HPC_ptr_to_inds_in_row[curlocalrow] = curindptr;
        int sx, sy, sz;
        for (sz=-1; sz<=1; sz++) {
          for (sy=-1; sy<=1; sy++) {
            for (sx=-1; sx<=1; sx++) {
              int curcol = currow+sz*nx*ny+sy*nx+sx;
              if ((ix+sx>=0) && (ix+sx<nx) && (iy+sy>=0) && (iy+sy<ny) && (curcol>=0 && curcol<total_nrow)) {
                if (!use_7pt_stencil || (sz*sz+sy*sy+sx*sx<=1)) { // This logic will skip over point that are not part of a 7-pt stencil
                  if (curcol==currow) {
                    HPC_ptr_to_diags[curlocalrow] = curvalptr;
                    *curvalptr++ = 27.0;
                    *curvalptr_s++ = 27.0;
                  }
                  else {
                    *curvalptr++ = -1.0;
                    *curvalptr_s++ = -1.0;
                  }
                  *curindptr++ = curcol;
                  nnzrow++;
                } 
              }
            } // end sx loop
          } // end sy loop
        } // end sz loop
        HPC_nnz_in_row[curlocalrow] = nnzrow;
        nnzglobal += nnzrow;
        (*x)[curlocalrow] = 0.0;
        (*b)[curlocalrow] = 27.0 - ((double) (nnzrow-1));
        (*xexact)[curlocalrow] = 1.0;
      } // end ix loop
    } // end iy loop
  } // end iz loop  
  HPC_start_row = start_row ; 
  HPC_stop_row = stop_row;
  HPC_total_nrow = total_nrow;
  HPC_total_nnz = total_nnz;
  HPC_local_nrow = local_nrow;
  HPC_local_ncol = local_nrow;
  HPC_local_nnz = local_nnz;
  return;
}

int HPCCG(const double * const b, double * const x,
    const int max_iter, const double tolerance, int &niters, double & normr,
    double * times, int num_splits, double *bounds, char **names)
{
  double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;
  int nrow = HPC_local_nrow;
  int ncol = HPC_local_ncol;
  float *rS =  new float[nrow];
  double *rD =  new double [nrow];
  float *p =  new float[ncol];
  double * ApD = new double[nrow];
  int status = DOUBLE;
  float *ApS = new float[nrow];

  normr = 0.0;
  double rtrans = 0.0;
  double oldrtrans = 0.0;

  int rank = 0; // Serial case (not using MPI)

  int print_freq = 5;//max_iter/10; 
  if (print_freq>50) print_freq=50;
  if (print_freq<1)  print_freq=1;

  // p is of length ncols, copy x to p for sparse MV operation
  waxpby(nrow, 1.0, x, 0.0, x, p); 
  HPC_sparsemv( p, ApD); 
  waxpby(nrow, 1.0, b, -1.0, ApD, rD);
  ddot(nrow, rD, rD, &rtrans, t4);
  normr = sqrt(rtrans);
  int k;
  int bin = 0;
  for (int b = 0; b < num_splits; b++){
    if ( normr > bounds[b] ){
      bin = b;
      break;
    }
  }
  for(k=1; k<max_iter && normr > tolerance; k++ )
  {

    if (k == 1)
    {
      waxpby(nrow, 1.0, rD, 0.0, rD, p);
    }
    else
    {
      oldrtrans = rtrans;
      //This is par_2
      if (status == DOUBLE){
        ddot (nrow, rD, rD, &rtrans, t4); 
        float beta = rtrans/oldrtrans;
        waxpby (nrow, 1.0, rD, beta, p, p); 
      }
      else{
        ddot (nrow, rS, rS, &rtrans, t4); 
        float beta = rtrans/oldrtrans;
        waxpby (nrow, 1.0f, rS, beta, p, p); 
      }

      //This is par_5
    }
    normr = sqrt(rtrans);

    for (int b = 0; b < num_splits; b++){
      if ( normr > bounds[b] ){
        if ( b != bin && b != 0){
          if (status == DOUBLE ){
            cast(nrow, ApD, ApS); 
            cast(nrow, rD, rS); 
            status = FLOAT;
          }
        }
        bin = b;
        break;
      }
    }


    if (rank==0 && (k%print_freq == 0 || k+1 == max_iter)){
      std::cout << "Iteration = "<< k << "   Residual = "<< normr << " BIN = " << bin << std::endl;
    }
    if (status == DOUBLE)
      HPC_sparsemv(p, ApD); // 2*nnz ops
    else 
      HPC_sparsemv_s(p, ApS); // 2*nnz ops
    float alpha = 0.0;
    //This is par(1)
    if (status == DOUBLE)
      ddot(nrow, p, ApD, &alpha, t4);
    else
      ddot(nrow, p, ApS, &alpha, t4);

    alpha = rtrans/alpha;
    waxpby(nrow, 1.0, x, alpha, p, x);// 2*nrow ops
    if (status == DOUBLE)
      waxpby(nrow, 1.0, rD, (double)-alpha, ApD, rD);  
    else
      waxpby(nrow, 1.0f, rS, -alpha, ApS, rS);  
    niters = k;
  }
  printf("I executed %d Iters %g\n", niters, normr);

  delete[] p;
  delete[] ApD;
  delete[] ApS;
  delete[] rD;
  delete[] rS;
  return(0);
}

void frexp10(double arg, int * exp)
{
  *exp =  (int) std::log2(std::fabs(arg));
  printf("Log 2 of %g is %g\n", arg, *exp);
}

int main(int argc, char *argv[])
{
  int niters = 0;
  double normr = 0.0;
  int max_iter = atoi(argv[6]);
  double tolerance = atof(argv[5]); 
  int num_splits = atoi(argv[7]);
  printf("Tolerance is %g iters are %d\n", tolerance, max_iter);
  int exponent;
  frexp10(tolerance, &exponent);
  printf("%d %g\n", exponent, tolerance); 

  double *bounds = (double *) malloc (sizeof(double)*num_splits);
  double step = (double) exponent / (double) num_splits; 
  for (int i = 1; i <= num_splits; i++){
    double tmp = step * i;
    bounds[i-1] = pow(2, tmp);
  }

  char **labels = (char**) malloc (sizeof(char*)*num_splits);
  for (int i = 0; i < num_splits; i++){
    labels[i] = (char*) malloc (sizeof(char)*100);
    snprintf(labels[i], 100, "Section"); 
  }

  double *x, *b, *xexact;
  double norm, d;
  int ierr = 0;
  int i, j;
  int ione = 1;
  double times[7];
  double t6 = 0.0;
  int nx,ny,nz;
  int size = 1; // Serial case (not using MPI)
  int rank = 0; 
  char *quality;
  nx = atoi(argv[1]);
  ny = atoi(argv[2]);
  nz = atoi(argv[3]);
  quality = argv[4];
  generate_matrix(nx, ny, nz, &x, &b, &xexact);
  startMeasure();
  ierr = HPCCG(  b, x, max_iter, tolerance, niters, normr, times, num_splits, bounds, labels);
  stopMeasure();
  size_t val = nx*ny*nz;
  double actualResidual;
  compute_residual(nx*ny*nz,x, xexact , &actualResidual);
  printf("Actual Error is %g\n", actualResidual);
  writeData(x,nx*ny*nz, DOUBLE,quality);   
  delete [] x;
  delete [] b;
  delete [] xexact;
  delete [] HPC_nnz_in_row;
  delete [] HPC_ptr_to_vals_in_row;
  delete [] HPC_ptr_to_inds_in_row;
  delete [] HPC_ptr_to_diags;
  delete [] HPC_list_of_vals;
  delete [] HPC_list_of_inds;
  return 0 ;
}
