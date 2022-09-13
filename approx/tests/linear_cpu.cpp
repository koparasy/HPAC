#include <iostream>
#include <cstdlib>

using namespace std;

void linear_model(double *x, double *y, double alpha, double beta, size_t elements){
#pragma approx in(x[0:elements]) out(y[0:elements])
#pragma omp parallel for firstprivate(alpha, beta) 
  for( int i = 0 ; i < elements; i++){
    y[i] = alpha*x[i] + beta;
  }
}

/* Programming model for approximations using HPAC on CPU */
template <class T>
void linear_model_cpu(T *x, T *y, T alpha, T beta, size_t elements){
  /* Challenges on offline training I need to somehow define that the
   * underlying operation is computing 'element' values by calling the same
   * computation multiple times. In the finest granularity the underlying 
   * computation computes 1 value. 
   * So that in the offline data base I give the following structure to the user:
    X   |   Y
    ----------
    x1  |   y1
    x2  |   y2
    ....
    xn  |   yn

   * This is important as I want to create a ML model that uses the same 
   * primitive input/output shape. 
   * A counter example is the following:
 '''  
  void linear_model_cpu(T *x, T *y, T *alpha, T *beta, size_t elements, int shape){
#pragma approx in(x[0:elements]) out(y[0:elements])
#pragma omp parallel for firstprivate(alpha, beta) 
    for( int i = 0 ; i < elements/shape; i++){
      for ( int j = 0; j < shape; j++) {
        y[i*shape+j] = alpha[j]*x[i*shape+j] + beta[j];
      }
    }
  '''
  * In this example all the data used are stil x[0:elements], y[0:elements]
  * However, the primitive function uses 'shape' elements to compute one fundamental
  * output value. So, in the finest granularity the underlying computation computes
  * 4 values. So, in the respective offline data base I need to give the following
  * structure to the user:
            X           |           Y
    -------------------------------------------
    x11, x12, x13, x14  |   y11, y12, y13, y14
    x21, x22, x23, x24  |   y11, y12, y13, y14
    ...                     
    xn1, xn2, xn3, xn4  |   y11, y12, y13, y14
  *
  
  How would we support for example stencil primitive operations? We need 
  to provide some mapping syntax that maps inputs into outputs.
*/
#pragma approx in(x[0:elements]) out(y[0:elements])
#pragma omp parallel for firstprivate(alpha, beta) 
  for( int i = 0 ; i < elements; i++){
    y[i] = alpha*x[i] + beta;
  }
}

int main(int argc, char *argv[]){
  if (argc != 2){
    cout << "Wrong command line\n";
    cout << argv[0] << " number of elements" << "\n";
    return -1;
  }
    
  size_t elements = atol(argv[1]);
  double *x = new double[elements];
  double *y = new double[elements];

  double alpha = 0.5;
  double beta = 1.5;

  for (int i = 0; i < elements; i++){
    x[i] = rand() /static_cast<double>( RAND_MAX ); 
  }

  linear_model(x, y, alpha, beta, elements);

  for (int i = 0; i < elements ; i++){
    cout << i << " " << x[i] << " " << y[i] << "\n";
  }

  return 0;
}

