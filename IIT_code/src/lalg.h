#ifndef BVSR_LALG
#define BVSR_LALG

#include "global.h"
#include "generic.h"

#define ZERO 1e-12

/*

allocate1D(length, val_to_fill);

allocate2D(p);
  create a pointer to a p*2 matrix; that is, return m (length=2; each is a pointer to a p-vector).

copy1D(x_pointer, y_pointer, p);
  copy p elements of y to x.
  note that p may be less than the length of x, and we don't have to start at x[0]. 

mat_vec_mul(X_pointer, y_pointer, p, n, result_pointer, if_transpose);
  if_transpose = 0: X is p * n matrix, y is n-vector, result is p-vec.
  if_transpose = 1: X is p * n matrix; y is p-vector, result is n-vec.

vec_vec_mul(x_pointer, y_pointer, length_x);
  x, y same length; return x^t y.

center_rows(X_pointer, p, n)
  X is p * n.

center_row(x_pointer, n, if_calc_var);
  x is n-vector; if_calc_var = 1: variance is saved to g_snp_mean and g_snp_var. 
  Thus, center_rows or center_row is used only when reading the data. 

center_array(x_pointer, n);
  return sample variance of x. 

array_sst(x_pointer, n);
  return total sum of squares; i.e. sample variance * n. 

L2_norm(x_pointer, n);
  return L2_norm of x.

array_del(indices_to_del, last_p, last_v, new_v);
  no intialization. 

matrix_del(covs_to_del, last_p, last_n, last_M, new_M);
  new_M is (p - k) * n. 


*/

double* allocate1D (int);
float*  allocate1D_float (int);
double* allocate1D (int, int);
double** allocate2D (int);
float** allocate1D_float_pointer (int);

template <typename T> void free1D(T*);
void free1D_pointer (float**);
void free2D (double**);

template <typename T, typename U> void copy1D (T*, const U*, int);
void copy1D (float**, float**, int);

template <typename T, typename U> void mat_vec_mul (T**, U*, int, int, double*, int);
template <typename T, typename U> double vec_vec_mul (T*, U*, int);

template <typename T> void center_rows(T*, int , int);
double center_array(double*, int);
double L2_norm (double*, int);
double array_sst(double*, int);
void array_del (const std::vector<int>&, int, double*, double*);
template <typename T> void matrix_del(const std::vector<int>&, int, int, T*, T*);

#endif

