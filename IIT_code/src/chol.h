#ifndef BVSR_CHOL
#define BVSR_CHOL

#include "xy.h"
#include "generic.h"
#include "global.h"

/*

Chol_add( new_p - last_p, last_p, new_XX, last_R, new_R);
  last_p is the dim of last_R; no initialization.
  
Chol_del( set_to_delete, last_p, last_R, new_R);	
  del must be sorted; no initialization.

Chol_init( new_p, new_XX, new_R);
  no initialization.

Chol_decomp_A(new_p, new_XX, new_R, vec_diag);	
  g_cov_R already computed; not used in the current version; no initialization.
  vec_diag has length new_p; thus, the first n_cov are not used.
  return log_determinant.

Chol_decomp_A(new_p, new_XX, new_R, diagonal);	
  g_cov_R already computed; no initialization.
  return log_determinant.

Chol_decomp_g_cov();
   use g_cov (cov data vector) and g_n_cov to compute g_cov_R.
   compute X_c^t X_c (which is not needed elsewhere). 

Chol_solve(p, R, z, b);		
   solve R^t R b = z using two steps: R^t b1 = z, and  R b = b1.
   solution saved in b; no initilization.
backwardk: R b = z; 

*/

void Chol_add (int, int, double*, double*, double*);
void Chol_del (std::vector<int>, int, double*, double*);	
void Chol_init (int, double*, double*);
double Chol_decomp_A (int, double*, double*, double);
void Chol_decomp_g_cov(void);
void Chol_solve (int, double*, double*, double*);
void Chol_solve_backward (int, double*, double*, double*);

#endif

