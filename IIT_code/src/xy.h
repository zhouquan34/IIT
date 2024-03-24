#ifndef BVSR_XY
#define BVSR_XY

#include "global.h"
#include "lalg.h"
#include "generic.h"

/*

No memory allocation. 

  X_init(covs, X); 
    Let k = length(covs). X dim will be (k + g_n_cov) * g_n_sub. 
  
  X_add(add_set, p_last, X_last, X_new);
    covs in add_set are appended to X.
  
  X_del(del_set, p_last, X_last, X_new);
    covs in del_set are deleted. For k in del_set, delete the (k + g_n_cov)-th row in X_last.
    del_set should be sorted.   
  
  Xy_init(X, p, Xy); 
    X is p * g_n_sub; compute Xy where y is g_pheno. Note that often p = msize = snp.size + g_n_cov. 

  Xy_add(add_set, p_last, Xy_last, Xy_new); 
    new elements are copied from g_xty (i.e. X y where X is the entire data) which is pre-computed.

  Xy_del(del_set, p_last, Xy_last, Xy_new); 
    del_set must be sorted. 

  XX_init(X, p, XX);
    compute XX = X^t X; diagonal elements (after g_n_cov) are added by g_diag_c for stability. 
    length of XX is 1 + 2 + ... + p = p(p+1)/2  (lower triangular). 
  
  XX_get(i, j);
    return X_i^t X_j 

g_xtx is a k * k matrix, where k = g_precomp + g_n_cov; it is also symmetric. 
  
  XX_init(set, X, p, XX); 
    set includes selected covariates; compute XX = X^t X. 
    this is different from XX_init(X, p, XX) in that we use g_xtx. 

  XX_add(add_set, p_last, X_new, XX_last, XX_new); 
  
  XX_add(set_new, add_set, p_last, X_new, XX_last, XX_new); 
    set_new is used so that we can directly use g_xtx; 

  XX_del(del_set, p_last, XX_last, XX_new);
    del_set does not need to be ordered.
*/

void X_init(const std::vector<int>&, float**);
void X_add(const std::vector<int>&, int, float**, float**);
void X_del(const std::vector<int>&, int, float**, float**);
void Xy_init(float**, int, double*);
void Xy_add(const std::vector<int>&, int, double*, double*);
void Xy_del(const std::vector<int>&, int, double*, double*);
void XX_init(float*, int, double*);
void XX_init(const std::vector<int>&, float**, int, double*);
void XX_add(const std::vector<int>&, const std::vector<int>&, int, float**, double*, double*);
void XX_del(const std::vector<int>&, int, double*, double*);
double XX_get(int, int);

#endif
