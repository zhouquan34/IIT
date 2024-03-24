#ifndef BVSR_RIDGE
#define BVSR_RIDGE

#include "global.h"
#include "generic.h"
#include "lalg.h"
#include "chol.h"

/*
 calc_globals();
    if g_n_cov > 0, regress out confounding covariates.
    compute g_yy, g_xty. 

 calc_xtx();
    precompute part of xtx.

 filter_snp(sigma); 
    filter identical snps if filtering is turned on.
    compute single SNP bayes factor and order BF. 
    compute g_single_id, g_single_ord.  

 order: calc_globals => filter_snp => calc_xtx. 
*/

void precalc(void);

#endif

