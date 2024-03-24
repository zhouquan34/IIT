#ifndef BVSR_READIN
#define BVSR_READIN


#include "global.h"
#include "generic.h"
#include "lalg.h"
#include "model.h"

/*
g_n_cov is forced to be zero. do not consider covariates in this version.  
g_data, g_cov are k x g_n_sub matrices; each row is a variable.
read_* also set g_n_sub and g_n_snp.

missing value: ? or NA. 
g_MISS = -20000: thus no input value can be less than that.

All input data will be centered (g_data, g_cov, g_pheno).
Input files are either comma or whitespace delimited. 
Intercept 1 is never used due to centering. 

meang (remove 3 cols): 
  each row: SNP1 A T 1 0 2 0.5 1 ...

plink (remove 1 row, 6 cols): 
  row 1: header
  row 2: 1 SNP222  0 2000 A T 0 0 1 1 0 0  

matrix (remove 1 row, transpose):
  row 1: header (SNP names)
  row 2: 0 1 2  0 0 

cov (transpose): same format as matrix, but without header.

pheno: g_n_sub rows, each row 1 number (value of y).

*/

void read_meang(std::string);
void read_plink(std::string);
void read_matrix(std::string);
void read_pheno(std::string);
void read_true(std::string);
void read_omit(std::string);

#endif

