#ifndef BVSR_ARGS
#define BVSR_ARGS

#define PARA_MG     1 
#define PARA_PLINK  2 
#define PARA_MAT    3
#define PARA_PHENO  4
#define PARA_OUT    5
#define PARA_ITER   6
#define PARA_BURN   7
#define PARA_RB     8
#define PARA_START  9
#define PARA_CHOL   10
#define PARA_TIME   11
#define PARA_LONG   12
#define PARA_LEX    13
#define PARA_BVSR   14
#define PARA_ICFE   15
#define PARA_LD	    16
#define PARA_GEOM   17
#define PARA_MINPI  18
#define PARA_MAXPI  19
#define PARA_MINH2  20
#define PARA_MAXH2  21
#define PARA_MAXJ   22
#define PARA_ICFL   23
#define PARA_HELP   24
#define PARA_MAN    25
#define PARA_BESTH  26
#define PARA_LOGH   27
#define PARA_HJUMP  28
#define PARA_UADD   29
#define PARA_UH2    30
#define PARA_NOF    31
#define PARA_PI_A   32
#define PARA_PI_B   33
#define PARA_SEED   34
#define PARA_COV    35 
#define PARA_PREC   36
#define PARA_SBFP   37
#define PARA_SBF    38
#define PARA_SIZEA  39
#define PARA_PI_A2  40
#define PARA_PI_B2  41
#define PARA_MAXS   42
#define PARA_MINPI2 43
#define PARA_MAXPI2 44
#define PARA_NSNP   45
#define PARA_NTRY   46
#define PARA_EXACT  47
#define PARA_MINS   48
#define PARA_PMAX   49
#define PARA_PMIN   50
#define PARA_RB1    51
#define PARA_NOT    52
#define PARA_NOT2   53
#define PARA_YHAT   54
#define PARA_GAMMA  55
#define PARA_MAXADD 56
#define PARA_MINADD 57
#define PARA_MAXDEL 58
#define PARA_MINDEL 59
#define PARA_G      60
#define PARA_RW     61
#define PARA_KAPPA  62
#define PARA_TRUE   63
#define PARA_SFIX   64
#define PARA_INIT   65
#define PARA_AI     66
#define PARA_OMIT   67
#define PARA_LB     68
#define PARA_H      69
#define PARA_H_A    70

#include "global.h"
#include "generic.h"
#include "readin.h"
#include "model.h"

void print_help(void);
void print_manual(void);
void print_continue_help(void);
void read_args(int, char**);

#endif



