#ifndef BVSR_OUTPUT
#define BVSR_OUTPUT

#include "generic.h"
#include "global.h"
#include "lalg.h"

void open_outputs(std::string);
void close_outputs(void);
void mcmc_output(void);
void mcmc_mean(void);

#endif

