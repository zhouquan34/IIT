#ifndef BVSR_GEN
#define BVSR_GEN

#include "global.h"

void gsl_init(void);
void gsl_exit(void);
bool compare_pair(std::pair<int, double>, std::pair<int, double>);
    
int weighted_sample (std::vector<double>&);
int log_weighted_sample (std::vector<double>&, double&);
double log_sum_log(double, double);
void get_row_col (const std::string, int&, int&);

template <typename T> void print_vector(std::vector<T>);
template <typename T> void print_array(const T*, int);
template <typename T> void print_matrix(const T*, int, int);
void write_log_cmds(int, char**);


/*
  write_time(seconds);
  
  output_time();
     output time used by each individual operation. 

  record_time();
     typical use:  record_time(); ... some operation ... ;  record_time("some operation")
*/

void write_time(double);   
void output_time(void);
void record_time(std::string="");


#endif

