#ifndef BVSR_GLOBAL
#define BVSR_GLOBAL

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <map>
#include <sys/time.h>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>

extern const double NINF;
/********** data **********/
extern const int g_n_cov;
extern int g_n_sub;
extern int g_n_snp;
extern int g_missing;
extern int g_true;
extern float* g_data;
extern float* g_cov;
extern float** g_data_mat;
extern float** g_cov_mat;
extern double* g_pheno;
extern double* g_xty;
extern double g_yy;
extern double* g_cov_R;
extern std::vector<std::string> g_snp_names;
extern std::vector<double> g_snp_var;
extern std::vector<double> g_snp_mean;
extern std::vector<int> g_snp_map;
extern std::vector<int> g_true_gamma;
extern std::vector<double> g_true_beta;
extern std::vector<int> g_init_model;
extern int g_precomp;
extern double* g_xtx;

/************ IO ************/
extern std::string g_out_prefix;
extern std::ofstream LOG;
extern std::ofstream PATH;
extern std::ofstream GAMMA;
extern std::ofstream MOD;
extern std::ofstream OUT;
extern std::ofstream YHAT;
extern std::ofstream ERR;
extern std::map<std::string, int> g_paras;

/************** time *********/
extern int g_output_time;
extern std::map<std::string, double> g_time;
extern std::map<std::string, int> g_time_call;
extern struct timeval g_t_beg;
extern struct timeval g_t_end;

/*********** model *********/
extern int g_use_kappa;
extern double g_kappa;
//extern double g_gg;
extern double g_gg_exp;
extern double g_gg_const;
extern double g_sigma;
extern double g_min_pi;
extern double g_max_pi;
extern int g_pi_range_set;
extern int g_log_uniform_h;
extern int g_pi_prior_alpha;
extern int g_pi_prior_beta;
extern int g_max_model_size;
extern int g_min_model_size;
extern int g_log_uniform_h;
extern double g_jump_h2;
extern int g_pmin;
extern int g_pmax;
extern int g_pset;

/************* ridge **********/
extern double g_diag_c;
extern double g_icf_abs_tol;
extern int g_icf_min_p;
extern std::map<int, int> g_icf_iter;
extern int g_chol_call;

/*********** mcmc settings **********/
extern int g_max_jump;
extern double g_long_range;
extern double g_prop_add;
extern std::string g_last_bvsr;
extern int g_ordered_start;
extern int g_start_size;
extern int g_rb_thin;
extern double g_jump_h2;
extern int g_rw;
extern int g_skip;
extern int g_lbsq;
extern int g_approx_informed; 
extern double g_upper_add;
extern double g_lower_add;
extern double g_upper_del;
extern double g_lower_del;
extern int g_h_func;
extern double g_h_A;
extern double g_log_IS;

/********* mcmc **********/
extern int g_mcmc_i;
extern int g_count;
extern int g_rb_count;
extern int g_mcmc_iter;
extern int g_last_iter;
extern int g_mcmc_warm;
extern int g_mcmc_stay;
extern int g_continue;
extern int g_exact_bf;
extern std::vector<double> g_single_bf;
extern std::vector<double> g_single_h2;
extern std::vector<double> g_pip;
extern std::vector<double> g_beta_pos;
extern std::vector<double> g_beta_neg;
extern std::vector<double> g_rb_pip;
extern std::vector<double> g_rb_beta_pos;
extern std::vector<double> g_rb_beta_neg;
extern std::vector<int> g_single_ord;
extern std::vector<int> g_single_id;
extern std::vector<int> g_propose;
extern std::vector<int> g_accept;
extern std::map<int, double> g_model_size;

/********* gsl ********/
extern const gsl_rng_type* gsl_type;
extern gsl_rng* gsl_r;
extern int gsl_seed;

#endif
