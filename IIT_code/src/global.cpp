#include "global.h"

using namespace std;
const double NINF = -1e15;

/*******  data *******/
const int g_n_cov = 0; // cannot be changed
int g_n_sub = 0;
int g_n_snp = 0;
int g_missing = 0;
float* g_data = NULL;
float* g_cov = NULL;
float** g_data_mat = NULL;
float** g_cov_mat = NULL; 
double* g_pheno = NULL;
double* g_xty = NULL;
double* g_cov_R = NULL;
double g_yy = 0;
vector<string> g_snp_names;
vector<double> g_snp_var;
vector<double> g_snp_mean;
vector<int> g_snp_map;
int g_precomp = 5000;
double* g_xtx = NULL;
vector<int> g_true_gamma;
vector<int> g_init_model;
vector<double> g_true_beta;
int g_true = 0;

/******** IO *********/
string g_out_prefix = "bvsr_test";
ofstream LOG;
ofstream PATH;
ofstream MOD;
ofstream YHAT;
ofstream OUT;
ofstream GAMMA;
ofstream ERR;
map<string, int> g_paras;

/******** time *********/
int g_output_time = 0;
struct timeval g_t_beg;
struct timeval g_t_end;
map<string, double> g_time;
map<string, int> g_time_call;

/******* model *******/
int g_use_kappa = 1; 
int g_rw = 0;
int g_skip = 0; 
int g_approx_informed = 0;
double g_upper_add = 2;
double g_lower_add = 0;
double g_upper_del = 1;
double g_lower_del = 0;
double g_min_pi = 0;
double g_max_pi = 1;
int g_lbsq = 0;
//double g_gg = 1.0;
double g_gg_exp = 3.0;
double g_gg_const = 0.0;
double g_sigma = 0.2; 
double g_kappa = 2;
int g_pi_range_set = 0;
int g_pi_prior_alpha = 0;
int g_pi_prior_beta = 1;
int g_max_model_size = 1000;
int g_min_model_size = 0;
int g_log_uniform_h = 1;

/******** ridge *******/
double g_diag_c = 0;
double g_icf_abs_tol = 1e-6;
int g_icf_min_p = 30;
map<int, int> g_icf_iter;
int g_chol_call = 0;

/******* mcmc settings *******/
int g_max_jump = 1;
double g_long_range = 0.2; 
double g_prop_add = 0.5; 
string g_last_bvsr;
int g_start_size = 1;
int g_rb_thin = 1;
double g_jump_h2 = 0.1;
int g_ordered_start = 0;
int g_h_func = 0;
double g_h_A = 0.5;

/******* mcmc ********/
int g_mcmc_i = 0;
int g_count = 0;
int g_rb_count = 0;
int g_mcmc_iter = 1000;
int g_last_iter = 0;
int g_mcmc_warm = 0;
int g_mcmc_stay = 0;
int g_continue = 0;
int g_exact_bf = 0;
vector<double> g_single_bf;
vector<double> g_single_h2;
vector<double> g_pip;
vector<double> g_beta_pos;
vector<double> g_beta_neg;
vector<double> g_rb_pip;
vector<double> g_rb_beta_pos;
vector<double> g_rb_beta_neg;
vector<int> g_single_ord;
vector<int> g_single_id;
vector<int> g_propose;
vector<int> g_accept;
map<int, double> g_model_size;
int g_pmin = -1;
int g_pmax = -1;
int g_pset = 0;
double g_log_IS = NINF; 

/****** gsl *******/
const gsl_rng_type* gsl_type;
gsl_rng* gsl_r;
int gsl_seed;

