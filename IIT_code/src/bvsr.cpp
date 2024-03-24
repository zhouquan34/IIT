#include "bvsr.h"

using namespace std;

void exit_bvsr(void){
	gsl_exit();
	free1D(g_data);
	free1D(g_cov);
	free1D(g_pheno);
	free1D(g_xty);
	free1D(g_cov_R);
	free1D(g_xtx);
	free1D_pointer(g_data_mat);
	free1D_pointer(g_cov_mat);
	return;
}

void calc_true(void){	
	model* m = new model(); 
	vector<int> s;
	for (int i=0; i<g_n_snp; i++){
		if (g_true_gamma[i] == 1){s.push_back(i);}
	}
	m->initialize(s);
	cerr << "True model log-likelihood: " << m->model_ll() << endl;	
	LOG << "True ll = " << m->model_ll() << "\n" << endl;	
}

int main(int argc, char** argv){
	cout.precision(15);

	read_args(argc, argv);
	gsl_init();
	open_outputs(g_out_prefix);
	write_log_cmds(argc, argv);
	precalc();
	
	if (g_true == 1){calc_true();}	

	clock_t begin_t = clock();
	//cout << "a: " << g_h_A << "; h " << g_h_func << endl; 
	mcmc();
	clock_t end_t = clock();
	double used_t =( (double) end_t - begin_t)/CLOCKS_PER_SEC;
	write_time(used_t);
	
	if (g_output_time == 1){ output_time();}
	
	close_outputs();
	exit_bvsr();
	return 0;
}


