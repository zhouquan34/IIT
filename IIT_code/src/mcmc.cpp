#include "mcmc.h"

const int g_n_prints = 20;
using namespace std;

void mcmc_sample(model*& last_model){
	model* proposed_model = new model();
	int k = last_model->propose(); 
	if (k < 0){
		proposed_model->del_snp(last_model, -1-k);
	}else{
		proposed_model->add_snp(last_model, k);	
	}	

	if (g_mcmc_i > 0 ){
		string paras = last_model->para_str();
		string model = last_model->model_str();
		PATH << g_mcmc_i - 1 <<  "\t"  << paras << endl;
		MOD << model << endl;
	}
	if (g_mcmc_i > g_mcmc_warm + 1) {
		last_model -> update_marginal();
	}
	delete(last_model);
	last_model = proposed_model;		
	return;
}


void mcmc_init(model*& model){
	if (g_max_model_size > g_n_snp) {g_max_model_size = g_n_snp;}
	if (g_min_model_size < 0){g_min_model_size = 0;}
	if ((g_upper_add - g_lower_add < 1e-3) && (g_upper_del - g_lower_del < 1e-3)){
		g_rw = 1; 
	}
	g_rb_thin = 1;
	if (g_rw == 1){
		g_rb_thin = g_mcmc_iter + g_mcmc_warm + g_last_iter + 100; // no RB
		cerr << "Use random walk proposals and no Rao-Blackwellization" << endl;	
	}

	if (g_skip == 1){
		g_rb_thin = g_mcmc_iter + g_mcmc_warm + g_last_iter + 100;
		cerr << "No Rao-Blackwellization" << endl;	
	}
	
	g_pip.clear(); g_pip.assign(g_n_snp, NINF);
	g_rb_pip.clear(); g_rb_pip.assign(g_n_snp, NINF);
	g_beta_pos.clear(); g_beta_pos.assign(g_n_snp, NINF);
	g_rb_beta_pos.clear(); g_rb_beta_pos.assign(g_n_snp, NINF);
	g_beta_neg.clear(); g_beta_neg.assign(g_n_snp, NINF);
	g_rb_beta_neg.clear(); g_rb_beta_neg.assign(g_n_snp, NINF);
	
	g_propose.assign(3,0);
	g_accept.assign(3,0);

	if (g_init_model.size() > 0){  
		model->initialize(g_init_model);
		LOG << "Starting model size = " << g_init_model.size() << endl;
	}else{
		vector<int> s;	
		if (g_start_size >= 0){
			vector<int> ids;
			for (int i=0; i<g_n_snp; i++){ ids.push_back(i); }
			if (g_ordered_start == 0){
				random_shuffle(ids.begin(), ids.end());
			}
			for (int i=0; i<g_start_size; i++){
			//	s.push_back(g_single_id[i]);
				s.push_back(ids[i]);
			}
		}else{
			for (int i=0; i<g_n_snp; i++){
				if (g_true_gamma[i] == 1){
					s.push_back(i);
				}
			}	
		}
		model->initialize(s);
		LOG << "Starting model size = " << s.size() << endl;
	}
	
	g_chol_call=0;

	return;
}

int calc_refresh(void){
	int u = 0;
	if (g_continue == 0){ u = (g_mcmc_warm + g_mcmc_iter) / 50; }
	else{ u = g_mcmc_iter / 50; }
	if (u <= 0){u = 1;}
	return u;
}

void calc_start_end(int& s, int& e){
	if (g_continue == 0){
		s = 1;
		e = g_mcmc_iter + g_mcmc_warm;
	}else{
		s = g_last_iter + g_mcmc_warm + 1;
		e = g_mcmc_iter + g_mcmc_warm + g_last_iter;
	}
}

void mcmc (void){
	cerr << "Initializing MCMC" << endl;
	g_diag_c = 0.0; 
	g_proposal_weights.assign(g_n_snp, 0);
	model* my_model = new model(); 
	mcmc_init(my_model);
	
	int print_unit = calc_refresh();
	g_mcmc_stay = 0;
	cerr << "MCMC sampling" << endl;
	if (g_true == 0){
		PATH << "No.\tModelSize\ttau\tlog-like\tpve\tIS\tMaxNB" << endl;
	}else{
		// se means squared error. 
		PATH << "No.\tModelSize\ttau\tlog-like\tpve\tIS\tMaxNB\tFP\tFN\tgamma-se\tbeta-se" << endl;
	}
	MOD << "#Iteration index matched with the PATH file" << endl;
	/*for (int i=0; i<g_snp_names.size(); i++){
		string snp = g_snp_names[i];
		OUT << snp << "\t";
		GAMMA << snp << "\t";
	}
	OUT << endl;
	GAMMA << endl;
	*/

	int start_iter = 0, total_iter = 0;
	calc_start_end(start_iter, total_iter);
	int print_int = (total_iter - start_iter + 1)/g_n_prints;
	for (g_mcmc_i = start_iter; g_mcmc_i <= total_iter; g_mcmc_i ++){
		if ( g_mcmc_i % print_unit == 0){cerr << '='; fflush(stdout);}
		mcmc_sample(my_model);
		if (g_mcmc_i > g_mcmc_warm + 1){
			if ( (g_mcmc_i - g_mcmc_warm - 1) % print_int == 0){ mcmc_mean();}
		}
	}
	
	string paras = my_model->para_str();
	string model = my_model->model_str();
	my_model -> update_marginal();
	mcmc_mean();	

	delete(my_model);	
	mcmc_output();	
	
	return;
}	

