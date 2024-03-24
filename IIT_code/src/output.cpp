#include "output.h"

using namespace std;

void open_outputs (string prefix){
	string log_file = prefix;
	log_file.append(".log.txt");
	string out_file = prefix;
	out_file.append(".beta.txt");
	string path_file = prefix;
	path_file.append(".path.txt");
	string model_file = prefix;
	model_file.append(".model.txt");	
	string gamma_file = prefix;
	gamma_file.append(".gamma.txt");
	string error_file = prefix;
	error_file.append(".error.txt");
		
	LOG.open(log_file.c_str(), ofstream::out);
	OUT.open(out_file.c_str(), ofstream::out);
//	GAMMA.open(gamma_file.c_str(), ofstream::out);
	PATH.open(path_file.c_str(), ofstream::out);
	MOD.open(model_file.c_str(), ofstream::out);
	ERR.open(error_file.c_str(), ofstream::out);
	LOG.precision(6);
	
	return;
}

void close_outputs(void){
	LOG.close();
	OUT.close();
	PATH.close();
	MOD.close();
//	GAMMA.close();
	ERR.close();
	return;
}

template <typename T> 
void post_quantile(map<T, double>& v){
	typename map<T, double>::iterator it;
	double n = 0;
	double ss = 0.0;
	for (it=v.begin(); it!=v.end(); it++){
		n += it->second;
		ss += (double) it->second * (double) it->first;
	}
	if (n > 0){
		double m = ss / n;
		vector<double> qs(5,0);
		double quant[] = {0.025, 0.05, 0.5, 0.95, 0.975};
		double k = 0.0;
		int qi = 0;
		for (it=v.begin(); it!=v.end(); it++){
			k += it->second;
			while (k/n > quant[qi] - ZERO){ 
				qs[qi] = it->first; 
				qi ++;
				if (qi == 5) {break;}
			}
			if (qi == 5) {break;}
		}
		LOG << "Mean = " << m << ";  Median = " << qs[2] << endl;
		LOG << "90% credible interval = (" << qs[1] <<  ", " << qs[3] << ")" << endl;
		LOG << "95% credible interval = (" << qs[0] <<  ", " << qs[4] << ")" << endl;	
	}
	return;
}

void mcmc_output(void){
	record_time();
	cerr << "\nOutput results" << endl;
	
	LOG << "\nModel size:" << endl;
	for (map<int, double>::iterator it=g_model_size.begin(); it != g_model_size.end(); ++it){
		it->second = exp(it->second - g_log_IS);
	//	cout << it->first << "\t" << it->second << "\n";
	}
	post_quantile(g_model_size);
	
	OUT << "SNP\tlog10(BF)\tPIP\tRB-PIP\tBeta\tRB-Beta\n";
	int total_iter = g_mcmc_iter  + g_last_iter; 
	cout << "total_iter = " << total_iter << endl;
	for (int i=0; i<g_snp_names.size(); i++){
		int smap = g_snp_map[i];
		double pip  = exp( g_pip[smap] - g_log_IS);
		double rb_pip  = exp( g_rb_pip[smap] - g_log_IS);
		double beta_p = exp( g_beta_pos[smap] - g_log_IS);
		double beta_n = exp( g_beta_neg[smap] - g_log_IS);
		double beta = beta_p - beta_n;
		double rb_beta_p = exp( g_rb_beta_pos[smap] - g_log_IS);
		double rb_beta_n = exp( g_rb_beta_neg[smap] - g_log_IS);
		double rb_beta = rb_beta_p - rb_beta_n;
		string snp = g_snp_names[i];
		OUT << snp << "\t" << g_single_bf[smap]/log(10.0) << "\t" << pip  << "\t" << rb_pip << "\t" << beta << "\t" << rb_beta << endl;
	}
	
	record_time("Output");
	cout << "Samples = " << g_count << "; RB samples = " << g_rb_count << endl;
	return;
}

void mcmc_mean(void){
	if (g_true == 0){return;}
	double e1 = 0;
	double e2 = 0; 
	double e3 = 0; 
	double e4 = 0; 
	for (int i=0; i<g_n_snp; i++){
		double pip  = exp( g_pip[i] - g_log_IS);
		double rb_pip  = exp( g_rb_pip[i] - g_log_IS);
		double beta_p = exp( g_beta_pos[i] - g_log_IS);
		double beta_n = exp( g_beta_neg[i] - g_log_IS);
		double beta = beta_p - beta_n;
		double rb_beta_p = exp( g_rb_beta_pos[i] - g_log_IS);
		double rb_beta_n = exp( g_rb_beta_neg[i] - g_log_IS);
		double rb_beta = rb_beta_p - rb_beta_n;
		e1 += (pip - g_true_gamma[i]) * (pip - g_true_gamma[i]);
		e2 += (rb_pip - g_true_gamma[i]) * (rb_pip - g_true_gamma[i]);
		e3 += (beta - g_true_beta[i]) * (beta - g_true_beta[i]);
		e4 += (rb_beta - g_true_beta[i]) * (rb_beta - g_true_beta[i]);
	}
	ERR << g_mcmc_i - 1 << "\t" << sqrt(e1/g_n_snp) << "\t" << sqrt(e2/g_n_snp) << "\t" << sqrt(e3/g_n_snp) << "\t" << sqrt(e4/g_n_snp)   << endl;
	return;
}

