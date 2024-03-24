#include "generic.h"

using namespace std;

bool compare_pair(pair<int, double> x, pair<int, double> y){ return x.second>y.second; }

double log_sum_log (double a, double b){
	if (a < b){
		double c = a; a = b; b = c;
	}
	return a + log(1.0 + exp(b-a));
}

void gsl_init (void){
	if (g_pset == 1){
		if (g_pmin > -1){
			g_min_pi = (double) g_pmin / (double) g_n_snp;	
			g_pi_range_set = 1;
		}
		if (g_pmax > -1){
			g_max_pi = (double) g_pmax / (double) g_n_snp;	
			g_pi_range_set = 1;
		}
		cerr << "Setting pi-max = " << g_max_pi << " and pi-min = " << g_min_pi << endl;
	}
	
	gsl_rng_env_setup();
	gsl_type = gsl_rng_default;
	gsl_r = gsl_rng_alloc(gsl_type);
	gsl_rng_set(gsl_r, gsl_seed);
	return;
}

void gsl_exit (void){
	gsl_rng_free(gsl_r);
	return;
}

void get_row_col (const string file, int& nrow, int& ncol ){
	ifstream DAT(file.c_str() ); 
	string s;
	int row = 0;
	while (getline(DAT, s))	{
		replace(s.begin(), s.end(), ',', ' ');
		row ++ ;
		if (s.empty()){row --;}
		if (row == 2){
			stringstream ss(s);
			int c = 0; string tmp;
			while (ss >> tmp){
				c++;
			}
			ncol = c;
		}
	}
	nrow = row;
	DAT.close();
	return;
}

template <typename T>
void print_array(const T* v, int p){
	for (int i = 0; i<p; i++){
		cout << v[i] << " ";
	}
	cout << endl;
	return; 
}
template void print_array <double> (const double* v, int p);
template void print_array <float> (const float* v, int p);

template <typename T>
void print_vector(vector<T> v){
	for (int i = 0; i<v.size(); i++){
		cout << v[i] << ", ";
	}
	//cout << endl;
	return; 
}
template void print_vector<int> (vector<int> v);
template void print_vector<float> (vector<float>);
template void print_vector<double> (vector<double>);

template <typename T>
void print_matrix(const T* v, int p, int k){
	for (int i = 0; i < p; i += k){
		for (int j=0; j<k; j++){
			cout << v[i + j] << " ";
		}
		cout << endl;
	}
	return;
}
template void print_matrix <double> (const double* v, int p, int k);
template void print_matrix <float> (const float* v, int p, int k);

void write_log_cmds(int argc, char** argv){
	LOG << "# Command: "; 
	for (int i=0; i<argc; i++){ LOG << argv[i] << " ";} 
	LOG  << endl;
	if (g_continue == 1){
		cerr << "# Continue MCMC sampling of dataset " << g_last_bvsr << endl;
		cerr << "# Please note that only '-b', '-s', '-o', '-r' are enabled." << endl;
		cerr << "# " << g_out_prefix << ".beta.txt will be computed using all runs; all the other information is only for this run." << endl;	
		LOG << "\n# This is the continuation of run " << g_last_bvsr << endl;
		LOG << "# " << g_out_prefix << ".beta.txt will be computed using all runs; all the other information is only for this run." << endl;	
		LOG << "# Try 'fastBVSR-post n_burn_in path_file_1 path_file_2 ...' to recompute the posterior for multiple runs." << endl;
	}
	LOG << endl;
	LOG << "N_subject = " << g_n_sub  << "; N_predictor = " << g_n_snp << endl;
}


int weighted_sample (vector<double>& p){
	double s = 0;
	for (int i=0; i<p.size(); i++){ s += p[i]; }

	double r = gsl_rng_uniform(gsl_r) ;
	double cum = 0;
	for (int i = 0; i<p.size(); i++){
		cum += p[i]/s;
		if (r < cum){ return i;}
	}
	cerr << "Wrong in sampling: " << r << ", " << cum  << endl;
	return 0;
}

int log_weighted_sample (vector<double>& log_p, double& IS){
	vector<double> cum; 
	double s = NINF; 
	double q = 0;
	for (int i=0; i<log_p.size(); i++){ 
		s = log_sum_log(s, log_p[i]); 
		cum.push_back(s); 
	}
	IS = s;	

	double log_r = log(gsl_rng_uniform(gsl_r));
	for (int i = 0; i<log_p.size(); i++){
		if ( log_r < cum[i] - s ){ 
			q = log_p[i] - s; 
			return i;
		}
	}
	cerr << "Wrong in sampling: " << log_r << ", " << s  << endl;
	return 0;
}

void write_time(double used_t){
	LOG << "Time used by MCMC = " ;
	if (used_t > 2000){LOG << used_t/60.0 << " m" << endl;}
	else{LOG << used_t << " s" << endl;}
	return;
}

void record_time (string x){
	if (x == ""){ 
		gettimeofday(&g_t_beg, NULL);
		return;
	}else{
		gettimeofday(&g_t_end, NULL);
		double dt = (double) (g_t_end.tv_usec - g_t_beg.tv_usec) / 1000.0 + (g_t_end.tv_sec - g_t_beg.tv_sec) * 1000.0;	
		if (g_time.find(x) != g_time.end()){
			g_time[x] += dt / 1000.0;
			g_time_call[x] ++;
		}else{
			g_time[x] = dt / 1000.0;	
			g_time_call[x] = 1;
		}
		return;
	}
}


void output_time(void){
	LOG << "\nTime usage details (in seconds):" << endl;
	map<double, string> tmp;
	map<string, double>::iterator it;
	for (it=g_time.begin(); it!=g_time.end(); it++){
		tmp[it->second] = it->first;
	}
	map<double, string>::iterator it2;
	for (it2=tmp.begin(); it2!=tmp.end(); it2++){
		LOG << it2->second << "\t" << it2->first << "\t" << g_time_call[it2->second] << endl;
	}
	return;
}


