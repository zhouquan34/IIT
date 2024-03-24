// let the column of g_data g_cov = g_n_sub
// both X, Y, Cov are centered

#include "readin.h"

using namespace std;
int g_MISS = -2e5; 

void set_g_data_pointer (void){
	g_data_mat = allocate1D_float_pointer(g_n_snp);	
	for (int i=0; i<g_n_snp; i++){
		g_data_mat[i] = &g_data[i * g_n_sub];
	}
	return;
}

void set_g_n_sub (int n, string file){
	if (g_n_sub == 0){ 
		g_n_sub = n;
	}else{
		if (g_n_sub != n){
			cerr << "Wrong number of subjects in file " << file << endl;
			exit(EXIT_FAILURE);
		}
	}
	return;
}

int set_g_n_snp (int n){
	if (g_n_snp == 0){ 
		g_n_snp = n;
		g_data = allocate1D_float(g_n_sub * g_n_snp); 
		return 0;
	}else{
		int k = g_n_snp;
		float* new_data = allocate1D_float(g_n_sub * (g_n_snp + n));
		copy1D(new_data, g_data, g_n_sub * g_n_snp); // Each row is a variable!
		g_n_snp += n;
		free1D(g_data);
		g_data = new_data;
		return k;
	}
}

int read_one_line (stringstream& ss, int c){
	vector<int> snp;
	int na = 0; double sum = 0.0;
	string v;
	while (ss >> v){
		if (v.find("?")!=string::npos || v.find("N")!=string::npos){
			snp.push_back(g_MISS - 1);
			g_missing ++;
		}else{
			double g = atof(v.c_str());
			snp.push_back(g);
			na ++;
			sum += g;
		}	
	}

// isinf method fails! not clear why.
	double mean = sum / na;
	double var = 0.0;
	for (int i=0; i<g_n_sub; i++){
		if (snp[i] <= g_MISS){
//		if (isinf(snp[i])){
			g_data[c] = 0.0;	
		}else{
			double g =  snp[i] - mean;
			g_data[c] = g;
			var += g*g;
		}
		c ++;
	}	
	if (var == 0.0){ cerr << "Wrong! Please remove SNPs with no minor alleles." << endl; }
	g_snp_var.push_back(var/g_n_sub);
	g_snp_mean.push_back(mean);

	return c;
}

void error_msg_one (string file){
	cerr << "Wrong! File " << file << " should have equal number of elements every row." << endl;
	exit(EXIT_FAILURE);
}

void missing_msg(void){
	if (g_missing == 0){return;}
	double all = g_n_sub * g_n_snp;
	double f = g_missing/all;
	cerr << "Overall missing rate = " << f << endl;
}

void read_meang(string file){ // n_cov x (g_n_sub + 3)
	record_time();
	cerr << "Reading mean genotype file " << file << endl;
	int nrow = 0, ncol = 0;
	get_row_col(file, nrow, ncol);
	ncol -= 3;	
	set_g_n_sub(ncol, file);
	int c = set_g_n_snp(nrow) * g_n_sub; // starting point in g_data	

	ifstream DAT(file.c_str()); 
	string s;  
	while (getline(DAT, s))	{
		replace(s.begin(), s.end(), ',', ' ');
		stringstream ss(s); 
		string name, item;
		ss >> name >> item >> item;	
		g_snp_names.push_back(name);
		c = read_one_line(ss, c);
	}
	DAT.close();	

	if (c != g_n_snp * g_n_sub){error_msg_one(file);}
	missing_msg();
	record_time("ReadData");
	
	set_g_data_pointer();
	return;
}


void read_plink(string file){ // (n_cov + 1) x (g_n_sub + 6)
	record_time();
	cerr << "Reading PLINK A-transpose file " << file << endl;
	int nrow = 0, ncol = 0;
	get_row_col(file, nrow, ncol);
	nrow --; ncol -= 6;	
	set_g_n_sub(ncol, file);
	int c = set_g_n_snp(nrow) * g_n_sub;

	ifstream DAT(file.c_str()); 
	string s; int r = 0; 
	while (getline(DAT, s))	{
		replace(s.begin(), s.end(), ',', ' ');
		stringstream ss(s); 
		if (r == 0){
			r ++;
		}else{
			string name, item;
			ss >> item >> name >> item >> item >> item >> item;	
			g_snp_names.push_back(name);
			c = read_one_line(ss, c);
		}
	}
	DAT.close();	

	if (c != g_n_snp * g_n_sub){error_msg_one(file);}
	missing_msg();
	record_time("ReadData");
	
	set_g_data_pointer();
	return;
}


void read_matrix(string file){ // (g_n_sub + 1) x n_cov
	record_time();
	cerr << "Reading matrix covariate file " << file << endl;
	int nrow = 0, ncol = 0;
	get_row_col(file, nrow, ncol);
	nrow --;	
	set_g_n_sub(nrow, file);
	int c = set_g_n_snp(ncol);
	int start = c * g_n_sub;

	ifstream DAT(file.c_str()); 
	string s; int row = 0;
	while (getline(DAT, s))	{
		replace(s.begin(), s.end(), ',', ' ');
		stringstream ss(s); 
		if (row == 0){
			string name;
			while (ss >> name){
				g_snp_names.push_back(name);
			}
		}else{
			float p;
			int col = 0;
			while (ss >> p){
				g_data[g_n_sub * (c + col) + row - 1] = p;
				col ++;
			}
			if (col != ncol){error_msg_one(file);}
		}
		row ++;
	}
	DAT.close();	

	center_rows(&g_data[start], ncol, g_n_sub);
	record_time("ReadData");
	
	set_g_data_pointer();
	return;
}

void read_pheno(string file){
	record_time();
	cerr << "Reading phenotype file " << file << endl;
	vector<double> tmp;
	ifstream DAT(file.c_str() ); 
	string s; double sum = 0.0;
	while (getline(DAT, s))	{
		stringstream ss(s); 
		double p;
		if (ss >> p){
			tmp.push_back(p);
			sum += p;
		}
	}
	DAT.close();	
	int n = tmp.size();
	double mean = sum / n;
	set_g_n_sub(n, file);
	
	// the phenotype is centered
	g_pheno = allocate1D(n);
	for (int i=0; i<n; i++){
		g_pheno[i] = tmp[i] - mean;
	}
	record_time("ReadData");
	return;
}

void read_true(string file){
	g_true_beta.assign(g_n_snp, 0);
	g_true_gamma.assign(g_n_snp, 0);
	if (file.empty()){return;}
	
	record_time();
	cerr << "Reading true beta file " << file << endl;
	ifstream DAT(file.c_str() ); 
	string s; 
	while (getline(DAT, s))	{
		stringstream ss(s); 
		int k; double b; 
		if (ss >> k >> b){
			g_true_gamma[k-1] = 1;
			g_true_beta[k-1] = b;
		}
	}
	DAT.close();	
	record_time("ReadData");
	return;
}


void read_omit(string file){
	cerr << "Reading skipping file " << file << endl;
	ifstream DAT(file.c_str() ); 
	string s; 
	while (getline(DAT, s))	{
		stringstream ss(s);  int k;  
		if (ss >> k){ g_pass.push_back(k);}
	}
	DAT.close();	
	
	if (g_pass.size() != g_n_snp){
		cerr << "The file " << file << " has the wrong number of SNPs!" << endl; 
	}	

	return;
}


