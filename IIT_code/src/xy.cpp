// No memory allocation in this file
#include "xy.h"

using namespace std;

// each row is a variable
void X_init(const vector<int>& s, float** X){
	int col = s.size() + g_n_cov;
	for (int i=0; i<g_n_cov; i++){
		X[i] = &g_cov[i * g_n_sub];
	}
	for (int i=g_n_cov; i<col; i++){
		X[i] = &g_data[ s[i - g_n_cov] * g_n_sub];
	}
	return;
}

void X_add (const vector<int>& add, int p_last, float** X_last, float** X){
	int s = add.size();
	copy1D(X, X_last, p_last);
	for (int i=0; i < s; i++){
		X[i + p_last] = &g_data[add[i] * g_n_sub];
	}
	return;
}


// del is sorted. del + g_n_cov is the true indices
void X_del (const vector<int>& del, int p_last, float** X_last, float** X){
//	record_time();
	int s = del.size();
	int c = 0;
	int d = 0;
	for (int i=0; i< s; i++){
		int r = del.at(i) + g_n_cov;
		if (r - c > 0){
			copy1D(&X[d], &X_last[c], r-c);
		}
		d = d + r - c;
		c = r + 1;		
	}
	if (c < p_last){	
		copy1D(&X[d], &X_last[c], p_last - c);
	}
//	record_time("UpdateX");
	return;
}


void Xy_init (float** X, int p, double* Xy){
//	record_time();
	mat_vec_mul(X, g_pheno, p, g_n_sub, Xy, 0);	
//	record_time("UpdateXy");
	return;
}

void Xy_add (const vector<int>& add, int p_last, double* Xy_last, double* Xy){
//	record_time();
	copy1D(Xy, Xy_last, p_last);
	int s = add.size();
	for (int i=0; i<s; i++){
		Xy[p_last + i] = g_xty[add[i]];
	}
//	record_time("UpdateXy");
	return;
}

void Xy_del (const vector<int>& del, int p_last, double* Xy_last, double* Xy){
//	record_time();
	int s = del.size();
	int c = 0;
	int d = 0;
	for (int i=0; i< s; i++){
		int r = del.at(i) + g_n_cov;
		if (r - c > 0){
			copy1D(&Xy[d], &Xy_last[c], r - c);
		}
		d = d + r - c;
		c = r + 1;		
	}
	if (c < p_last){	
		copy1D(&Xy[d], &Xy_last[c], p_last - c);
	}
//	record_time("UpdateXy");
	return;
}

// length of XX is 1 + 2 + ... + p = p(p+1)/2  (lower triangle)
void XX_init (float* X, int p, double* XX){
	record_time();
	int c = 0;
	for (int i=0; i<p; i++){
		for (int j=0; j<i; j++){
			XX[c] = vec_vec_mul(&X[i * g_n_sub], &X[j * g_n_sub], g_n_sub);
			c ++;
		}
		XX[c] = vec_vec_mul(&X[i * g_n_sub ], &X[i * g_n_sub], g_n_sub);
		if (i >= g_n_cov){
			XX[c] += g_diag_c;
		}
		c++;
	}
	record_time("UpdateXX");
	return;
}


double XX_get (int i, int j){
	if (i == j){
		return g_snp_var[i] * g_n_sub + g_diag_c;		
	}
	int oi = g_single_ord[i];
	int oj = g_single_ord[j];
	if (oi < g_precomp && oj < g_precomp){
		return g_xtx[ (g_precomp + g_n_cov)*(oi + g_n_cov) + oj + g_n_cov];
	}else{
		return vec_vec_mul(&g_data[i * g_n_sub], &g_data[j * g_n_sub], g_n_sub);
	}
}

// assuming no g_cov; otherwise use the commented block below. 
void XX_init (const vector<int>& subset, float** X, int p, double* XX){
	record_time();
	int c = 0;
	for (int i=0; i<p; i++){
		for (int j=0; j<=i; j++){
			XX[c] = XX_get(subset[i], subset[j]);
			c++;
		}
	}
	record_time("UpdateXX");
	return;
}


/*
void XX_init (const vector<int>& subset, float** X, int p, double* XX){
	record_time();
	int c = 0;
	int pre = g_precomp + g_n_cov; // g_xtx is pre * pre
	for (int i=0; i<g_n_cov; i++){
		for (int j=0; j<=i; j++){
			XX[c] = g_xtx[pre * i + j]; 
			c ++;
		}
	}	
	for (int i=g_n_cov; i<p; i++){
		int o = g_single_ord[subset[i - g_n_cov]];
		if (o < g_precomp){
			for (int j=0; j<g_n_cov; j++){
				XX[c] = g_xtx[pre * (o + g_n_cov) + j];
				c ++;
			}
			for (int j=g_n_cov; j<i; j++){
				if (g_single_ord[subset[j - g_n_cov]] < g_precomp){
					XX[c] = g_xtx[pre * (o + g_n_cov) + g_single_ord[subset[j - g_n_cov]] + g_n_cov];
				}else{
					XX[c] = vec_vec_mul(X[i], X[j], g_n_sub);
				}
				c ++;
			}
		}else{
			for (int j=0; j<i; j++){
				XX[c] = vec_vec_mul(X[i], X[j], g_n_sub);
				c ++;
			}
		}
		XX[c] = g_snp_var[subset[i-g_n_cov]] * g_n_sub + g_diag_c;
		c++;
	}
	record_time("UpdateXX");
	return;
}
*/


// subset is the subset of XX (not last)
void XX_add (const vector<int>& subset, const vector<int>& add, int p_last, float** X, double* XX_last, double* XX){
	record_time();
	int s = add.size();
	int c = p_last * (p_last + 1) / 2;
	copy1D(XX, XX_last, c);
	for (int i=0; i<s; i++){
		int o = g_single_ord[add[i]];
		if (o < g_precomp){
			for (int j=0; j<g_n_cov; j++){
				XX[c] = g_xtx[(g_precomp + g_n_cov) * (o + g_n_cov) + j];
				c ++;
			}
			for (int j = g_n_cov; j < p_last + i; j++){
				if (g_single_ord[subset[j - g_n_cov]] < g_precomp){
					XX[c] = g_xtx[(g_precomp + g_n_cov) * (o + g_n_cov) + g_single_ord[subset[j - g_n_cov]] + g_n_cov];
				}else{
					XX[c] = vec_vec_mul(X[p_last + i], X[j], g_n_sub);
				}
				c ++;
			}
		}else{
			mat_vec_mul(X, X[p_last + i], p_last + i, g_n_sub, &XX[c], 0);	
			c += (p_last + i);		
		}
		XX[c] = g_snp_var[add[i]] * g_n_sub + g_diag_c;
		c ++;
	}
	record_time("UpdateXX");
}


// again del + g_n_cov is the true index
void XX_del (const vector<int>& del, int p_last,  double* XX_last, double* XX){
	record_time();
	int s = del.size();
	int nxx = p_last * (p_last + 1)/2;
	set<int> to_del; // automatically sorted since it is a set  
//	cerr << "Del\t" << c << "\t" << p_last << "\t" << nxx << endl;
	for (int i=0; i<s; i++){
		int d = del[i] + g_n_cov;
		int start = d*(d+1)/2;
		for (int j=start; j<=start+d; j++){
			to_del.insert(j);	
		}
		int j = start + d;
		d ++ ;
		j += d;
		while (j < nxx){
			to_del.insert(j);
			d ++;
			j += d;
		}	
	}
	
	int p = p_last - s;	
	int c = 0;
	set<int>::iterator it = to_del.begin();
	int next_d = *it;
	for (int i=0; i<nxx; i++){
		if (i == next_d){
			it ++;
			if (it != to_del.end()){
				next_d = *it;
			}else{
				next_d = -1;
			}
		}else{
			XX[c] = XX_last[i];
			c ++;
		}
	}
	if ( c != p * (p+1)/2 ) {
		cerr << "Wrong in XX_del! " << endl; 
	}
	
	record_time("UpdateXX");
	return;
}



