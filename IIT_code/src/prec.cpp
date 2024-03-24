#include "prec.h"

using namespace std;

void calc_globals (void){	
	record_time();
	if (g_n_cov > 0){	
		Chol_decomp_g_cov();	
		
		double* wy = allocate1D(g_n_cov);
		double* cb = allocate1D(g_n_cov);
		mat_vec_mul(g_cov_mat, g_pheno, g_n_cov, g_n_sub, wy, 0);
		Chol_solve(g_n_cov, g_cov_R, wy, cb);
		double* hwy = allocate1D(g_n_sub);
		mat_vec_mul(g_cov_mat, cb, g_n_cov, g_n_sub, hwy, 1);
		for (int i=0; i<g_n_sub; i++){
			g_pheno[i] -= hwy[i];
		}
		free1D(hwy);
		free1D(wy);
		free1D(cb);
	
	}
	
	g_yy = vec_vec_mul(g_pheno, g_pheno, g_n_sub);
	g_xty = allocate1D(g_n_snp);
	mat_vec_mul(g_data_mat, g_pheno, g_n_snp, g_n_sub, g_xty, 0);		

	record_time("CalcGlobal");
	return;
}

void calc_single_bf(vector< pair<int,double> >& bfs){
	record_time();		
	for (int i=0; i<g_n_snp; i++){
		double xy = g_xty[i];
		double xx = g_n_sub * g_snp_var.at(i);
		double yhy = xy * xy / xx;	
		double log_bf = -0.5 * g_n_sub * log(1.0 - yhy/g_yy + g_gg_const) - (0.5 * g_gg_exp + g_kappa) * log(g_n_snp);
		g_single_bf.push_back(log_bf);
		g_single_h2.push_back(yhy);
		pair<int, double> sing = make_pair(i,  log_bf);
		bfs.push_back(sing);
	}
	sort(bfs.begin(), bfs.end(), compare_pair);
	record_time("SingleBF");
	return;
}

void sort_single_bf(void){
	record_time();
	vector< pair<int, double> > single_prob; 
	for (int i=0; i<g_n_snp; i++){
		pair<int, double> sing = make_pair(i, g_single_bf[i]);
		single_prob.push_back(sing);
	}	
	sort(single_prob.begin(), single_prob.end(), compare_pair);
	g_single_id.clear();
	g_single_ord.clear();
	g_single_ord.assign(g_n_snp, 0);
	for (int i=0; i<g_n_snp; i++){
		g_single_id.push_back(single_prob[i].first);	
		g_single_ord[single_prob[i].first] = i;
	}
	record_time("SingleBF");
	return;
}

// no longer performs filtering. only calculates single snp bf. 
void filter_snp(void){
	vector< pair<int, double> > bfs;
	calc_single_bf(bfs);
	g_snp_map.assign(g_n_snp, -1);	
	for (int i=0; i<g_n_snp; i++){ g_snp_map[i] = i; }	
	sort_single_bf();
	return;
}

void calc_xtx (void){
	record_time();
	cerr << "Precalculation of part of the gram matrix" << endl;
	if (g_precomp > g_n_snp){g_precomp = g_n_snp;}
	int p = g_precomp + g_n_cov;
	g_xtx = allocate1D( p * p , 0 );
	
	// calculate X_c^t X_c. 
	for (int i = 0; i < g_n_cov; i++){
		int r = i * p;
		for (int j = 0; j <= i; j++){
			g_xtx[r + j] = vec_vec_mul(&g_cov[i * g_n_sub], &g_cov[j * g_n_sub], g_n_sub);
		}
	}

 	// if g_n_cov > 0, calculate X_c^t X_a (X_a = X_precomp)
	if (g_n_cov > 0){
		for (int i = g_n_cov; i < p; i++){
			int r = i * p;
			for (int j = 0; j < g_n_cov; j++){
				g_xtx[r + j] = vec_vec_mul(&g_data[g_single_id[i - g_n_cov] * g_n_sub], &g_cov[j * g_n_sub], g_n_sub);
			}
		}
	}

	// calculate X_a^t X_a. 
	for (int i = g_n_cov; i < p; i++){
		int r = i * p;
		for (int j = g_n_cov; j < i; j++){
			g_xtx[r + j] = vec_vec_mul(&g_data[g_single_id[i - g_n_cov] * g_n_sub], &g_data[g_single_id[j - g_n_cov] * g_n_sub], g_n_sub);
		}
	}

	// make g_xtx symmetric. 
	for (int i=0; i<p; i++){
		for (int j=i+1; j<p; j++){
			g_xtx[i * p + j] = g_xtx[j * p + i];
		}
	}	
	record_time("PrecalcXX");
	return;
}

void precalc (void){
	g_gg_const = pow(g_n_snp, -g_gg_exp);
	if (g_gg_const > 1e-6){
		g_gg_const = 1.0 / ( pow(g_n_snp, g_gg_exp) - 1.0 ); 
	}
	if (g_gg_const < 0){cerr << "Parameter g is wrong" << endl; g_gg_const = 0.0;}	
	calc_globals();
	filter_snp();
	calc_xtx();
}


