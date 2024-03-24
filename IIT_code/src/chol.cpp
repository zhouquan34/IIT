// R is upper triangle

#include "chol.h"
using namespace std;

// p is the dim of next_R; k is the col id we need to figure out; diag is diag^2
void forward_sub (int k, int p, double* XX, double* next_R, double diag){
	int xx_id = k * (1 + k) / 2;

	for (int i = 0; i < k; i ++){
		if (fabs(next_R[i * p + i]) < ZERO){
			next_R[i * p + k] = 0;
		}else{
			double s = XX[xx_id];	
			for (int j=0; j<i; j++){
				s -= next_R[j * p + i] * next_R[j * p + k];
			}
			next_R[i * p + k] = s / next_R[i * p + i];
		}
		xx_id ++;
	}
	double s  = XX[xx_id] + diag;
	for (int i=0; i<k; i++){
		s -= next_R[i * p + k] * next_R[i * p + k];
	}
	if (s < ZERO){s = 0;}
	next_R[k * p + k] = sqrt(s);
	return;
}

// p is the dim of R
void Chol_add (int add, int p, double* XX, double* R, double* next_R){
	record_time();
	
	int next_p = p + add;
	int i = 0, j = 0, k = 0, row = 0;
	while (i < p * p){
		next_R[k] = R[i];
		i ++ ; j ++ ; k ++ ;
		if (j == p){
			row ++ ;
			j = row;
			k += (add + row);
			i += row;
		}
	}

	for (i = 0; i < add; i ++){
		forward_sub(p + i, p + add,  XX, next_R, 0);			
	}
	
	record_time("CholAdd");
	return;	
}

void Chol_add_no_time (int add, int p, double* XX, double* R, double* next_R, double diag){
	int next_p = p + add;
	
	int i = 0, j = 0, k = 0, row = 0;
	while (i < p * p){
		next_R[k] = R[i];
		i ++ ; j ++ ; k ++ ;
		if (j == p){
			row ++ ;
			j = row;
			k += (add + row);
			i += row;
		}
	}

	for (i = 0; i < add; i ++){
		forward_sub(p + i, p + add,  XX, next_R, diag * diag);			
	}
	return;	
}


// p is the dim of R
void Chol_del (vector<int> del, int p, double* R, double* next_R){	
	record_time();

	int dsize = del.size(); int d = 0;
	int next_p = p - dsize; 
	vector<int> index;
	
	for (int i=0; i<g_n_cov; i++){
		index.push_back(i);
	}
	for (int i=g_n_cov; i<p; i++){
		if (i == (del[d] + g_n_cov) ){
			if (d == (dsize-1)){
				for (int j = i+1; j<p; j++){index.push_back(j);}
				break;
			}else{
				d++;	
				continue;
			}
		}else{
			index.push_back(i);
		}
	}
	
	double* R1 = allocate1D( dsize * next_p);
	for (int i=0; i<next_p; i++){
		int r1 = i * next_p; 
		int r2 = i * p;
		for (int j=0; j<next_p; j++){
			next_R[r1 + j] = R[r2 + index[j]];
		}
	}
	for (int i=next_p; i<p; i++){
		int r1 = (i-next_p) * next_p; 
		int r2 = i * p;
		for (int j=0; j<next_p; j++){
			R1[r1 + j] = R[r2 + index[j]];
		}
	}
	
	// consider three cases. in next_R, in R1, or in both
	for (int i=0; i<index.size(); i++){
		for (int j = index[i]; j > i; j -- ){
			double* row_1;  double* row_2;
			if (j > next_p){
				row_1 = &R1[(j-1-next_p)*next_p];
				row_2 = &R1[(j-next_p)*next_p];
			}else if (j == next_p){
				row_1 = &next_R[(next_p-1)*next_p];
				row_2 = &R1[0];
			}else{
				row_1 = &next_R[(j-1)*next_p];
				row_2 = &next_R[j*next_p];
			}

			double a = row_1[i]; 
			double b = row_2[i];
			if (b == 0) {continue;} // should add this
			double r = sqrt(a*a + b*b);
			double c = a/r;
			double s = -b/r;
			row_1[i] = r;
			row_2[i] = 0;
			for (int k = i + 1; k < index.size(); k++){
				double a1 = row_1[k];
				double b1 = row_2[k];
				row_1[k] = c*a1 - s*b1;
				row_2[k] = s*a1 + c*b1;
			}
		}		
	}
	free1D(R1);

	record_time("CholDel");
	return;
}

// p is the dim of XX
void Chol_init (int p, double* XX, double* R){
	double* R0 = allocate1D(1); 
	R0[0] = sqrt(XX[0]);	
	Chol_add(p-1, 1, XX, R0, R);
	free1D(R0);
}


void Chol_decomp_g_cov (void){
	if (g_n_cov == 0){return;}
	g_cov_R = allocate1D(g_n_cov * g_n_cov);
	int nc = g_n_cov * (g_n_cov + 1) / 2;
	double* cc = allocate1D(nc);
	XX_init(g_cov, g_n_cov, cc);

	if (g_n_cov == 1){
		g_cov_R[0] = sqrt(cc[0]);
	}else{
		double* R0 = allocate1D(1);
		R0[0] = sqrt(cc[0]);	
		Chol_add_no_time(g_n_cov-1, 1, cc, R0, g_cov_R, (double) 0);
		free1D(R0);
	}
	free1D(cc);
	return;
}

double Chol_decomp_A (int p, double* XX, double* R, double diag){
	record_time();
	if (g_n_cov == 0){
		double* R0 = allocate1D(1); 
		R0[0] = sqrt(XX[0] + diag * diag);	
		Chol_add_no_time(p - 1, 1, XX, R0, R, diag);
		free1D(R0);
	}else{
		Chol_add_no_time(p - g_n_cov, g_n_cov, XX, g_cov_R, R, diag);	
	}
	double log_det = 0.0;
	for (int i=0; i<p; i++){
		log_det += log(R[i * p + i]);
	}

	record_time("CholDecomp");
	return 2.0 * log_det;
}

// Rt R b = z; Rt b1 = z; R b = b1
void Chol_solve (int p, double* R, double* z, double* b){
	record_time();
	
	for (int i = 0; i < p; i ++){
		double s = z[i];
		for (int j = 0; j < i; j++){
			s -= b[j] * R[j * p + i]; // transpose
		}
		b[i] = s / R[i * p + i]; 
	}

	for (int i = p - 1; i >= 0; i --){
		double s = b[i];
		int ri = i * p;
		for (int j = p - 1; j > i; j --){
			s -= b[j] * R[ri + j]; 
		}
		b[i] = s / R[i * p + i]; 
	}

	record_time("CholSolve");
	return;
}

// R b = z 
void Chol_solve_backward (int p, double* R, double* z, double* b){
	record_time();
	
	for (int i = p - 1; i >= 0; i --){
		double s = z[i];
		int ri = i * p;
		for (int j = p - 1; j > i; j --){
			s -= b[j] * R[ri + j]; 
		}
		b[i] = s / R[i * p + i]; 
	}

	record_time("CholSolve");
	return;
}



