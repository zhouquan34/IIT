#include "lalg.h"

using namespace std;

double* allocate1D (int dim){
	double* m = (double*) malloc((size_t) (dim * sizeof(double)));
	return m;
}

float* allocate1D_float (int dim){
	float* m = (float*) malloc((size_t) (dim * sizeof(float)));
	return m;
}

float** allocate1D_float_pointer (int dim){
	float** m = (float**) malloc((size_t) (dim * sizeof(float*)));
	return m;
}

double* allocate1D (int dim, int v){
	double* m = (double*) malloc((size_t) (dim * sizeof(double)));
	fill_n(m, dim, (double) v);
	return m;
}

double** allocate2D (int dim){
	double** m;
	m = (double**) malloc((size_t) (2 * sizeof(double*))); // every row is double*
	m[0] = (double*) malloc((size_t) (2 * dim * sizeof(double)));
	fill_n(m[0], 2 * dim * sizeof(double), 0);
	m[1] = m[0] + dim;
	return m;
}

template <typename T> void free1D(T* m){
	if (m == NULL){return;}
	free(m);
	m = NULL;
	return;
}
template void free1D <double>(double* m); 
template void free1D <float>(float* m); 

void free1D_pointer(float** m){
	if (m == NULL){return;}
	free(m);
	m = NULL;
	return;
}

void free2D (double** m){
	if (m == NULL) {return;}
	free(m[0]);
	free(m);
	m = NULL;
	return;
}

template <typename T, typename U> void copy1D(T* m1, const U* m2, int p){
	for (int i=0; i<p; i++){
		m1[i] = m2[i];
	}
	return;
}
template void copy1D <double, double>(double* m1, const double* m2, int p); 
template void copy1D <double, float>(double* m1, const float* m2, int p); 
template void copy1D <float, float>(float* m1, const float* m2, int p); 
template void copy1D <float, double>(float* m1, const double* m2, int p); 

// const float** m2 means values are constant. pointer to pointer to constant float. 
// cannot assign pointer to constant to pointer to non_constant. 
void copy1D (float** m1, float** m2, int p){
	for (int i=0; i<p; i++){
		m1[i] = m2[i];
	}
	return;
}

// M (p x n) V (n x 1) = Z
void mat_times_vec (double* m, double* v, int p, int n, double* z){	
	int i = 0, j = 0, r = 0;
	int np = n*p;
	double s = 0.0;
	for (i = 0; i < np; i++){
		s += m[i] * v[j];
		j ++; 
		if ( j == n){
			z[r] = s;
			r ++;
			s = 0.0;
			j = 0;
		} 
	}
	return;
}


double L2_norm (double* v, int p){
	double s = 0.0;
	for (int i=0; i<p; i++){
		s += v[i] * v[i];
	}
	return sqrt(s);
}


template <typename T, typename U>
double vec_vec_mul (T* v1, U* v2, int p){
	T* vp1 = &v1[0]; U* vp2 = &v2[0];
	double s = 0.0;
	for (int i=0; i<p; i++){
		s += (*vp1) * (*vp2);
		vp1 ++;
		vp2 ++;
	}
	return s;
}

template double vec_vec_mul<double, double> (double* v1, double* v2, int p);
template double vec_vec_mul<double, float> (double* v1, float* v2, int p);
template double vec_vec_mul<float, double> (float* v1, double* v2, int p);
template double vec_vec_mul<float, float> (float* v1, float* v2, int p);

template <typename T>
void center_row(T* m, int n, int calc_var){
	double sum = 0.0;
	for (int i=0; i<n; i++){
		sum += m[i];
	}
	double mean = sum / n;
	if (calc_var == 1){
		double var = 0.0;
		for (int i=0; i<n; i++){
			m[i] = m[i] - mean;
			var += m[i] * m[i];
		}
		var = var/n;
		g_snp_mean.push_back(mean);
		g_snp_var.push_back(var);
	}else{
		for (int i=0; i<n; i++){
			m[i] = m[i] - mean;
		}
	}
	return;
}

template <typename T>
void center_rows(T* m, int p, int n){
	record_time();
	for (int i=0; i<p; i++){
		center_row(&m[i * n], n, 1);
	}
	record_time("Centering");
	return;
}	
template void center_rows <double> (double* m, int p, int n);
template void center_rows <float> (float* m, int p, int n);


double center_array(double* v, int p){
	double ss1 = 0.0;
	double ss2 = 0.0;
	for (int i=0; i<p; i++){
		ss1 += v[i];
		ss2 += v[i] * v[i];
	}
	double sst =  ss2 - ss1 * ss1 / (double) p;
	double mean = ss1 / p;
	for (int i=0; i<p; i++){
		v[i] -= mean;	
	}
	return sst / p;
}


void array_del (const vector<int>& del, int p_last, double* v_last, double* v){
	int s = del.size();
	int c = 0;
	int d = 0;
	for (int i=0; i< s; i++){
		int r = del.at(i);
		if (r - c > 0){
			copy1D(&v[d], &v_last[c], r - c);
		}
		d = d + r - c;
		c = r + 1;		
	}
	if (c < p_last){	
		copy1D(&v[d], &v_last[c], p_last - c);
	}
	return;
}

template <typename T> 
void matrix_del (const vector<int>& del, int p_last, int n_last, T* m_last, T* m){
	int s = del.size();
	int c = 0;
	int d = 0;
	for (int i=0; i< s; i++){
		int r = del.at(i);
		if (r - c > 0){
			copy1D(&m[d * n_last], &m_last[c * n_last], (r-c) * n_last);
		}
		d = d + r - c;
		c = r + 1;		
	}
	if (c < p_last){	
		copy1D(&m[d * n_last], &m_last[c * n_last], (p_last - c) * n_last);
	}
	return;
}
template void matrix_del <float> (const vector<int>& del, int p_last, int n_last, float* m_last, float* m);
template void matrix_del <double> (const vector<int>& del, int p_last, int n_last, double* m_last, double* m);
 
double array_sst(double* v, int p){
	double ss1 = 0.0;
	double ss2 = 0.0;
	for (int i=0; i<p; i++){
		ss1 += v[i];
		ss2 += v[i] * v[i];
	}
	return ss2 - ss1 * ss1 / (double) p;
}

template <typename T, typename U> 
void mat_times_vec_opt (T** m, U* v, int p, int n, double* z){
	double* zpos = &z[0];
	int k = p/4;
	for (int i=0; i<k; i++){
		T* mp1 = m[i*4];
		T* mp2 = m[i*4 + 1];
		T* mp3 = m[i*4 + 2];
		T* mp4 = m[i*4 + 3];	
		double z1 = 0, z2 = 0, z3 = 0, z4 = 0;
		U* vpos = &v[0];
		for (int j=0; j<n; j++){
			z1 += (*mp1) * (*vpos);
			z2 += (*mp2) * (*vpos);
			z3 += (*mp3) * (*vpos);
			z4 += (*mp4) * (*vpos);
			vpos ++;
			mp1 ++;
			mp2 ++;
			mp3 ++;
			mp4 ++;
		}
		*zpos = z1;  zpos ++;
		*zpos = z2;  zpos ++;
		*zpos = z3;  zpos ++;
		*zpos = z4;  zpos ++;
	}
	
	for (int i=k*4; i<p; i++){
		double z0 = 0;
		U* vpos = &v[0];
		T* mp1 = m[i];
		for (int j=0; j<n; j++){
			z0 += (*mp1) * (*vpos);
			vpos ++;
			mp1 ++;
		}
		*zpos = z0; zpos ++;
	}

}

template <typename T, typename U>
void mat_t_times_vec_opt (T** m, U* v, int p, int n, double* z){
	U* vpos = &v[0];
	int k = p/4;
	for (int i=0; i<k; i++){
		double* zpos = &z[0];
		U v1 = *vpos; vpos ++; 
		U v2 = *vpos; vpos ++; 
		U v3 = *vpos; vpos ++; 
		U v4 = *vpos; vpos ++; 
		T* mp1 = m[i*4];
		T* mp2 = m[i*4 + 1];
		T* mp3 = m[i*4 + 2];
		T* mp4 = m[i*4 + 3];		
		for (int j=0; j<n; j++){
			*zpos += (*mp1) * v1;
			*zpos += (*mp2) * v2;
			*zpos += (*mp3) * v3;
			*zpos += (*mp4) * v4;
			zpos ++;
			mp1 ++;
			mp2 ++;
			mp3 ++;
			mp4 ++;
		}
	}
	
	for (int i=k*4; i<p; i++){
		double* zpos = &z[0];
		U v0 = *vpos; vpos ++;
		T* mp1 = m[i];
		for (int j=0; j<n; j++){
			*zpos += (*mp1) * v0;
			zpos ++;
			mp1 ++;
		}
	}

	return;
}

// t(M (p x n)) V (p x 1) = Z
template <typename T, typename U>
void mat_vec_mul (T** m, U* v, int p, int n, double* z, int trans){
	if (trans == 0){
		mat_times_vec_opt (m, v, p, n, z);
	}else{
		fill_n(z, n, 0);
		mat_t_times_vec_opt (m, v, p, n, z);
	}
	return;
}
template void mat_vec_mul<double, double> (double** m, double* v, int p, int n, double* z, int trans);
template void mat_vec_mul<double, float> (double** m, float* v, int p, int n, double* z, int trans);
template void mat_vec_mul<float, double> (float** m, double* v, int p, int n, double* z, int trans);
template void mat_vec_mul<float, float> (float** m, float* v, int p, int n, double* z, int trans);


