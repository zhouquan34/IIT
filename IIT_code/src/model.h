#ifndef BVSR_MODEL
#define BVSR_MODEL

#define FAIL_ALPHA -1e15
#define UNDEF_LL  -1e14

#include "chol.h"
#include "prec.h"
#include "generic.h"
#include "global.h"
#include "lalg.h"

class model{
	private:
		double pi;
		double sse;
		double tau;
		double like;
		double pve_rb;
		double gamma_rb_mse;
		double beta_rb_mse;
		double IS;
		double max_nb;
		int msize;
		int mode;
		std::vector<int> subset;
		std::vector<int> selected;
		std::vector<double> neighbor_ll; 
		std::vector<double> neighbor_beta; 
		
		// allocate(subset_size + g_n_cov) allocate all pointers; this function 
		void allocate(int); 
		
		// copy everything but pointers and neighbor_ll 
		void copy(class model*);

		// calculate neighboring model ll and set it in neighbor_ll; return truncated value. 
		double add_ll (int); // add_ll(k): adding the k-th snp (in all snps). 
		double del_ll (int); // del_ll(k): deleting the k-th element in subset. 
		
		double calc_det(void);	
		void precalc_add(void);
	public:
		double* R;
		float** X;
		double* XX;
		double* Xy;
		double* Beta;
		double* Rinv; 		
		double* Rinv_xy;		
	
		model(); // assign neighbor_ll to UNDEF_LL
		~model(); // free all pointers 
		void clean(void); // free all pointers
		int model_size(void); // return number of snps
		double model_ll(void);		
	
		// find_snp(k) return the index of the subset element that equals k  
		int find_snp(int); 
		std::string para_str(void); // return the model name in string
		std::string model_str(void); // return the parameter vector in string
		
		void update_marginal(void);	
	
		// initialize using a given set of snps; allocate all and initialize X, XX, Xy, R; calculate ll. 
		void initialize(std::vector<int>&);
	
		// initialize by adding or deleting a set of snps from last_model
		// copy from last_model, allocate all, compute X, XX, Xy, R, calculate ll. 
		// set msize, selected, subset but not neighbor_ll. 
		void init_add(class model*, std::vector<int>); 
		void init_add(class model*, int); 
		void init_del(class model*, std::vector<int>); 
		void init_del(class model*, int); 

		// sample from posterior gamma distr; should be called after calc_rss. 
		void sample_tau(void); 
		
		// Calculate Beta using Chol_solve, calculate ll. 
		void calc_rss(void);
		
		int propose(void);
	
		// add_snp(last_model, add_index);
		// del_snp(last_model, del_index) 
		void add_snp(model*, int);
		void del_snp(model*, int); 

		// update g_rb_pip g_rb_beta 
		void rao_blackwell(void);
};

extern std::vector<double> g_proposal_weights; 
extern std::vector<int> g_pass; 

#endif

