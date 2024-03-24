#include "args.h"

using namespace std;

string true_file = "";

void print_help(void){
	cout << "Example command: " << endl;
	cout << "\t\t./iit -m test.mat -p test.ph -o out -s 10000" << endl;
	cout << endl;
	cout << "Arguments:" << endl;
	cout << "-m: n x p design matrix file" << endl;
	cout << "-p: response file" << endl;
	cout << "-o: output prefix" << endl;
	cout << "-w: number of burn-in iterations" << endl;
	cout << "-s: number of sampling iterations of this run" << endl;
//	cout << "-R: do Rao-Blackwellization every R iterations" << endl;
	cout << "-r: random seed" << endl;
	cout << endl;
//	cout << "By default, gimh uses upper_add = 2, upper_del = 1, lower_add = lower_del = 0" << endl;
	exit(EXIT_FAILURE);
	return;
}


void print_manual(void){
	cout << "All arguments:" << endl;
	cout << "--mg:        input mean genotype file" << endl;
	cout << "--plink:     input plink A-transpose genotype file" << endl;
	cout << "--mat:       input matrix genotype file" << endl;
	cout << "--pheno:     input phenotype file" << endl;
	cout << "--out:       output prefix" << endl;
//	cout << "--last:      output prefix of the dataset to continue" << endl;
	cout << "--mcmc:      MCMC iteration" << endl;
	cout << "--burn:      burn-in iteration" << endl;
	cout << "--rb:        Rao-Blackwell frequency" << endl;
	cout << "--start:     (-nstart) starting model size" << endl;
//	cout << "--lunif-h2:  use uniform prior on log-h2" << endl;
	cout << "--time:      output time usage details" << endl;
//	cout << "--long:      long-range proposal proportion" << endl;
//	cout << "--jump:      max jump size" << endl;
//	cout << "--h2-walk:   max h2 walk distance" << endl;
//	cout << "--icf-err:   ICF max error" << endl;
//	cout << "--pi-alpha:  prior alpha of pi" << endl;
//	cout << "--pi-beta:   prior beta of pi" << endl;
//	cout << "--pi-min:    min pi" << endl;
//	cout << "--pi-max:    max pi" << endl;
	cout << "--prec:      dim of XtX to be precalculated" << endl;
	cout << "--max-size:  (-smax) max model size allowed" << endl;
	cout << "--min-size:  (-smin) min model size allowed" << endl;
//	cout << "--exact:     exact calculation of BF" << endl;
	cout << "--seed:      GSL seed" << endl;
	cout << "--help:      print help" << endl;
	cout << "--HELP:      print full help" << endl;
//	cout << "--add-max:   cutoff value for addition weights" << endl;
//	cout << "--add-min:   cutoff value for addition weights" << endl; 
//	cout << "--del-max:   cutoff value for deletion weights" << endl;
//	cout << "--del-min:   cutoff value for deletion weights" << endl;
//	cout << "--rw:        no weighting and use random walk proposals" << endl;
//	cout << "--g-exp:     g prior with g = p^c" << endl; 
	cout << "--ostart:    starting model size" << endl; 
	cout << "--init:      initial model file" << endl; 
//	cout << "--approx:    approximated informed" << endl; 
//	cout << "--omit:      file listing snps to be skipped" << endl; 
//	cout << "--lb:        use square-root locally balanced proposals" << endl; 
	cout << "--hf         function h in IIT" << endl; 
	cout << "--ha         constant a in h(u)=u^a" << endl;
	exit(EXIT_FAILURE);
	return;
}


void read_start(string file){
	cerr << "Reading initial model file " << file << endl;
	ifstream DAT(file.c_str() ); 
	string s; 
	while (getline(DAT, s))	{
		stringstream ss(s);  int k;  
		if (ss >> k){ g_init_model.push_back(k - 1);}
	}
	DAT.close();	
	return;
}



void init_paras(void){
	g_paras["-g"] = PARA_MG;
	g_paras["--mg"] = PARA_MG;
	g_paras["-a"] = PARA_PLINK;
	g_paras["--plink"] = PARA_PLINK;
	g_paras["-m"] = PARA_MAT;
	g_paras["--mat"] = PARA_MAT;
	g_paras["-p"] = PARA_PHENO;
	g_paras["--pheno"] = PARA_PHENO;
	g_paras["-o"] = PARA_OUT;
	g_paras["--out"] = PARA_OUT;
//	g_paras["-c"] = PARA_COV;
//	g_paras["--cov"] = PARA_COV;
	g_paras["-s"] = PARA_ITER;
	g_paras["--mcmc"] = PARA_ITER;
	g_paras["-w"] = PARA_BURN;
	g_paras["--burn"] = PARA_BURN;
	g_paras["-R"] = PARA_RB;
	g_paras["--rb"] = PARA_RB;
	g_paras["-k"] = PARA_START;
	g_paras["--start"] = PARA_START;
	g_paras["--ostart"] = PARA_SFIX;
	g_paras["--init"] = PARA_INIT;
	g_paras["-t"] = PARA_TIME;
	g_paras["--time"] = PARA_TIME;
	g_paras["-b"] = PARA_BVSR;
	g_paras["--last"] = PARA_BVSR;
	g_paras["-z"] = PARA_SEED;
	g_paras["-r"] = PARA_SEED;
	g_paras["--seed"] = PARA_SEED;
	g_paras["-h"] = PARA_HELP;
	g_paras["--help"] = PARA_HELP;	
	g_paras["--long"] = PARA_LONG;
	g_paras["--geom"] = PARA_GEOM;
	g_paras["--icf-err"] = PARA_ICFE;
	g_paras["--pi-alpha"] = PARA_PI_A;
	g_paras["--pi-beta"] = PARA_PI_B;
	g_paras["--pi-min"] = PARA_MINPI; 
	g_paras["--pi-max"] = PARA_MAXPI; 
	g_paras["--add-max"] = PARA_MAXADD; 
	g_paras["--add-min"] = PARA_MINADD; 
	g_paras["--del-max"] = PARA_MAXDEL; 
	g_paras["--del-min"] = PARA_MINDEL; 
	g_paras["--rw"] = PARA_RW; 
	g_paras["--g-exp"] = PARA_G; 
	g_paras["--h2-min"] = PARA_MINH2;
	g_paras["--h2-max"] = PARA_MAXH2;
	g_paras["--jump"] = PARA_MAXJ;
	g_paras["--lunif-h2"] = PARA_LOGH;
	g_paras["--h2-walk"] = PARA_HJUMP;
	g_paras["--prec"] = PARA_PREC;
	g_paras["--max-size"] = PARA_MAXS;
	g_paras["--min-size"] = PARA_MINS;
	g_paras["--exact"] = PARA_EXACT;
	g_paras["--kappa"] = PARA_KAPPA;
	g_paras["--HELP"] = PARA_MAN;
	g_paras["--true"] = PARA_TRUE;
	g_paras["-pmax"] = PARA_PMAX;
	g_paras["-pmin"] = PARA_PMIN;
	g_paras["--approx"] = PARA_AI;		
	g_paras["--omit"] = PARA_OMIT; 	
	g_paras["--lb"] = PARA_LB; 	
	g_paras["--hf"] = PARA_H;	
	g_paras["--ha"] = PARA_H_A;	
	
	// PIMASS parameters //
	g_paras["-nstart"] = PARA_START;
	
	// PIMASS parameteres not defined //
	g_paras["-num"] = PARA_RB1;
	g_paras["-hmin"] = PARA_MINH2;
	g_paras["-hmax"] = PARA_MAXH2;
	g_paras["-smax"] = PARA_MAXS;
	g_paras["-smin"] = PARA_MINS;
	g_paras["-pos"] = PARA_NOT;
	g_paras["-cc"] = PARA_NOT;
	g_paras["-v"] = PARA_NOT2;
	g_paras["-exclude-maf"] = PARA_NOT2;
	g_paras["-exclude-nopos"] = PARA_NOT2;
	g_paras["-silence"] = PARA_NOT2;
	return;
}

void print_continue_help(void){
	return;
/*	cerr << "The continuation mode only accepts the following arguments: \n" << endl;
	cerr << "-b: prefix of the MCMC to continue" << endl;
	cerr << "-o: output prefix" << endl;
	cerr << "-s: number of MCMC iterations of this run" << endl;
	cerr << "-z: random seed" << endl;
	cerr << "\nAll the other arguments are forced to take their values in the first run" << endl;
	exit(EXIT_FAILURE);
*/
}

int set_args_part_one (string para, char* val){
	if (g_paras.find(para) == g_paras.end()){ return 0;}
	switch(g_paras[para]){
		case PARA_BVSR:
			g_last_bvsr = val; return 1;
		case PARA_SEED:
		//	cerr << "WARN: -R is for Rao-Blackwellization and -r is for random seed" << endl;
			gsl_seed = atoi(val); return 1;
		case PARA_OUT:
			g_out_prefix = val; return 1;
		case PARA_ITER:
			g_mcmc_iter = atoi(val); return 1;	
		default:
			return 0;
	}
	return 0;
}

int set_args_part_two (string para, char* val){
	if (g_paras.find(para) == g_paras.end()){ return 0;}
	switch(g_paras[para]){
		case PARA_MG:
			read_meang(val); return 1;
		case PARA_PLINK:
			read_plink(val); return 1;
		case PARA_MAT:
			read_matrix(val); return 1;
		case PARA_PHENO:
			read_pheno(val); return 1;
		case PARA_TRUE:
			true_file = val; g_true = 1; return 1;
		case PARA_ITER:
			g_mcmc_iter = atoi(val); return 1;  // may be g_last_iter
		case PARA_BURN:
			g_mcmc_warm = atoi(val); return 1;
		case PARA_RB:
			g_rb_thin = atoi(val); return 1;
		case PARA_RB1:
			g_rb_thin = atoi(val) * 10; return 1;
		case PARA_START:
			g_start_size = atoi(val); return 1; 
		case PARA_SFIX:
			g_start_size = atoi(val); g_ordered_start = 1; return 1; 	
		case PARA_TIME:
			g_output_time = 1; return 1;
		case PARA_LOGH: // turn on
			g_log_uniform_h = 1; return 1;
		case PARA_LONG:	
			g_long_range = atof(val); return 1;
		case PARA_ICFE:
			g_icf_abs_tol = atof(val); return 1; 
		case PARA_MINPI:
			g_min_pi = atof(val); g_pi_range_set = 1; return 1;
		case PARA_MAXPI:
			g_max_pi = atof(val); g_pi_range_set = 1; return 1;
		case PARA_G:
			g_gg_exp = atof(val); return 1;
		case PARA_MAXADD:
			g_upper_add = atof(val); g_rw = 0; return 1;
		case PARA_MINADD:
			g_lower_add = atof(val); g_rw = 0; return 1;
		case PARA_MAXDEL:
			g_upper_del = atof(val); g_rw = 0; return 1;
		case PARA_MINDEL:
			g_lower_del = atof(val); g_rw = 0; return 1;
		case PARA_RW:
			g_rw = 1; return 1;
		case PARA_LB:
			g_lbsq = 1; return 1;
		case PARA_AI:
			g_approx_informed = 1; return 1;
		case PARA_KAPPA:
			g_kappa = atof(val); return 1;
		case PARA_PMAX:
			g_pmax = atof(val); g_pset = 1; return 1;
		case PARA_PMIN:
			g_pmin = atof(val); g_pset = 1; return 1;		
		case PARA_MAXJ:
			g_max_jump = atoi(val); return 1;
		case PARA_HJUMP:
			g_jump_h2 = atof(val); return 1;
		case PARA_PI_A:
			g_pi_prior_alpha = atoi(val); return 1;
		case PARA_PI_B:
			g_pi_prior_beta = atoi(val); return 1;
		case PARA_PREC:
			g_precomp = atoi(val); return 1;
		case PARA_MAXS:
			g_max_model_size = atoi(val); return 1;
		case PARA_MINS:
			g_min_model_size = atoi(val); return 1;
		case PARA_EXACT:
			g_exact_bf = 1; return 1;
		case PARA_INIT:
			read_start(val); return 1; 
		case PARA_OMIT:
			g_skip = 1; read_omit(val); return 1; 
		case PARA_H:
			g_h_func = atoi(val);
			return 1;
		case PARA_H_A:
			g_h_A = atof(val);
			return 1;
		case PARA_HELP:
			print_help(); return 1;
		case PARA_MAN:
			print_manual(); return 1;	
		case PARA_NOT:
			cerr << "NOTE: -pos, -cc are not enabled in fastBVSR" << endl; return 1;
		case PARA_NOT2:
			cerr << "NOTE: -v, -silence, -exclude-maf, -exclude-nopos are not enabled in fastBVSR"<< endl; return 1;
		default:
			return 0;
	}
	return 0;
}


int scan_continue(int argc, char** argv){
	for (int i=1; i<argc; i++){
		if (argv[i][0] == '-'){
			string para = argv[i];
			if (g_paras.find(para) != g_paras.end()){
				if (g_paras[para] == PARA_BVSR){
					return 1;
				}	
			}
		}
	}
	return 0;
}

void read_continue_args (int argc, char** argv){
	for (int i=1; i<argc; i++){
		if (argv[i][0] == '-'){
			string para = argv[i];
			char* val = NULL;
			if (i < argc - 1) { val = argv[i+1]; }
			if (set_args_part_one(para, val) == 0){
				cerr << "Argument " << para << " is invalid!" << endl;
				print_help();
			}
		}
	}
	return;
}

void read_last_args (void){
	string log_file = g_last_bvsr;
	log_file.append(".log.txt");		
	ifstream flog(log_file.c_str());
	if (flog.fail()){
		cerr << "Log file of last run does not exist!" << endl;
		exit(EXIT_FAILURE);
	}

	string s;
	int mcmc_iter_copy = g_mcmc_iter; // already set
//	int gsl_seed_copy  = gsl_seed;
	g_mcmc_iter = 1000; // default value
	while (getline(flog, s)){
		stringstream ss(s); 
		vector<string> args;	
		char pound; ss >> pound;  
		if (pound != '#'){break;} // scan all the previous commands
		if (s.empty()){break;}
		string tmp; ss >> tmp;
		while (ss >> tmp){ args.push_back(tmp); }
		LOG << "# Previous: "; 
		for (int i=0; i<args.size(); i++){ LOG << args[i] << " ";} 
		LOG << endl;
		
		for (int i=1; i<args.size(); i++){
			if (args[i].at(0) != '-'){ continue; }
			char* val = NULL;
			if (i < args.size() - 1) { val = const_cast<char*> (args[i+1].c_str()); }
			set_args_part_two(args[i], val);
		}
		g_last_iter += g_mcmc_iter; 
	}
	flog.close();
	g_mcmc_iter = mcmc_iter_copy; // recover
//	gsl_seed = gsl_seed_copy;
	return;
}


void read_args (int argc, char** argv){
	if (argc <= 1){ print_help(); }
	init_paras();
	if (scan_continue(argc, argv) == 1){
		g_continue = 1;
		read_continue_args(argc, argv);
		read_last_args();
		return;
	}

	for (int i=1; i<argc; i++){
		if (argv[i][0] == '-' && (!isdigit(argv[i][1])) ){
			string para = argv[i];
			char* val = NULL;
			if (i < argc - 1) { val = argv[i+1]; }
			int defined = 0; 
			defined += set_args_part_one(para, val);
			defined += set_args_part_two(para, val);
			if (defined == 0){		
				print_help();
			}
		}
	}
	read_true(true_file);	
	return;
}


