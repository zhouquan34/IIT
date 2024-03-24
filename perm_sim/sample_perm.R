source("logsum.R")
args = commandArgs(TRUE)
seed = as.numeric(args[1]) 
out.file = args[2]  
p = as.numeric(args[3])
n.mcmc = as.numeric(args[4])
h.ind = as.numeric(args[5])
snr = as.numeric(args[6])
scheme = as.numeric(args[7])
A = 0.5
log.h = function(u){stop("h is not defined!")}
if (length(args) > 7){
	A = as.numeric(args[8])
}

######### set function h #########
if (h.ind == 1){
	log.h = function(u){
		if (u < 0){return(log(1 + exp(u)))}	
		return(u + log(1 + exp(-u)))
	}
}
if (h.ind == 2){
	log.h = function(u){
		if (u < 0){return(u)}	
		return(0)
	}
}
if (h.ind == 3){
	log.h = function(u){return(u * A)}
}
##################################

######## simulate posterior ######
log.f <- function(k, phi){ g = -abs((1:p)-k) * log(p) * phi }
set.seed(34)  # this seed is always fixed 
W = matrix(0, p, p)
for (i in 1:p){
	if (scheme < 3){
		mean = runif(1, max(1, i-0.1), min(p, i+0.1) )
		sigma = runif(1, 0.5, 1) * snr
		W[i, ] = log.f(mean, sigma)	
	}else{
		mean = runif(1, max(1, i-0.5), min(p, i+0.5) )
		sigma = runif(1, 0.1, 1) * snr
		W[i, ] = log.f(mean, sigma)
	}
}
if (scheme == 2){
	W[1, ] = log.f(1, p*snr)
	W[2, ] = W[1, ]
}
#################################


###### informed proposal ###########
comb = combn(p, 2)
nb = ncol(comb)
# s[k]: pos of k; inverse of tau
calc_proposal = function(s){
	prop = numeric(nb)
	ll = numeric(nb)
	for (j in 1:nb){
		a = comb[1, j]
		b = comb[2, j]
		ll[j] = -W[a, s[a]] - W[b, s[b]] + W[a, s[b]] + W[b, s[a]]
		prop[j] = log.h(ll[j])
	}
	log.Z = log_sum(prop) 
	prop = exp(prop - max(prop))
	prop = prop/sum(prop)
	return(list("prop"=prop, "ll"=ll, "Z"=log.Z))
}

########### MCMC #############
set.seed(seed)
s = sample(1:p, p) 
log.pi = 0
for (i in 1:p){ log.pi = log.pi + W[i, s[i]] }
match = sum(s == (1:p))
K = calc_proposal(s)

out = matrix(0, nrow=n.mcmc, ncol=p+6)
sum_w = -1e15
sum_wx = rep(-1e15, p)
start_t = Sys.time()
for (t in 1:n.mcmc){
	k = sample(1:nb, 1, prob=K$prop)
	pb = K$prop[k]
	a = comb[1, k]
	b = comb[2, k]
	match = match + (s[a] == b) + (s[b] == a) - (s[a] == a) - (s[b] == b)
	tmp = s[b]; s[b] = s[a]; s[a] = tmp
	
	log.pi = log.pi + K$ll[k]
	K = calc_proposal(s)
	omega = (1 - 2*A) * log.pi - K$Z
	for (i in 1:p){ sum_wx[i] = log_sum_two(sum_wx[i], log(s[i]) + omega)}
	sum_w = log_sum_two(sum_w, omega)
	## output: iteration #, proposal prob, # of matches, log(pi), log(Z), log(omega), mean rank.
	out[t, ] = signif(c(t, pb, match, log.pi, K$Z, omega, exp(sum_wx - sum_w)), 12)
}

end_t = Sys.time()
time = difftime(end_t, start_t, units="secs")
cat(args, time, "\n", append=TRUE, file = 'sim.log')
write.table(out, file=out.file, quote=F, row.names=F, col.names=F)

