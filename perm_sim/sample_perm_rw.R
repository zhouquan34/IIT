args = commandArgs(TRUE)
seed = as.numeric(args[1])
out.file = args[2]
p = as.numeric(args[3])
n.mcmc = as.numeric(args[4])
thin = as.numeric(args[5]) ### different from IIT
snr = as.numeric(args[6])
scheme = as.numeric(args[7])


## all on log-scale
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


set.seed(seed)
s = sample(1:p, p)
log.pi = 0
for (i in 1:p){
	log.pi = log.pi + W[i, s[i]]
}
match = sum(s == (1:p))
comb = combn(p, 2)
nb = ncol(comb)

out = matrix(0, nrow=n.mcmc/thin, ncol=p+4)
sum_wx = rep(0, p)
nacc = 0 
start_t = Sys.time()
for (t in 1:n.mcmc){
	k = sample(1:nb, 1)
	a = comb[1, k]
	b = comb[2, k]
	s1 = s
	s1[a] = s[b]
	s1[b] = s[a]
	ll = -W[a, s[a]] - W[b, s[b]] + W[b, s[a]] + W[a, s[b]]
	acc = 0 
	if (ll > 0){acc = 1}
	else{
		if (runif(1) < exp(ll)){acc = 1}
	}
	if (acc == 1){
		log.pi = log.pi + ll
		match = match + (s[a] == b) + (s[b] == a) - (s[a] == a) - (s[b] == b)
		s = s1		
	}	
	nacc = nacc + acc	
	for (i in 1:p){
		sum_wx[i] = sum_wx[i] + s[i]
	}
	if (t %% thin == 0){
		out[t/thin, ] = signif(c(t, k, match, log.pi, sum_wx/t),10)
	}
}
end_t = Sys.time()
time = difftime(end_t, start_t, units="secs")
cat(args, nacc, time, "\n", append=TRUE, file = 'sim_rw.log')
write.table(out, file=out.file, quote=F, row.names=F, col.names=F)


