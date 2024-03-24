args <- commandArgs(TRUE)

A0 = 2e6
A1 = 5000
qs = c(0.025, 0.975)
null = -4.00002e-9

qstring <- function(x){
	m = round(median(x))
	if (m == A0){m = "$2e6+$"}
	qx = round(quantile(x, 0.95))
	if (qx == A0){qx = "$2e6+$"}
	if (qx == A1){qx = "$5000+$"}
	ss = paste(m, " (", qx,  ")", sep="")
	return(ss)
}

for (j in 1:3){
	out=matrix(0, nrow=100, ncol=6)
	file = paste('summary/A', j,  '.hit',  sep="")
	d <- read.table(file)
	true <- d[,2]	
	pick = seq(3,13,by=2)
	ns = length(pick)
	d0 = d[,pick]
	high = as.numeric(apply(d0, 1, max))
	stat = character(ns)
	for (s in 1:ns){
		fail = (d0[,s] < high)
		iter = d[,pick[s] + 1]
		if (s == 1){
			iter[fail] = A0
		}else{
			iter[fail] = A1
		}
		stat[s] = qstring(iter)
		out[,s] = iter
	}
	write.table(out, file = paste('summary/A', j, '_hit.txt', sep=''), quote=F, row.names=F, col.names=F)
	cat(j, sum(abs(high-true)<1e-10), sum(abs(high-null)<1e-10), stat, sep=" ")
	cat("\n")
}



