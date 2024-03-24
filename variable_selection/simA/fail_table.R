qstring <- function(x){
	m = signif(mean(x), 5)
	qx = signif(quantile(x, qs), 5)
	ss = paste(m, " (", qx[1], ", ", qx[2], ")", sep="")
	return(ss)
}

for (j in 1:3){
	file = paste('summary/A', j,  '.hit',  sep="")
	d <- read.table(file)
	pick = c(3,5,7,9, 11, 13)
	ns = length(pick)
	d0 = d[,pick]
	high = as.numeric(apply(d0, 1, max))
	fail = numeric(ns)
	for (s in 1:ns){
		fail[s] = sum(d0[,s] < high)
	}
	cat(j, fail, sep="\t")
	cat("\n")	
}

