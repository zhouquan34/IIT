#! perl -w

args = commandArgs(TRUE)
dir = args[1]
s = as.numeric(args[2])
n = as.numeric(args[3])
step = as.numeric(args[4])
p = as.numeric(args[5])
pick = seq(step, n, by = step)
np = length(pick)
sum_x = matrix(0, nrow=np, ncol=p)
sum_xx = matrix(0, nrow=np, ncol=p)
nr = 100
cat("")
cols = (1:p)+6
if (s == 0) {cols = (1:p)+4}

for (i in 1:nr){
	f = paste(dir, "/r", i, "_s", s, ".txt", sep="")
	if (s == 0){
		f = paste(dir, "/r", i, "_rw.txt", sep="")
	}
	cat("Reading", f, "\r")
	d = as.matrix(read.table((f)))
	d1 = d[pick, cols]
	sum_x = sum_x + d1
	sum_xx = sum_xx + d1^2
}
cat("\n")

V = sum_xx / nr  - (sum_x / nr)^2
mean = signif(sum_x[np, ]/nr, 6)
cat(dir, s, mean, "\n", append=TRUE, file='means.txt')
cat(signif(as.numeric(summary(V[np, ])),3), "\n")
f = paste("summary/", dir, "_s", s, ".txt", sep="")
write.table(V, file=f, quote=F, col.names=F, row.names=F)

