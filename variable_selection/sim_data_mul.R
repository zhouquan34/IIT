args <- commandArgs(TRUE)
seed <- as.numeric(args[1])
snr <- as.numeric(args[2])
set.seed(seed)
n <- as.numeric(args[3])
p <- as.numeric(args[4])
dir <- args[5]
m <- as.numeric(args[6])
s = 20
b0 = runif(s, 2, 3) * sample(c(-1,1), s, replace=TRUE)
beta = c(b0, rep(0, p - s)) * sqrt( log(p)/n ) * snr

C = diag(p)
for (i in 1:p){
	for (j in 1:p){
		C[i,j] = exp(-abs(i - j))
	}
}
R = chol(C)

X = matrix(rnorm(n*p), nrow=n) %*% R
X = scale(X, center=TRUE)
y = X %*% beta + rnorm(n)
X = round(X, 6)

pre = paste(dir, "/mul", args[1], "_", args[6], sep="")
true = cbind(1:s, beta[1:s])

write.table(X, file=paste(pre, ".mat", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=" ")
write.table(y, file=paste(pre, ".ph", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(true, file=paste(pre, ".beta", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

