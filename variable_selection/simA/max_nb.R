options(digits=12)
args <- commandArgs(TRUE)
pre0 <- args[1]

lab = c('t0',  't1', 't2', 't3', 't4')
l = length(lab)
outd = c()
nmode = numeric(100)
nu = 0
for (j in 1:100){
	dat = c()
	for (k in 1:l){
		ss = lab[k]
		pre = paste(pre0, j, "_", ss, sep="")
		err = paste(pre, ".path.txt", sep="")
		d = read.table(err, header=F, skip=1)
		d = d[,c(4,7,8,9)]
		d = d[-nrow(d), ]
		score = d[,3] * 1e8 + d[,4] * 1e6 + d[,1] + d[,2] * 1e4
		d = cbind(d, score)
		dat = rbind(dat, d)	
	}
	s = dat[,5]
	uni = split(seq_along(s), s)
	nu = nu + length(uni)
	u = sapply(uni, "[[", 1)
	td = dat[u, 2:4]
	nmode[j] = sum(td[,1] < 0)
	outd = rbind(outd, td)
}

out = args[2]
nmode = as.matrix(nmode)
write.table(outd, file = paste(out, "_dn.txt", sep=""), quote=F, col.names=F, row.names=F)
write.table(nmode, file = paste(out, "_nms.txt", sep=""), quote=F, col.names=F, row.names=F)

