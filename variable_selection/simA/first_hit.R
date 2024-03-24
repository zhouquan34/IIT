args <- commandArgs(TRUE)
pre0 <- args[1]

for (j in 1:100){
	stat = c()
	true <- 0
	for (ss in c('r',  't0',  't1', 't2', 't3', 't4')){
		pre = paste(pre0, j, "_", ss, sep="")
		log <- paste(pre, ".log.txt", sep="")
		path <- paste(pre, ".path.txt", sep="")

		a <- readLines(log)
		for (i in 1:length(a)){
			line = a[[i]]	
			ch = unlist(strsplit(line, "\\s+"))
			if (length(ch) < 2){next}
			if (ch[1] == "True" && ch[2] == "ll"){
				true = as.numeric(ch[4])
				break
			}
		}
		sel = 4
		if (ss == 'r'){sel = 5}
		d <- read.table(path, header=FALSE, skip=1)
		m = max(d[,sel])
		diff = m - true 
		k = min(which(d[,sel] >= m - 1e-6)) - 1
		stat = c(stat, m, d[k, 1] + 1) ## add 1 for correction 
	}
	stat = c(j, true, stat)
	cat(stat, "\n")
}


