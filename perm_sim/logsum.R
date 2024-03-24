log_sum_two <- function(a, b){
	if (a < b){
		return(b + log(1 + exp(a - b)))
	}else{
		return(a + log(1 + exp(b - a)))
	}
}

log_sum <- function(x){
	if (length(x) == 1){return(x)}
	s = log_sum_two(x[1], x[2])
	if (length(x) == 2){return(s)}
	for (i in 3:length(x)){
		s = log_sum_two(s, x[i])
	}
	return(s)
}

