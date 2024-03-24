#! /bin/bash

grep "Time used" out*/*log.txt >time.log
perl parse_time.pl

for i in {1..3}
	do
	echo $i
	Rscript max_nb.R out$i/s summary/A$i 
	Rscript first_hit.R out$i/s >summary/A${i}.hit
done


