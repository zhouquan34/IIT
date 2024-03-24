#! /bin/bash

dir=$1
pre=$2
n=$3
p=$4
k=$5

snrs=(3 2 1) # Setting1: SNR=3; Setting3: SNR=1

for m in {1..3}
	do
	snr=${snrs[$m-1]}
	echo $m $snr
	Rscript sim_data_${pre}.R $k $snr $n $p $dir $m
	./gimh -m $dir/${pre}${k}_$m.mat -p $dir/${pre}${k}_${m}.ph --true $dir/${pre}${k}_${m}.beta -r $k  --kappa 2 --g-exp 3  --start 10 -w 2500000 -s 2500000 -o $dir/out$m/s${k}_r  --rw 
	./iit -m $dir/${pre}${k}_$m.mat -p $dir/${pre}${k}_${m}.ph --true $dir/${pre}${k}_${m}.beta -r $k  --kappa 2 --g-exp 3   --start 10 -s 5000   -o $dir/out$m/s${k}_t0  --ha 0.5
	./iit -m $dir/${pre}${k}_$m.mat -p $dir/${pre}${k}_${m}.ph --true $dir/${pre}${k}_${m}.beta -r $k  --kappa 2 --g-exp 3   --start 10 -s 5000   -o $dir/out$m/s${k}_t1  --hf 1 
	./iit -m $dir/${pre}${k}_$m.mat -p $dir/${pre}${k}_${m}.ph --true $dir/${pre}${k}_${m}.beta -r $k  --kappa 2 --g-exp 3   --start 10 -s 5000   -o $dir/out$m/s${k}_t2  --hf 2
	./iit -m $dir/${pre}${k}_$m.mat -p $dir/${pre}${k}_${m}.ph --true $dir/${pre}${k}_${m}.beta -r $k  --kappa 2 --g-exp 3   --start 10 -s 5000   -o $dir/out$m/s${k}_t3  --hf 3 
	./iit -m $dir/${pre}${k}_$m.mat -p $dir/${pre}${k}_${m}.ph --true $dir/${pre}${k}_${m}.beta -r $k  --kappa 2 --g-exp 3   --start 10 -s 5000   -o $dir/out$m/s${k}_t4  --ha 0.3
	rm $dir/${pre}${k}_${m}.*
done

