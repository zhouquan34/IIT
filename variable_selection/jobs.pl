#! perl -w

my $dir = 'simA';
my $pre = 'mul';
my $n = 1000; 
my $p = 5000;

my $cpu = 20;
my $each = 5;
foreach my $j (1 .. $cpu){
	open OUT, ">jobs/s$j.sh";
	print OUT  <<HERE; 
#! /bin/bash

HERE
		
	foreach my $m (1 .. $each){
		my $k = ($j-1)*$each + $m; 
		print OUT "bash simulation.sh $dir $pre $n $p $k \n"; 
	}

	close OUT;
		
	system("nohup bash jobs/s$j.sh &");	
}



