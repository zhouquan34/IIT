#! perl -w

my %count; 
my %sum;
open IN, "time.log";
while (<IN>){
	if (/out(\d+)\/s(\d+)_([risqt]+\d?)\.log\.txt.+=\s+([\d\.]+)\s+s/){
		my $group = "G" . $1 . $3;
		$count{$group} ++;
		$sum{$group} += $4;
		my $t = $4; 
	}	
}
close IN;

my @suff = qw/r t0 t1 t2 t3 t4/; 

	foreach my $k (1 .. 3){
		print "$k";
		foreach my $s (@suff){
			my $a = "G$k$s"; 
			printf "\t%.4f", $sum{$a}/$count{$a};
		}
		print "\n";
	}



