#! perl -w

my @set;
open IN, "setting";
while (<IN>){
	if (/^#/){next;}
	my @col = split /\s+/;
	push @set, \@col; 
}
close IN;
my $ns = $#set;

foreach my $k (1 .. 20){
	my $s = ($k-1)*5+1;
	my $e = $s + 4;
	open OUT, ">jobs/j$k.sh";
	print OUT '#! /bin/bash' . "\n";
	print OUT 'for i in {' . $s . '..' . $e . "} \n";
	print OUT "\tdo \n";
	
	foreach my $m (0 .. $ns){
		my $p = $set[$m]->[0];
		my $n1 = $set[$m]->[1];
		my $n2 = $set[$m]->[2];
		my $scheme = $set[$m]->[3];
		my $snr = $set[$m]->[4];
		my $thin = $set[$m]->[5];
		my $head = "Rscript sample_perm.R \$i out$m/r\$\{i\}_";
		my $head_rw = "Rscript sample_perm_rw.R \$i out$m/r\$\{i\}_";
		my @cmd;
		push @cmd,  $head . "s1.txt $p $n1 1 $snr $scheme";
		push @cmd,  $head . "s2.txt $p $n1 2 $snr $scheme";
		push @cmd,  $head . "s3.txt $p $n1 3 $snr $scheme 0.5";
		push @cmd,  $head . "s4.txt $p $n1 3 $snr $scheme 0.4";
		push @cmd,  $head . "s5.txt $p $n1 3 $snr $scheme 0.3";
		#push @cmd,  $head . "s6.txt $p $n1 3 $snr $scheme 0.6";
		push @cmd,  $head_rw . "rw.txt $p $n2 $thin $snr $scheme";
		foreach my $cmd (@cmd){
			print OUT "\t$cmd\n";
		}		
	}
	print OUT "done\n";
	close OUT;
	#system("nohup bash jobs/j$k.sh &");
}

