#!perl -w

use strict;

my $pos = shift(@ARGV) or die;#ALL_SNPs_table - produced by Identify_SNPs.pl
my $exon = shift(@ARGV) or die; #A file with the positions of exons in genes - ex: FBgn0031081     X       19965197        19965511
my $output = shift(@ARGV) or die;
unlink(qq{$output});

my %names = ();
my %pos = ();
my %exon = ();

open(P, "<$pos");
while(my $line = <P>){
    chomp $line;
    my @p = split(/\t/, $line);
    $pos{$p[0] . "\t" . $p[1]} = 1;
}
close P;

my %used = ();
open(E, "<$exon");
while(my $line = <E>){
    chomp $line;
    my @e = split(/\t/, $line);
    for(my $i = $e[2]; $i < $e[3]; $i++){
	if(exists($pos{$e[1] . "\t" . $i})){
	    if(!(exists($used{$e[1] . "\t" . $i . "\t" . $e[0]}))){
		push(@{$names{$e[1] . "\t" . $i}}, $e[0]);
		$used{$e[1] . "\t" . $i . "\t" . $e[0]} = 1;
	    }
	}
    }
}
close E;

open(B, ">>$output");
while((my $k, my $v) = each(%pos)){

    if(!(exists($names{$k}))){
	print B $k, "\t", "None\n";
    }elsif(exists($names{$k})){
	print B $k, "\t", join(",", @{$names{$k}}), "\n";
    }
}
close B;



