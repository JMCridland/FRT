#!perl -w

use strict;

my $snps = shift(@ARGV) or die;#DGRP.r6.snps.gz - a tab separated file with the DGRP SNP information just for the lines included in the analysis
my $pos = shift(@ARGV) or die;#A tab separated file with the positions of each gene in the genome
my $output = shift(@ARGV) or die;
unlink(qq{$output});

my %keeppos = ();
my $head = "";
my $total = 0;
open(C, "gunzip -c $snps | ");
while(my $line = <C>){
    chomp $line;
    my @c = split(/\t/, $line);
    if($line =~ m/^chr/){
	$head = $line .  "\t" . "GeneID" . "\t" . "GeneName";	
    }else{
	##first find only the ones that are distinguishing
	my %tmp = ();
	my $count = 0;
	for(my $i = 4; $i <= 9; $i++){
	    if(!(exists($tmp{$c[$i]}))){
		$tmp{$c[$i]} = 1;
		$count++;
	    }		
	}
	if($count > 1){ #keep this one
	    $total++;
	    $keeppos{$c[0] . "\t" . $c[1]} = $line;	    
	}
    }
}
close C;
print "Total distinguishing SNPs = ", $total, "\n";

my $bed = "ALL_SNPs.bed"; ##output bed formatted file
unlink(qq{$bed});
open(A, "<$pos");
open(D, ">>$output");
print D $head, "\n";

my %used = ();
open(E, ">>$bed");
while(my $line = <A>){
    chomp $line;
    my @b = split(/\t/, $line);
    for(my $i = $b[2]; $i <= $b[3]; $i++){
	if(exists($keeppos{$b[1] . "\t" . $i})){
	    if(!(exists($used{$b[1] . "\t" . $i}))){
		print D $keeppos{$b[1] . "\t" . $i}, "\t", $b[0], "\t", $b[5], "\n";
		print E $b[1], "\t", ($i - 1), "\t", $i, "\n";
		$used{$b[1] . "\t" . $i} = 1;
	    }
	}
    }
}
close D;
close E;
close A;

