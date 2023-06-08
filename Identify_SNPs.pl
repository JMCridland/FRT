#!perl -w

use strict;

my $snps = shift(@ARGV) or die;#/data/julie/Denovo/loose_files/DGRP.r6.snps.gz
my $pos = shift(@ARGV) or die;#/data/FlyRef/Dmel_r6.41_gene_ranges
my $output = shift(@ARGV) or die;
unlink(qq{$output});

my %keeppos = ();
my $head = "";
my $total = 0;
open(C, "gunzip -c $snps | ");
#chr     pos     ref     alt     line_304        line_307        line_357        line_360        line_399        line_517
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

my $bed = "ALL_SNPs.bed";
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

