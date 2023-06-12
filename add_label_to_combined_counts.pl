#!perl -w

use strict;


my $input = shift(@ARGV) or die;#All_SNPs_with_genes_and_category
my $indir = shift(@ARGV) or die;

my %used = ();
my %gene = (); #genes for a given SNP
my %type = (); #category for a gene at a SNP
my %usedtype = ();

open(A, "<$input");
##X       19961695        FBgn0031081     FBgn0031081     five_prime_UTR
while(my $line = <A>){
    chomp $line;
    my @a = split(/\t/, $line);
    if(!(exists($used{$a[0] . "\t" . $a[1] . "\t" . $a[2]}))){
	$used{$a[0] . "\t" . $a[1] . "\t" . $a[2]} = 1;
	$usedtype{$a[0] . "\t" . $a[1] . "\t" . $a[2] . "\t" . $a[4]} = 1;

	push(@{$gene{$a[0] . "\t" . $a[1]}}, $a[2]);
	push(@{$type{$a[0] . "\t" . $a[1] . "\t" . $a[2]}}, $a[4]);
    }
    if(!(exists($usedtype{$a[0] . "\t" . $a[1] . "\t" . $a[2] . "\t" . $a[4]}))){
	$usedtype{$a[0] . "\t" . $a[1] . "\t" . $a[2] . "\t" . $a[4]} = 1;
	push(@{$type{$a[0] . "\t" . $a[1] . "\t" . $a[2]}}, $a[4]);
    }    
}
close A;


opendir DIR, "$indir";
my @counts = grep{/\.combined\.counts$/} readdir DIR;
closedir DIR;

foreach my $counts (@counts){
    my $out = $counts;
    $out = $out . ".gene_info";
    unlink(qq{$out});
    
    open(C, "<$counts");
    open(D, ">>$out");
    while(my $line = <C>){
	chomp $line;
	my @c = split(/\t/, $line);
	if($line =~ m/pos/){
	    print D $line, "\t", "GeneInfo", "\n";
	}else{
	    if(!(exists($gene{$c[0] . "\t" . $c[1]}))){
		print D $line, "\t", "NA\n";
	    }elsif(exists($gene{$c[0] . "\t" . $c[1]})){
		print D $line, "\t";

		my @tmp = ();
		foreach my $g (@{$gene{$c[0] . "\t" . $c[1]}}){
		    my $set = $g . ":" . join(",", @{$type{$c[0] . "\t" . $c[1] . "\t" . $g}});
		    push(@tmp, $set);
		}		
		print D join(",", @tmp), "\n";
	    }
	}
    }
    close C;
    close D;
}
