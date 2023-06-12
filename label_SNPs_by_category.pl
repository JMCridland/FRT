#!perl -w

use strict;

my $input = shift(@ARGV) or die; ##All_SNPs_with_gene_list
my $flydir = shift(@ARGV) or die; #/data/FlyRef/Drosophila_melanogaster/6.41/fasta/
my $output = shift(@ARGV) or die;
unlink(qq{$output});


my %pos = ();

##2R      4993342 None
open(A, "<$input");
while(my $line = <A>){
    chomp $line;
    my @a = split(/\t/, $line);
    $pos{$a[0] . "\t" . $a[1]} = $a[2];
}
close A;

opendir DIR, "$flydir";
my @range = grep{/\.ranges/} readdir DIR;
closedir DIR;


open(B, ">>$output");
foreach my $range (@range){
    $range = $flydir . $range;
    open(R, "<$range");
    while(my $line = <R>){
	chomp $line;
	##FBgn0031081     FBgn0031081     intron  X       19961846        19963954
	my @x = split(/\t/, $line);

	for(my $i = $x[4]; $i <= $x[5]; $i++){
	    if(exists($pos{$x[3] . "\t" . $i})){
		print B $x[3], "\t", $i, "\t", $x[0], "\t", $x[1], "\t", $x[2], "\n";
	    }
	}
    }
    close R;
}
close B;
