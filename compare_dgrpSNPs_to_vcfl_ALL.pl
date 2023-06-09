#!perl -w

use strict;

my $input = shift(@ARGV) or die; #Output file from Identify_SNPs.pl
my $dir = shift(@ARGV) or die; #specify the directory with the vcf files

my %diffsnps = ();
my %conv = ();
my %id = ();
my %snps = (); 
open(A, "<$input");
while(my $line = <A>){
    chomp $line;
    ##chr     pos     ref     alt     line_304        line_307        line_357        line_360        line_399        line_517        GeneID     GeneName
    my @a = split(/\t/, $line);
    if($line =~ m/^chr/){
	for(my $i = 4; $i <= 9; $i++){
	    $a[$i] =~ s/line_//;
	    $conv{$i} = $a[$i]; 
	    print $conv{$i}, "\n";
	}
    }else{
	for(my $j = 4; $j <= 8; $j++){
	    if(($a[$j] ne "-") and ($a[9] ne "-")){
		if($a[$j] ne $a[9]){ #there is a difference between the line and 517
		    $diffsnps{$conv{$j} . "\t" . $conv{9} . "\t" . $a[0] . "\t" . $a[1]} = $a[$j] . "\t" . $a[9];
		    $id{$a[0] . "\t" . $a[1]} = $a[10] . "\t" . $a[11];
		   # $snp{$a[0] . "\t" . $a[1]} = $a[2] . "\t" . $a[3];
		}
	    }
	}
    }
}
close A;

opendir DIR, "$dir";
my @vcf = grep{/_ALL\.vcf$/} readdir DIR;
closedir DIR;

foreach my $vcf (@vcf){
    my @v = split(/\./, $vcf);
    
    my $out = $v[0] . "_distinguishingSNPs";
    unlink(qq{$out});
    
    my $id = $v[0];
    $id =~ s/_ALL//;
    $id =~ s/_\D+//; 
    
    open(X, ">>$out");
    print X $id, "\t", "517", "\t", "chr", "\t", "pos", "\t", $id, "_SNP", "\t", "517_SNP", "\t", "vcfSNP", "\t", "Ref_cov", "\t", "Alt_cov", "\t", "RefSNP", "\t", "AltSNP", "\n"; # "\t", "GeneID", "\t", "GeneName", "\n";
    open(V, "<$vcf");
    while(my $line = <V>){
	chomp $line;
	if($line !~ m/^#/){
	    my @x = split(/\t/, $line);
	    if(exists($diffsnps{$id . "\t" . "517" . "\t" . $x[0] . "\t" . $x[1]})){ 
		$snps{$id . "\t" . "517" . "\t" . $x[0] . "\t" . $x[1]} = $x[3] . "\t" . $x[4];
		my @y = split(/;/, $x[7]);
		
		my $which = ();
		for(my $i = 0; $i < scalar(@y); $i++){
		    if($y[$i] =~ m/DP4=/){
			$which = $y[$i];
		    }
		}
		if($which =~ m/DP4=/){
		    $which =~ s/DP4=//;
		    my @w = split(/,/, $which);
		    
		    my @tmp = ();
		    my $ref = $w[0] + $w[1];
		    my $alt = $w[2] + $w[3];
		    
		    if($ref > 0){
			push(@tmp, $x[3]);
		    }
		    if($alt > 0){
			push (@tmp, $x[4]);
		    }
		    if(($alt > 0) or ($ref > 0)){
			print X $id, "\t", "517", "\t", $x[0], "\t", $x[1], "\t", $diffsnps{$id . "\t" . "517" . "\t" . $x[0] .  "\t" . $x[1]}, "\t", join(",", @tmp), "\t", $ref, "\t", $alt, "\t", $snps{$id . "\t" . "517" . "\t" . $x[0] . "\t" . $x[1]}, "\n"; #"\t", $id{$x[0] . "\t" . $x[1]}, "\n";
		    }
		}
	    }
	}
    }
    close V;
    close X;
}
