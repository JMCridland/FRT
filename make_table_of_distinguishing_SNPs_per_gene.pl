#!perl -w

use strict;

my $input = shift(@ARGV) or die; #output from Identify_genes_for_each_site.pl
my $indir = shift(@ARGV) or die;

my %pos = ();
##2L      3151519 None

open(A, "<$input");
##X       14920903        FBgn0030598
while(my $line = <A>){
    chomp $line;
    my @a = split(/\t/, $line);
    $pos{$a[0] . "\t" . $a[1]} = $a[2];
}
close A;

###For each cross I want a list of genes and the number of potential distinguishing SNPs

opendir DIR, "$indir";
my @snps = grep {/ALL_distinguishingSNPs$/} readdir DIR;
closedir DIR;

opendir DIR2, "$indir";
my @masked = grep{/_Diff_masked$/} readdir DIR2;
closedir DIR2;

opendir DIR3, "$indir";
my @same = grep{/_Same_masked$/} readdir DIR3;
closedir DIR3;

my %maskpos = ();###all positions per cross that are masked

foreach my $mask (@masked){
    my @x = split(/_/, $mask);
    ##399     517     X       316070  Masked
    open(M, "<$mask");
    while(my $line = <M>){
	chomp $line;
	my @m = split(/\t/, $line);
	if($m[4] eq "Masked"){
	    $maskpos{$x[0] . "\t" . $m[2] . "\t" . $m[3]} = 1;
	}
    }
    close M;
}
foreach my $same (@same){
    my @x = split(/_/, $same);
    ##399     517     X       316070  Masked
    open(S, "<$same");
    while(my $line = <S>){
	chomp $line;
	my @s = split(/\t/, $line);
	if($s[4] eq "Masked"){
	    $maskpos{$x[0] . "\t" . $s[2] . "\t" . $s[3]} = 1;
	}
    }
    close S;
}

my %cross = ();
my @cross = (304,307,360,399);
my %used = ();
my %gene = ();
foreach my $snps (@snps){
    print $snps, "\n";
    my @s = split(/_/, $snps);

    my $id = $s[0]; ##the cross

    open(A, "<$snps");
    while(my $line = <A>){
	chomp $line;
	my @a = split(/\t/, $line);
	if($line !~ m/chr/){
	    #304     517     chr     pos     304_SNP 517_SNP vcfSNP  Ref_cov Alt_cov RefSNP  AltSNP
	    if(($a[7] + $a[8]) > 0){
		my @gene = split(/,/, $pos{$a[2] . "\t" . $a[3]});

		foreach my $g (@gene){
		    $gene{$g} = 1;
		    if(!(exists($maskpos{$id . "\t" . $a[2] . "\t" . $a[3]}))){ ##if it isn't masked
			if(!(exists($used{$id . "\t" . $a[2] . "\t" . $a[3] . "\t" . $g}))){ #for a particular cross position gene set
			    $used{$id . "\t" . $a[2] . "\t" . $a[3] . "\t" . $g} = 1;
			    if(!(exists($cross{$id . "\t" . $g}))){ #add to the number of SNPs for that gene
				$cross{$id . "\t" . $g} = 1;
			    }elsif(exists($cross{$id . "\t" . $g})){
				$cross{$id . "\t" . $g}++;
			    }
			}
		    }
		}
	    }
	}
    }
    close A;
}

foreach my $cross (@cross){
    
    my $output = $cross . ".genes_with_distinguishingSNPs";
    unlink(qq{$output});
    open(X, ">>$output");
    while((my $k, my $v) = each(%gene)){
	if(!(exists($cross{$cross . "\t" . $k}))){
	    $cross{$cross . "\t" . $k} = 0;
	}
	print X $k, "\t", $cross{$cross . "\t" . $k}, "\n";
    }
    close X;
}
