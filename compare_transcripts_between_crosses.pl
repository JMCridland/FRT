#!perl -w

use strict;

##I want to make a table that combines the information across all crosses and tissues for each gene/SNP

##I need to do this at the SNP level
#I'd like to compare within a cross - 304 (three tissues) 307 (three tissues) 360 (1 tissue) 399 (1 tissue)
##I'd also like to compare within a tissue SR, ST and PV

my $SNPs = shift(@ARGV) or die; #All_SNPs_with_genes_and_category - 
my $indir = shift(@ARGV) or die;

my %keep = ();
my %klist = ();
#X       19961695        FBgn0031081     FBgn0031081     five_prime_UTR
open(S, "<$SNPs"); ##look through all the SNPs 
while(my $line = <S>){
    chomp $line;
    my @s = split(/\t/, $line);
    if($s[4] !~ m/intron/){ #can skip these
	if(!(exists($keep{$s[0] . "\t" . $s[1]}))){	    
	    $klist{$s[0] . "\t" . $s[1]} = $s[2];  ##chrom pos  = FBgn
	    $keep{$s[0] . "\t" . $s[1]} = 1;
	}elsif(exists($klist{$s[0] . "\t" . $s[1]})){
	    if($klist{$s[0] . "\t" . $s[1]} ne $s[2]){ #different gene - don't keep - I only want unique to a single gene SNPs
		$keep{$s[0] . "\t" . $s[1]}++;
	    }
	}
    }
}
close S;

opendir DIR, "$indir";
my @both = grep {/ALL_distinguishingSNPs$/} readdir DIR; ##has all distinguishing SNPs per cross - this is potential useful info
closedir DIR;

opendir DIR2, "$indir";
my @same = grep {/Same_masked$/} readdir DIR2; ##All SNPs masked that are the same between parents of a cross
closedir DIR2;

opendir DIR3, "$indir";
my @diff = grep {/Diff_masked$/} readdir DIR3;##All SNPs masked that are the different between parents of a cross
closedir DIR3;


my %masked = (); # I want to screen out any masked SNPs - by cross parents - so anything masked in 517 will always be masked, others may vary

foreach my $mask (@same){
    my @m = split(/_/, $mask);
    open(M, "<$mask");
    while(my $line = <M>){
	chomp $line;
	my @x = split(/\t/, $line);
	if($x[4] eq "Masked"){
	    $masked{$m[0] . "\t" . $x[2] . "\t" . $x[3]} = 1; ##p1  chrom pos
	}
    }
    close M;
}

foreach my $mask2 (@diff){
    my @m2 = split(/_/, $mask2);
    open(M2, "<$mask2");
    while(my $line = <M2>){
	chomp $line;
	my @x = split(/\t/, $line);
	if($x[4] eq "Masked"){
	    $masked{$m2[0] . "\t" . $x[2] . "\t" . $x[3]} = 1; ##p1  chrom pos
	}
    }
    close M2;
}

my @line = (304,307,360,399);
my @tissue = ("PV","SR","ST");

my %FsnpsLine = ();
my %MsnpsLine = ();

my %FsnpsTissue = ();
my %MsnpsTissue = ();

my %FsnpsBoth = ();
my %MsnpsBoth = ();
my %posLine = (); #stores all non-masked positions by line
my %posTissue = (); #stores all non-masked positions by tissue
my $header = "";

#304_SR_ALL_distinguishingSNPs
foreach my $both (@both){  #going through the distinguishing SNPs now
    my @x = split(/_/, $both);

    $both = $indir . $both;

    open(A, "<$both");
    while(my $line = <A>){
	chomp $line;
	my @a = split(/\s+/, $line);
#304     517     chr     pos     304_SNP 517_SNP vcfSNP  Ref_cov Alt_cov RefSNP  AltSNP  
#	print $a[11], "\t", $a[12], "\n";
	if($line =~ m/pos/){
	    $header = $a[2] . "\t" . $a[3] . "\t" . "FemaleCov" . "\t" . "MaleCov";
	}else{
	    if($keep{$a[2] . "\t" . $a[3]} == 1){ #only use SNPs unique to a gene
		if(!(exists($masked{$x[0] . "\t" . $a[2] . "\t" . $a[3]}))){ ##masking is by line

		    ##If I get to here then the data from the SNPs file is not masked - it can be added regardless of Line or Tissue
		    
		    $posLine{$x[0] . "\t" . $a[2] . "\t" . $a[3]} = 1; #what position belongs to what gene - for each cross - UPDATED 

		    $posTissue{$x[1] . "\t" . $a[2] . "\t" . $a[3]} = 1; #what position belongs to what gene - for each cross - UPDATED 
		    
		    my $pos = $a[2] . "\t" . $a[3];
		    #By line
		    ##to get the sex specific coverage
		    my $ftmp = 0;
		    my $mtmp = 0;	    
		    if($a[4] eq $a[9]){##female is reference
			$ftmp = $a[7];
		    }
		    if($a[4] eq $a[10]){#female is alt
			$ftmp = $a[8];
		    }
		    if($a[5] eq $a[9]){#male is ref
			$mtmp = $a[7];
		    }
		    if($a[5] eq $a[10]){#male is alt
			$mtmp = $a[8];
		    }
		    
		    if(!(exists($FsnpsLine{$x[0] . "\t" . $pos}))){
			$FsnpsLine{$x[0] . "\t" . $pos} = $ftmp;
		    }elsif(exists($FsnpsLine{$x[0] . "\t" . $pos})){
			$FsnpsLine{$x[0] . "\t" . $pos} = $FsnpsLine{$x[0] . "\t" . $pos} + $ftmp;
		    }
		    if(!(exists($MsnpsLine{$x[0] . "\t" . $pos}))){
			$MsnpsLine{$x[0] . "\t" . $pos} = $mtmp;
		    }elsif(exists($MsnpsLine{$x[0] . "\t" . $pos})){
			$MsnpsLine{$x[0] . "\t" . $pos} = $MsnpsLine{$x[0] . "\t" . $pos} + $mtmp;
		    }
		    
		    #By tissue
		    if(!(exists($FsnpsTissue{$x[1] . "\t" . $pos}))){
			$FsnpsTissue{$x[1] . "\t" . $pos} = $ftmp;
		    }elsif(exists($FsnpsTissue{$x[1] . "\t" . $pos})){
			$FsnpsTissue{$x[1] . "\t" . $pos} = $FsnpsTissue{$x[1] . "\t" . $pos} + $ftmp;
		    }
		    if(!(exists($MsnpsTissue{$x[1] . "\t" . $pos}))){
			$MsnpsTissue{$x[1] . "\t" . $pos} = $mtmp;
		    }elsif(exists($MsnpsTissue{$x[1] . "\t" . $pos})){
			$MsnpsTissue{$x[1] . "\t" . $pos} = $MsnpsTissue{$x[1] . "\t" . $pos} + $mtmp;
		    }
		    
		    
		    ##By both
		    if(!(exists($FsnpsBoth{$x[0] . "\t" . $x[1] . "\t" . $pos}))){
			$FsnpsBoth{$x[0] . "\t" . $x[1] . "\t" . $pos} = $ftmp;
		    }elsif(exists($FsnpsBoth{$x[0] . "\t" . $x[1] . "\t" . $pos})){
			$FsnpsBoth{$x[0] . "\t" . $x[1] . "\t" . $pos} = $FsnpsBoth{$x[0] . "\t" . $x[1] . "\t" . $pos} + $ftmp;
		    }
		    if(!(exists($MsnpsBoth{$x[0] . "\t" . $x[1] . "\t" . $pos}))){
			$MsnpsBoth{$x[0] . "\t" . $x[1] . "\t" . $pos} = $mtmp;
		    }elsif(exists($MsnpsBoth{$x[0] . "\t" . $x[1] . "\t" . $pos})){
			$MsnpsBoth{$x[0] . "\t" . $x[1] . "\t" . $pos} = $MsnpsBoth{$x[0] . "\t" . $x[1] . "\t" . $pos} + $mtmp;
		    }
		}
	    }
	}
    }
    close A;
}


foreach my $line (@line){
    my $output = $indir . $line . ".combined.counts";
    unlink(qq{$output});

    open(X, ">>$output");
    print X $header, "\n";
    while((my $k, my $v) = each(%posLine)){
	my @t = split(/\t/, $k);
	if($t[0] eq $line){	##only pick the ones for the right cross
	    if(!(exists($FsnpsLine{$k}))){
		$FsnpsLine{$k} = 0;
	    }
	    if(!(exists($MsnpsLine{$k}))){
		$MsnpsLine{$k} = 0;
	    }	
	    print X $t[1], "\t", $t[2], "\t", $FsnpsLine{$k}, "\t", $MsnpsLine{$k}, "\n";
	}
    }
    close X;
}

foreach my $tissue (@tissue){
    my $output2 = $indir . $tissue . ".combined.counts";
    unlink(qq{$output2});

    open(Y, ">>$output2");
    print Y $header, "\n";
    while((my $k, my $v) = each(%posTissue)){
	my @t = split(/\t/, $k);
	if($t[0] eq $tissue){	##only pick the ones for the right cross
	    if(!(exists($FsnpsTissue{$k}))){
		$FsnpsTissue{$k} = 0;
	    }
	    if(!(exists($MsnpsTissue{$k}))){
		$MsnpsTissue{$k} = 0;
	    }
	    
	    print Y $t[1], "\t", $t[2], "\t", $FsnpsTissue{$k}, "\t", $MsnpsTissue{$k}, "\n";
	}
    }
    close Y;
}

foreach my $l (@line){
    foreach my $t (@tissue){
	my $snpfile = $l . "_" . $t . "_ALL_distinguishingSNPs";

	if(-e($snpfile)){
	    my $output3 = $indir . $l . "_" . $t . ".combined.counts";
	    unlink(qq{$output3});
	    
	    open(Z, ">>$output3");
	    print Z $header, "\n";
	    while((my $k, my $v) = each(%posLine)){
		my @t = split(/\t/, $k);
		if($t[0] eq $l){	##only pick the ones for the right cross
		    if(!(exists($FsnpsBoth{$l . "\t" . $t . "\t" . $t[1] . "\t" . $t[2]}))){
			$FsnpsBoth{$l . "\t" . $t . "\t" . $t[1] . "\t" . $t[2]} = 0;
		    }
		    if(!(exists($MsnpsBoth{$l . "\t" . $t . "\t" . $t[1] . "\t" . $t[2]}))){
			$MsnpsBoth{$l . "\t" . $t . "\t" . $t[1] . "\t" . $t[2]} = 0;
		    }
		    print Z $t[1], "\t", $t[2], "\t", $FsnpsBoth{$l . "\t" . $t . "\t" . $t[1] . "\t" . $t[2]}, "\t", $MsnpsBoth{$l . "\t" . $t . "\t" . $t[1] . "\t" . $t[2]}, "\n";
		}
	    }
	    close Z;
	}
    }
}
