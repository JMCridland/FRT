#!perl -w

use strict;

my $indir = shift(@ARGV) or die;#/data/julie/FemaleRT/
my $tag = shift(@ARGV) or die;##using the transcript_coverage files 
my $snpmin = shift(@ARGV) or die;##the number of SNPs
my $covmin = shift(@ARGV) or die;##the minimum total coverage
my $tpercent = shift(@ARGV) or die; ##moving to 75%

###I am doing this in a by gene way using the longest transcript (CDS)
#my $percentmaletrans = 1; ##percentage of the observed female transcript length required for the male transcript
#my $percentmaleSNPs = 0.9; #the percentage I require of male SNPs for total female SNPs
#my $output = shift(@ARGV) or die;
#unlink(qq{$output});

opendir DIR, "$indir";
my @trans = grep{/\.$tag/} readdir DIR;
closedir DIR;
   
foreach my $transcript (@trans){ #these are the files with the male and female coverage per gene
    ###I want to impose a minimum number of SNPs per transcript (male or female)
    ###I also need a minimum transcript coverage - 10? over SNPs
    ###to consider a male transcript I want it to have a minimum percentage of the SNPs I see in the female transcript
    my @f = split(/\./, $transcript);
    
    if($f[0] =~ m/_/){ ##I only want to look at the ones that are RAL and Tissue i.e. 304_SR	
	my $transcript = $indir . $f[0] . ".transcript_coverage";
	my $output = $indir . $f[0] . "." . $snpmin  . ".SNPs." . $covmin . ".cov.high_confidence.transcripts";
	unlink(qq{$output});
	open(B, ">>$output");
	
	open(T, "<$transcript"); ##this is the file that looks at transcript coverage
	while(my $line = <T>){
	    chomp $line;
	    my @t = split(/\t/, $line);
	    #GeneID  FemaleNum       F_first F_last  Fpercent        FCov    MaleNum M_first M_last  Mpercent        MCov    CDSdist	    
	    if($line =~ m/GeneID/){
		print B $line, "\t", "Tag\n";
	    }elsif($line !~ m/GeneID/){
		my $ftrans = $t[4]; ##percent coverage
		my $mtrans = $t[9]; ##percent coverage
		my $dist = $t[11]; ##distance
	        my $fnumSNPs = $t[1]; ##number of female SNPs used in transcript calculation
		my $mnumSNPs = $t[6]; #number of male SNPs used in transcript calculaton
		my $fcov = $t[5]; ##Female Coverage
		my $mcov = $t[10];##Male coverage

		my $male = 0;
		my $female = 0;
		my $keep = 0;
		if(($fcov >= $covmin) and ($fnumSNPs >= $snpmin) and ($ftrans >= $tpercent)){ ##start with 10 and 3
		    $female = 1;
		    $keep = 1;
		 #   if(($mnumSNPs >= ($fnumSNPs * $percentmaleSNPs)) and ($mcov >= $covmin) and ($mtrans >= $tpercent) and ($mtrans >= ($ftrans * $percentmaletrans))){##is it also expressed in the male?
		#    if(($mnumSNPs >= $snpmin) and ($mcov >= $covmin) and ($mtrans >= $tpercent) and ($mtrans >= ($ftrans * $percentmale))){
			###currently requiring at least some percentage of the SNPs in female and the same transcript coverage
		#	$male = 1;
		  #  }	    
		}#else{##Did not meet female transcript minimum - is it in the male?
		if(($mcov >= $covmin) and ($mnumSNPs >= $snpmin) and ($mtrans >= $tpercent)){##is it also expressed in the male?
		    $male = 1;
		    $keep = 1;
		}
		#}
		###Now determine if I am calling things Female, Male, Both or Neither
		my $tag = "";
		if(($female == 1) and ($male == 1)){
		    $tag = "Both";
		}elsif(($female == 1) and ($male == 0)){
		    $tag = "Female";
		}elsif(($female == 0) and ($male == 1)){
		    $tag = "Male";
		}elsif(($female == 0) and ($male == 0)){
		    $tag = "Neither";
		}

		#if($keep == 1){
		    print B $line, "\t", $tag, "\n";
		#}
	    }
	}
	close T;
	close B;
    }
}
    
