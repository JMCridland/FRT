#!perl -w

use strict;

my $indir = shift(@ARGV) or die;#path to the transcript files
my $tag = shift(@ARGV) or die;##the label at the end of the files to grep for
my $snpmin = shift(@ARGV) or die;##the number of SNPs
my $covmin = shift(@ARGV) or die;##the minimum total coverage
my $tpercent = shift(@ARGV) or die; ##the percent coverage of the transcript

opendir DIR, "$indir";
my @trans = grep{/\.$tag/} readdir DIR;
closedir DIR;
   
foreach my $transcript (@trans){ #these are the files with the male and female coverage per gene
    ###I want to impose a minimum number of SNPs per transcript (male or female)
    ###I also need a minimum transcript coverag
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
		}
		if(($mcov >= $covmin) and ($mnumSNPs >= $snpmin) and ($mtrans >= $tpercent)){##is it also expressed in the male?
		    $male = 1;
		    $keep = 1;
		}
		
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
		print B $line, "\t", $tag, "\n";
	    }
	}
	close T;
	close B;
    }
}
    
