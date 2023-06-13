#!perl -w

use strict;

my $indir = shift(@ARGV) or die;#path to the directory with the gene.info files output from add_label_to_combined_counts.pl
my $tag = shift(@ARGV) or die;# *gene.info  
my $cov = shift(@ARGV) or die; # the minimum coverage
my $cdsdir = shift(@ARGV) or die; #/directory with the output files from get_intron_exon_longest_transcript.pl

opendir DIR2, "$cdsdir";
my @lookup = grep {/\.lookup/} readdir DIR2;
closedir DIR2;

my %position = ();
my %trans = ();
foreach my $lookup (@lookup){ ##I have to save all these positions
#    print $lookup, "\n";
    $lookup = $cdsdir . $lookup;
    open(X, "<$lookup");
    while(my $Cline = <X>){
	chomp $Cline;
	my @x = split(/\t/, $Cline);		
	##FBgn0041626     1164    X       0       20141819
	$position{$x[0] . "\t" . $x[4]} = $x[3];
	$trans{$x[0]} = $x[1];
    }
}
close X;
print "Lookup stored\n";
	    
opendir DIR, "$indir";
my @counts = grep{/\.$tag/} readdir DIR;
closedir DIR;

foreach my $counts (@counts){
    my @f = split(/\./, $counts);
    if($f[0] =~ m/_/){ ##I only want to look at the ones that are RAL and Tissue i.e. 304_SR
	my $counts = $indir . $counts;
	
	my %male = ();
	my %female = ();
	my %mfirst = ();
	my %mlast = ();
	my %ffirst = ();
	my %flast = ();
	my %keep = ();
	my %fcov = ();
	my %mcov = ();
	open(A, "<$counts");###First store the gene info
	while(my $line = <A>){
	    ##chr     pos     FemaleCov       MaleCov GeneInfo
	    #2R      13838828        1       0       FBgn0033874:exon,CDS,gene
	    if($line !~ m/GeneInfo/){
		my @a = split(/\t/, $line);
		my @b = split(/:/, $a[4]);

		if(($line =~ m/exon/) or ($line =~ m/UTR/)){

		    if(exists($position{$b[0] . "\t" . $a[1]})){ #position is part of the longest transcript
			if($a[2] >= $cov){
			    $keep{$b[0]} = 1;
			    if(!(exists($female{$b[0]}))){
				$fcov{$b[0]} = $a[2];
				$female{$b[0]} = 1;
				$ffirst{$b[0]} = $a[1];
				$flast{$b[0]} = $a[1];
			    }elsif(exists($female{$b[0]})){ ##update
				$female{$b[0]}++;
				$fcov{$b[0]} = $fcov{$b[0]} + $a[2];
				if($a[1] < $ffirst{$b[0]}){
				    $ffirst{$b[0]} = $a[1];
				}
				if($a[1] > $flast{$b[0]}){
				    $flast{$b[0]} = $a[1];
				}
			    }
			}
			###male
			if($a[3] >= $cov){
			    $keep{$b[0]} = 1;
			    if(!(exists($male{$b[0]}))){
				$male{$b[0]} = 1;
				$mcov{$b[0]} = $a[3];
				$mfirst{$b[0]} = $a[1];
				$mlast{$b[0]} = $a[1];
			    }elsif(exists($male{$b[0]})){ ##update
				$male{$b[0]}++;
				$mcov{$b[0]} = $mcov{$b[0]} + $a[3];
				if($a[1] < $mfirst{$b[0]}){
				    $mfirst{$b[0]} = $a[1];
				}
				if($a[1] > $mlast{$b[0]}){
				    $mlast{$b[0]} = $a[1];
				}
			    }
			}
		    }
		}
	    }
	} ###done storing
	
	my $output = $f[0] . ".transcript_coverage";
	unlink(qq{$output});
	open(B, ">>$output");
	print B "GeneID\tFSNPs\tF_first\tF_last\tFpercentTrans\tFCov\tMSNPs\tM_first\tM_last\tMpercentTrans\tMCov\tCDSdist\n";

	while((my $k, my $v) = each (%keep)){ ###This is now by gene - use the lookup table to calculate real transcript length
	    
	    if(!(exists($male{$k}))){
		$male{$k} = 0;
		$mfirst{$k} = "NA";
		$mlast{$k} = "NA";
		$mcov{$k} = 0;
	    }
	    if(!(exists($female{$k}))){
		$female{$k} = 0;
		$ffirst{$k} = "NA";
		$flast{$k} = "NA";
		$fcov{$k} = 0;
	    }
	    my $mlen = 0;
	    my $flen = 0;
	    
	    if($flast{$k} !~ m/NA/){
		$flen = abs($position{$k . "\t" . $flast{$k}} - $position{$k . "\t" . $ffirst{$k}}) + 1;		
	    }
	    if($mlast{$k} !~ m/NA/){	  
		$mlen = abs($position{$k . "\t" . $mlast{$k}} - $position{$k . "\t" . $mfirst{$k}}) + 1;
	    }
	    my $fcov = sprintf("%.3f", ($flen / $trans{$k}));
	    my $mcov = sprintf("%.3f", ($mlen / $trans{$k}));
	    
	    print B $k, "\t", $female{$k}, "\t", $ffirst{$k}, "\t", $flast{$k}, "\t", $fcov, "\t", $fcov{$k}, "\t", $male{$k}, "\t", $mfirst{$k}, "\t", $mlast{$k}, "\t", $mcov, "\t", $mcov{$k}, "\t", $trans{$k}, "\n";
    }
	close B;
	%keep = ();	
	%male = ();
	%female = ();
	%mfirst = ();
	%mlast = ();
	%ffirst = ();
	%flast = ();
	%fcov = ();
	%mcov = ();
    }
}
