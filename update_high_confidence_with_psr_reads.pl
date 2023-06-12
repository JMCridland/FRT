#!perl -w

###The idea here is to update my calls for some of the transcripts based on information from the transcript in other tissues in the same cross and the parent specific read information
use strict;

my $indir = shift(@ARGV) or die; #/data/julie/FemaleRT/
my $indir2 = shift(@ARGV) or die; #/data/julie/FemaleRT/PSR/

opendir DIR, "$indir";
my @high = grep{/\.high_confidence\.transcripts$/} readdir DIR;
closedir DIR;

my @lines = (304,307,360);
my %tissue = ();
@{$tissue{304}} = ("PV","SR","ST");
@{$tissue{307}} = ("PV","SR","ST");
@{$tissue{360}} = ("SR","ST");

my $header = "";
my %outhigh = ();
my %male_gene = (); ##is it male?
my %female_gene = (); #is it female
my %male_data = ();
my %female_data = ();
my %genes = (); #for each output
my %table = (); #store the file records to print back out
foreach my $high (@high){

    my @h = split(/\./, $high);
    #304_PV.5.SNPs.10.cov.high_confidence.transcripts
    
    my @i = split(/_/, $h[0]);
    my $out = $high . ".updated";
    $outhigh{$i[0] . "\t" . $i[1]} = $out;
    
    $high = $indir . $high;
    open(H, "<$high");
    while(my $line = <H>){
	chomp $line;
	my @s = split(/\t/, $line);
	if($line =~ m/GeneID/){
	    $header = $line;
	}else{
	    $table{$s[0] . "\t" . $i[0] . "\t" . $i[1]} = $line;
	    push(@{$genes{$i[0] . "\t" . $i[1]}}, $s[0]);
	    ##FBgn0024558    5       22582237        22583211        0.802   410     5       22582237        22583211        0.802   321     1128    Both
	    if($s[12] =~ m/Both/){
		$male_gene{$s[0] . "\t" . $i[0]}++; #already called as male - key is gene line
		$female_gene{$s[0] . "\t" . $i[0]}++; #already called as male - key is gene line
		
	    }elsif($s[12] =~ m/Male/){  
		$male_gene{$s[0] . "\t" . $i[0]}++; #already called as male - key is gene line
		
	    }elsif($s[12] =~ m/Female/){
		$female_gene{$s[0] . "\t" . $i[0]}++; #already called as female - key is gene line		    
	    }
	    ##key is gene line tissue
	    $male_data{$s[0] . "\t" . $i[0] . "\t" . $i[1]} = $s[6] . "\t" . $s[10] . "\t" . $s[11]; ###Number of SNPS, SNP coverage, Transcript length	
	    $female_data{$s[0] . "\t" . $i[0] . "\t" . $i[1]} = $s[1] . "\t" . $s[5] . "\t" . $s[11]; ###Number of SNPS, SNP coverage, Transcript length
	}
    }
    close H;
}

opendir DIR2, "$indir2";
my @psr = grep{/\.ps\.read\.counts/} readdir DIR2;
closedir DIR2;

my %male_psr_trans = ();
my %male_psr_count = ();
my %female_psr_trans = ();
my %female_psr_count = ();
foreach my $psr (@psr){
    my @pr = split(/x/, $psr);
    my @pr2 = split(/_/, $pr[0]);
    
    $psr = $indir2 . $psr;
    open(P, "<$psr");
    ##GeneID  Start1  Stop1   Start2  Stop2   Total_reads1    Total_reads2
    while(my $line = <P>){
	chomp $line;
	my @p = split(/\t/, $line);
	##key is gene line tissue
	if($line !~ m/GeneID/){
	    if($p[3] != 9999999999999){
		$male_psr_trans{$p[0] . "\t" . $pr2[0] . "\t" . $pr2[1]} = abs($p[4] - $p[3])+1;
		$male_psr_count{$p[0] . "\t" . $pr2[0] . "\t" . $pr2[1]} = $p[6];
	    }
	     if($p[1] != 9999999999999){
		$female_psr_trans{$p[0] . "\t" . $pr2[0] . "\t" . $pr2[1]} = abs($p[2] - $p[1])+1;
		$female_psr_count{$p[0] . "\t" . $pr2[0] . "\t" . $pr2[1]} = $p[5];
	    }
	}
    }
    close P;
}

my %male_update = (); #hold the updated records
my %female_update = (); #hold the updated records
foreach my $l (@lines){
    while((my $k, my $v) = each(%male_data)){
	my @k = split(/\t/, $k);###key is gene line tissue
	if($l eq $k[1]){ ###on the right RAL line
	    if(exists($male_gene{$k[0] . "\t" . $k[1]})){ #could be updated - means that at least one is already identified as male
		if($male_gene{$k[0] . "\t" . $k[1]} < scalar(@{$tissue{$l}})){
		    #they are all already male - no need to check
		    my @y = split(/\t/, $v);
		    my $newsnp = $y[0];
		    my $newcov = $y[1];		    
		    #   print $k[0], "\t", $newsnp, "\t", $newcov, "\n";
		    if(exists($male_psr_trans{$k})){
			my $newtrans = sprintf("%.3f",($male_psr_trans{$k} / $y[2]));
			#	print $k[0], "\t", $newtrans, "\n";
			if(($newsnp >= 5) and ($newcov >= 10) and ($newtrans >= 0.75)){
			    #    print $k[0], "\t", "Update\n";
			    $male_update{$k} = $newtrans . "\t" . $newcov . "\t" . $newsnp; ##to update!!!
			}			
		    }
		}
	    }
	}
    }
}
foreach my $l (@lines){
    while((my $k, my $v) = each(%female_data)){
	my @k = split(/\t/, $k);###key is gene line tissue
	if(exists($female_gene{$k[0] . "\t" . $k[1]})){ #could be updated - means that at least one is already identified as male
	    if($female_gene{$k[0] . "\t" . $k[1]} < scalar(@{$tissue{$l}})){
		#they are all already male - no need to check
		my @y = split(/\t/, $v);
		my $newsnp = $y[0];
		my $newcov = $y[1];
		#   print $k[0], "\t", $newsnp, "\t", $newcov, "\n";
		if(exists($female_psr_trans{$k})){
		    my $newtrans = sprintf("%.3f", ($female_psr_trans{$k} / $y[2]));
		    #	print $k[0], "\t", $newtrans, "\n";
		    if(($newsnp >= 5) and ($newcov >= 10) and ($newtrans >= 0.75)){
			#    print $k[0], "\t", "Update\n";
			$female_update{$k} = $newtrans . "\t" . $newcov . "\t" . $newsnp; ##to update!!!
		    }			
		}
	    }
	}
    }
}	

while((my $k2, my $v2) = each(%outhigh)){
    unlink(qq{$v2});

    open(X, ">>$v2");
    print X $header, "\n";
    foreach my $gene (@{$genes{$k2}}){
	my @t = split(/\t/, $table{$gene . "\t" . $k2});
	##GeneID  FSNPs   F_first F_last  FpercentTrans   FCov    MSNPs   M_first M_last  MpercentTrans   MCov    CDSdist Tag

	my $tag = $t[12];
	my $male = 0;
	my $female = 0;

	if($tag eq "Both"){
	    $male++;
	    $female++;
	}
	if($tag eq "Male"){
	    $male++;
	}
	if($tag eq "Female"){
	    $female++;
	}

	my $fsnp = $t[1];
	my $ftrans = $t[4];
	my $fcov = $t[5];
	my $msnp = $t[6];
	my $mtrans = $t[9];
	my $mcov = $t[10];
		
	
	if(exists($male_update{$gene . "\t" . $k2})){
	    my @m = split(/\t/, $male_update{$gene . "\t" . $k2});
	    $mtrans = $m[0];
	    $mcov = $m[1];
	    $msnp = $m[2];
	#    print "Male\t", $gene, "\t", $k2, "\n";
	    $male++;	    
	}
	if(exists($female_update{$gene . "\t" . $k2})){
	    my @f = split(/\t/, $female_update{$gene . "\t" . $k2});
	    $ftrans = $f[0];
	    $fcov = $f[1];
	    $fsnp = $f[2];
	  #  print "Female\t", $gene, "\t", $k2, "\n";
	    $female++;	    
	}
	
	if(($female > 0) and ($male > 0)){
	    $tag = "Both";
	}elsif(($female > 0) and ($male == 0)){
	    $tag = "Female";
	}elsif(($female == 0) and ($male > 0)){
	    $tag = "Male";
	}else{
	 $tag = "Neither";   
	}
	##GeneID  FSNPs   F_first F_last  FpercentTrans   FCov    MSNPs   M_first M_last  MpercentTrans   MCov    CDSdist Tag
	print X $t[0], "\t", $fsnp, "\t", $t[2], "\t", $t[3], "\t", $ftrans, "\t", $fcov, "\t", $msnp, "\t", $t[7], "\t", $t[8], "\t", $mtrans, "\t", $mcov, "\t", $t[11], "\t", $tag, "\n";
    }
    close X;
}

