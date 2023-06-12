#!perl -w

use strict;

my $transcript = shift(@ARGV) or die; #/data/FlyRef/Drosophila_melanogaster/6.41/fasta/dmel-all-exon...
my $outdir = shift(@ARGV) or die;#/data/julie/Dmel_CDS_lookup/
###For each gene I need to make a lookup table with the position in the CDS and the position of each snp
####Pos    ChPos
####1      400000
####2      400001



##Nep3-PA type=CDS; loc=X:join(19963955..19964071,19964782..19964944,19965006..19965126,19965197..19965511,19965577..19966071,19966183..19967012,19967081..19967223,19967284..19967460); name=Nep3-RA; dbxref=FlyBase:FBpp0070000,FlyBase_Annotation_IDs:CG9565-PA,REFSEQ:NP_523417,GB_protein:AAF45370,UniProt/Swiss-Prot:Q9W5Y0,FlyMine:FBpp0070000,modMine:FBpp0070000; MD5=5c3c0b466b23e32d99dfff926a5e8c6b; length=2361; parent=FBgn0031081,FBtr0070000; release=r6.41; species=Dmel; 

my %pID = (); ##FBgn
my %length = (); #longest CDS
my %loc = (); #location string

open(T, "gunzip -c $transcript | "); ##First go through the CDS file and find the longest 
while(my $line = <T>){
    chomp $line;
    if($line =~ m/^>/){
	$line =~ s/\s+//go;
	my @t = split(/;/, $line);

	my $length = "";
	my $loc = "";
	my $p = "";

	foreach my $t (@t){
	    if($t =~ m/parent=/){
		$p = $t;
	    }
	    if($t =~ m/length=/){
		$length = $t;
	    }
	    if($t =~ m/loc=/){
		$loc = $t;
	    }	    
	}

	$p =~ s/parent=//;
	my @p = split(/,/, $p);
	my $parent = $p[0];

	$length =~ s/length=//;

	$loc =~ s/loc=//go;
	$loc =~ s/join//go;
	$loc =~ s/complement//go;
	$loc =~ s/\(//go;
	$loc =~ s/\)//go;
	$loc =~ s/\.\./\./go;
	
	if(!(exists($pID{$parent}))){
	    $pID{$parent} = $parent;
	    $length{$parent} = $length;
	    $loc{$parent} = $loc;
	}elsif(exists($pID{$parent})){
	    if($length > $length{$parent}){
		$pID{$parent} = $parent;
		$length{$parent} = $length;
		$loc{$parent} = $loc;
	    }
	}
    }
}
   
    
###Now make a table with the positions for the longest CDS

while((my $k, my $v) = each(%pID)){
  #  print $k, "\n";
    ##I need to sort the location string and print out a table for each gene
    my $output = $outdir . $k . ".transcript.lookup";
    unlink(qq{$output});

    my @l = split(/:/, $loc{$k});
    my $ch = $l[0];
    my @l2 = split(/,/, $l[1]);
    my @l3 = sort @l2;
    
    open(X, ">>$output");
    #print X $k, "\t", $length{$k}, "\t", $loc{$k}, "\n";

    my $start = 1;
    foreach my $x (@l3){
	my @y = split(/\./, $x);
	for(my $i = $y[0]; $i <= $y[1]; $i++){
	    print X $k, "\t", $length{$k}, "\t", $ch, "\t", $start, "\t", $i, "\n";
	    $start++;
	}
    }
    close X;
}
