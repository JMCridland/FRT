#!perl -w

use strict;

my $transcript = shift(@ARGV) or die; #dmel-all-transcript-r6.41.fasta - downloaded from FlyBase
my $outdir = shift(@ARGV) or die;#path to output directory
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
