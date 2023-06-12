#!perl -w

use strict;

my $transcript = shift(@ARGV) or die;#/data/FlyRef/Drosophila_melanogaster/6.41/fasta/dmel-transcript-reformatted.fasta.gz
my $ncRNA = shift(@ARGV) or die; #/data/FlyRef/Drosophila_melanogaster/6.41/fasta/dmel-ncRNA-reformatted.fasta.gz
my $snps = shift(@ARGV) or die; #ALL_SNPs_table
my $lookupdir = shift(@ARGV) or die; #the directory with the lookup tables

#my $out = shift(@ARGV) or die;
#unlink(qq{$out});

###I am going to need the transcript fastas - I want the longest transcript for each gene
my %pID = (); ##FBgn
my %fasta = ();##store longest fasta
my %length = ();
my %strand = (); #+ or -
my $id = "";

my @tr = ($transcript,$ncRNA);
foreach my $tr (@tr){
    open(T, "gunzip -c $tr | "); ##First go through the CDS file and find the longest
    ###The minus strand transcripts must be reverse complemented!!
  LOOP:while(my $line = <T>){
      chomp $line;
      if($line =~ m/^>/){
	  $line =~ s/\s+//go;
	  my @t = split(/;/, $line);
	  my $length = "";
	  my $loc = "";
	  my $p = "";
	  
	  my $strand = "+";
	  if($line =~ m/complement/){
	      $strand = "-";
	  }
	  foreach my $t (@t){
	      if($t =~ m/parent=/){
		  $p = $t;
	      }     
	  }
	  $p =~ s/parent=//;
	  my @p = split(/,/, $p);
	  my $parent = $p[0];
	  $id = $p[0]; ##this is the current ID
	  
	  if(!(exists($pID{$parent}))){
	      $pID{$parent} = $parent;
	      $strand{$parent} = $strand;
	  }
      }else{
	  if(!(exists($fasta{$id}))){
	      if($strand{$id} eq "-"){
		  $line = reverse $line;
		  $line =~ tr/ACGT/TGCA/;
	      }	    
	      $fasta{$id} = $line;
	      $length{$id} = length($line);
	  }elsif(exists($fasta{$id})){
	      if($strand{$id} eq "-"){
		  $line = reverse $line;
		  $line =~ tr/ACGT/TGCA/;
	      }
	      if(length($line) > length($fasta{$id})){
		  $fasta{$id} = $line;
		  $length{$id} = length($line);
	      }
	  }
      }
  }
    close T;
    print "Transcripts found\n";
}

my %variant = (); #store the SNPs I need to update - only keep distinguishing SNPs 
my %conv = ();
my %usedpos = (); ##only store the positions I care about

open(S, "<$snps");
while(my $line = <S>){
    #chr     pos     ref     alt     line_304        line_307        line_357        line_360        line_399        line_517        GeneID  GeneName
    chomp $line;
    my @s = split(/\t/, $line);
    if($line =~ m/pos/){
	for(my $i = 4; $i <= 9; $i++){
	    $s[$i] =~ s/line_//;
	    $conv{$i} = $s[$i];
	}
    }elsif($line !~ m/pos/){
	if((length($s[2]) == 1) and (length($s[3]) == 1)){ #I don't want to deal with indels
	    for(my $j = 4; $j <= 8; $j++){
		if(($s[$j] !~ m/-/) and ($s[9] !~ m/-/)){		    
		    if($s[$j] ne $s[9]){
			### Line chrom pos gene = SNP
			$usedpos{$s[10] . "\t" . $s[0] . "\t" . $s[1]} = 1; ##chrom pos gene
			if($conv{$j} !~ m/357/){			    
			    $variant{$conv{$j} . "\t" . $s[0] . "\t" . $s[1] . "\t" . $s[10]} = $s[$j];
			    $variant{$conv{9} . "\t" . $s[0] . "\t" . $s[1] . "\t" . $s[10]} = $s[9];	
			}
		    }
		}
	    }
	}
    }
}
close S;
print "SNPs stored\n";
my @lines = (304,307,360,399,517);
my %trans_out = (); #stores the finished transcripts - I need a full set for each line
while((my $k, my $v) = each(%fasta)){ #initializing these with a reference copy of the longest transcript
    foreach my $l (@lines){
	$trans_out{$k . "\t" . $l} = $v; ### GeneID   Line
    }
}
print "Line specific transcritps stored\n";
###Now that I have all the SNPs stored and have the transcripts stored I need to update them using the *transcript.lookup
opendir DIR, "$lookupdir";
my @look = grep{/\.lookup$/} readdir DIR;
closedir DIR;

my %genepos = ();
foreach my $file (@look){##this stores the lookup - only need positions that are to be updated
    $file = $lookupdir . $file;
    open(X, "<$file");
    while(my $line = <X>){
	###I need chrom pos gene to within transcript position
	#FBgn0038272     1510    3R      1       14892392 - I had to update these to 1 base start positions!!!!!
	my @x = split(/\s+/, $line);
	if(exists($usedpos{$x[0] . "\t" . $x[2] . "\t" . $x[4]})){
	    $genepos{$x[0] . "\t" . $x[2] . "\t" . $x[4]} = $x[3]; #to keep track of which transcript position to update
	}
    }
    close X;
}
print "Lookup positions stored\n";
#open(Z, ">>$out");
##The transcript sequence is always given in the plus strand order - so I just need to update by position
while((my $k, my $v) = each(%variant)){
    # Line chrom pos gene = SNP
    my @var = split(/\t/, $k);
    if(exists($trans_out{$var[3] . "\t" . $var[0]})){
	if(exists($genepos{$var[3] . "\t" . $var[1] . "\t" . $var[2]})){
#	    print Z $var[3], "\t",  $k, "\t", $v, "\n";
	    
#	    print Z  $var[3], "\t", $genepos{$var[3] . "\t" . $var[1] . "\t" . $var[2]}, "\n";
	    ##Get the right transcript
#	    print Z  $var[3], "\t", $trans_out{$var[3] . "\t" . $var[0]}, "\n";
	    
	    ##Expr   Offset   Length   Replacement
	    substr($trans_out{$var[3] . "\t" . $var[0]}, ($genepos{$var[3] . "\t" . $var[1] . "\t" . $var[2]} - 1), 1, $v);
	  #  substr($trans_out{$var[3] . "\t" . $var[0]}, 5, 1, "X");
	    
#	    print Z $var[3], "\t", $trans_out{$var[3] . "\t" . $var[0]}, "\n";    
	    
	    #print X $k, "\t", $v, "\n";
	}
    }
}
#close Z;

print "Transcripts updated\n";

foreach my $id (@lines){
    my $output = $id . "specific.transcripts";
    unlink(qq{$output});

    open(T, ">>$output");
    while((my $k, my $v) = each(%pID)){
	
	print T ">", $k, "_", $id, "\n", $trans_out{$k . "\t" . $id}, "\n";
    }
    close T;
}



