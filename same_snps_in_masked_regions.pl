#!perl -w

use strict;

my $input = shift(@ARGV) or die; #Lack_tableS3.txt
my $same = shift(@ARGV) or die;
my $output = shift(@ARGV) or die;
unlink(qq{$output});

my @s = split(/_/, $same);
my $id = "RAL-" . $s[0];
my %lines = ();

#$lines{"RAL-304"} = 304;
#$lines{"RAL-307"} = 307;
#$lines{"RAL-357"} = 357;
#$lines{"RAL-360"} = 360;
#$lines{"RAL-399"} = 399;
$lines{$id} = $s[0];
$lines{"RAL-517"} = 517;

my %remove = (); #residual heterozygosity
open(A, "<$input");
#this is that weird file that is actually three tables, but poorly formatted so I really only care about the first four columns
#Line	Chrom	Start	Stop			Line	Chrom	Start	Stop			Line	True Heterozygosity	Masked
while(my $line = <A>){
    chomp $line;
    my @a = split(/\s+/, $line);
    #   print $a[0], "\n";
    $a[1] =~ s/Chr//;
    if(exists($lines{$a[0]})){
#	print $line, "\n";
	for(my $i = $a[2]; $i <= $a[3]; $i++){
	    $remove{$lines{$a[0]} . "\t" . $a[1] . "\t" . $i}= 1; #line chrom pos to remove
	}
    }
}
close A;

open(B, "<$same"); #304_ST_ALL_Same_snps
open(C, ">>$output");
while(my $line = <B>){
    chomp $line;
    my @a = split(/\t/, $line);
    ##304     517     X       467515  T       T       T,.     219     3       T       .

    if(exists($remove{$a[0] . "\t" . $a[2] . "\t" . $a[3]}) or exists($remove{$a[1] . "\t" . $a[2] . "\t" . $a[3]})){
	print C $a[0], "\t",  $a[1], "\t", $a[2], "\t", $a[3], "\t", "Masked", "\n";
    }else{
	print C $a[0], "\t",  $a[1], "\t", $a[2], "\t", $a[3], "\t", "Notmasked", "\n";
    }
}
close B;
close C;


