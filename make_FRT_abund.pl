#!perl -w

use strict;

my $indir = shift(@ARGV) or die;
my $output = shift(@ARGV) or die;
unlink(qq{$output});


opendir DIR, "$indir";
my @abund = grep {/\.abund\.tab$/} readdir DIR;
closedir DIR;

my %data = ();
my %genes = ();
my @ids = ();
foreach my $file (@abund){
    my @s = split(/\./, $file);
    my $id = $s[0] . "_" . $s[1];
    push(@ids, $id);

    $file = $indir . $file;
    open(A, "<$file");
    while(my $line = <A>){
        chomp $line;
        if($line !~ m/TPM/){
            my @a = split(/\t/, $line);
            $genes{$a[0]} = 1;
            $data{$id . "\t" . $a[0]} = $a[8]; ##Sex&tissue info \t Gene  =  TPM
        }
    }
    close A;
}

my @sorted = sort (@ids);

open(B, ">>$output");
print B "GeneID";
foreach my $id (@sorted){
    print B "\t", $id;
}
print B "\n";
while((my $k, my $v) = each(%genes)){
    print B $k;
    foreach my $id (@sorted){
	print B "\t", $data{$id . "\t" . $k};
    }
    print B "\n";
}
close B;
