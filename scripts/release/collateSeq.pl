#!/usr/local/bin/perl -w

use strict;

my $fafile = shift;
my @seqs;
open(OUT, ">$fafile.noprot") || die "cant open the out file\n";
open(FA, "<$fafile") || die "cant open the fasta file\n";
@seqs=<FA>;
close(FA);

foreach (my $i=0; $i<=scalar(@seqs)-1; $i++){
    if ($seqs[$i]=~/protein/){++$i; next;}
    if ($seqs[$i]=~/mol\:na/){
        print OUT $seqs[$i], $seqs[$i+1];
        ++$i;
    }
}
