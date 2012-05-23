#!/usr/bin/perl -w
use strict;

my $oldfile = 'data/old_K27_list/K27_ncra_10percent.tab';
my $newRSEG = 'data/RSEGdomain_genes/N_crassa/Nc_H3K27me3_genes.txt';

open(my $fh => $oldfile) || die $!;
my %old;
while(<$fh>) {
    my ($ncu) = split;
    $old{$ncu}++;
}

open($fh=>$newRSEG) || die $!;
my $shared;
my $missing;
while(<$fh>) {
    my ($ncu) = split;
    if( $old{$ncu} ) {
	$shared++;
    } else {
	$missing++;
    }
}

print "shared=$shared missing=$missing\n";
