#!/usr/bin/perl -w
use strict;
my $data_chunk_size = 6;
my @names = qw(Nc Nt Nd);
my $file = shift || die $!;
open(my $fh =>$file) || die "$file: $!";
my $header =<$fh>;
my @header = split(/\s+/,$header);
my $i = 0;
my %lookup = map { $_=> $i++ } @header;
my %stats;

while(<$fh> ){
    chomp;
    my @row = split(/\t/,$_);
    my $single_match = pop @row;

    my $count = 0;
    my %indata;
    while( @row ) {
	if( $row[0] ne '' ) {
	    my (undef,$name) = split(/_/,$header[$count*$data_chunk_size]);
	    $indata{$name} = [splice(@row,0,$data_chunk_size)];
	} else {
	    splice(@row,0,$data_chunk_size);
	}
	$count++;
    }
    $stats{join('-',sort keys %indata)}++;
    
}

for my $stat ( sort { $stats{$b} <=> $stats{$a} } keys %stats ) {
    print join("\t", $stat, $stats{$stat}),"\n";
}
