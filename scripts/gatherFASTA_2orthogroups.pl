#!/usr/bin/perl -w
use strict;
my $dir = shift;
my $db = shift || 'cleanseq';
opendir(DIR, $dir);
my %table;
my %sps;
for my $file ( readdir(DIR) ) {
    if( $file =~ /^(OG\d+\_\d+)-/ ) {
	my $og = $1;
	open(my $fh => "$dir/$file") || die $!;
	while(<$fh>) {
	    next if(/^\#/);
	    my @row = split;	    
	    my ($q,$h) = @row;
	    my ($sp) = split(/\|/,$h);
	    $sps{$sp}++;
	    my $bits = pop @row;
	    my $evalue = pop @row;
	    if( ! exists $table{$og} ||
		! exists $table{$og}->{$sp} ||
		$table{$og}->{$sp}->[1] > $evalue ) {
		$table{$og}->{$sp} = [ $h, $evalue];
	    }
	}
    }
}

my @splist = sort keys %sps;
for my $og ( keys %table ) {
    for my $sp ( @splist ) {
	print join("\t", $og, $sp, @{$table{$og}->{$sp} || []}),"\n";
    }
}
