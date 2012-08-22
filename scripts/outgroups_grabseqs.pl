#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;

my $dir = shift;
my $db = Bio::DB::Fasta->new($dir);

my %tbl;
while(<>) {
    my ($orth,$sp,$name,$evalue) = split;
    if( ! defined $name ) {
	next;
    }
    push @{$tbl{$orth}}, $name;
}

for my $s ( keys %tbl ) {
    my $out = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$s.orth.fasta");
    for my $s ( @{$tbl{$s}} ) {
	$out->write_seq($db->get_Seq_by_acc($s));
    }
}
