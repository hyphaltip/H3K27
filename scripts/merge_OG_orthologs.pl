#!/usr/bin/perl -w
use Bio::SeqIO;
my $dir1 = 'orthomclI20.outgroup_orthologs';
my $dir2 = 'orthomcl_orthologs_I20';
my $odir = 'combined_I20';

mkdir($odir);

my %og;

opendir(D1,$dir1) || die $!;
for my $file ( readdir(D1) ) {
    if( $file =~ /(\S+)\.orth\.fasta$/) {
	push @{$og{$1}}, "$dir1/$file";
    }
}

opendir(D2,$dir2) || die $!;
for my $file ( readdir(D2) ) {    
    if( $file =~ /(\S+)\.fa$/) {
	if( ! exists $og{$1} ) {
	    warn("skipping '$1', it didn't have an outgroup file\n");
	    next;
	}	    
	push @{$og{$1}}, "$dir2/$file";
    }
}

for my $set (keys %og ) {
    my $out = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$odir/$set.fas");
    for my $set ( @{$og{$set}} ) {
	my $in = Bio::SeqIO->new(-format => 'fasta',
				 -file   => $set);
	while( my $s = $in->next_seq ) {
	    $out->write_seq($s);
	}
    }
}
