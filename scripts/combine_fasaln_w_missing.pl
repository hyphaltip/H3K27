#!/usr/bin/perl -w
use strict;
use List::Util qw(shuffle);
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Getopt::Long;
my $iformat = 'fasta';
my $oformat = 'nexus';
my $outfile = 'allseq.nex';
my $ext = 'fasaln.trim';
my $dir;
my $max = 100;
my $GAP = '-';
my %skip = ('amac' => 1, 'mcin' => 1, 'abra' => 1, 'snod' => 1, 'post_PC15' => 1, 'tmes' => 1, 'calb_sc5314' => 1, 'cgui' => 1, 'cglo' => 1, 'BdenJEL423' => 1,'pbla' => 1,
	    'rgra' => 1, 'spun' => 1);
GetOptions('d|dir:s'   => \$dir,
	   'ext:s'     => \$ext,
	   'o|out:s'   => \$outfile,
	   'm|max:i'   => \$max,
	   );

die("need a dir") unless $dir && -d $dir;

opendir(DIR, $dir) || die"$dir: $!";

my (%matrix);
my $count = 0;
my @files;
for my $file (shuffle readdir(DIR) ) {
    next if $file eq $outfile;
    next unless ($file =~ /(\S+)\.\Q$ext\E$/);
    my $in = Bio::AlignIO->new(-format => $iformat,
			       -file   => "$dir/$file");
    if( my $aln = $in->next_aln ) {
	next if $aln->length == 0;
	for my $seq ( $aln->each_seq ) {
	    $matrix{$seq->id} = '';
	}
    }
    push @files, $file;
}
my @all = keys %matrix;

for my $file (@files ) {
    next if $file eq $outfile;
    next unless ($file =~ /(\S+)\.\Q$ext\E$/);
    my $in = Bio::AlignIO->new(-format => $iformat,
			       -file   => "$dir/$file");
    if( my $aln = $in->next_aln ) {
	next if $aln->length == 0;
	my %seen;
	for my $seq ( $aln->each_seq ) {
	    $matrix{$seq->id} .= $seq->seq;
	    $seen{$seq->id}++;
	}
	for my $id ( @all ) {
	    if( ! $seen{$id} ) {
		$matrix{$id} .= $GAP x $aln->length;
		$seen{$id}++;
	    }
	}
    }
    last if $count++ > $max;
}

my $bigaln = Bio::SimpleAlign->new;
while( my ($id,$seq) = each %matrix ) {
    next if $skip{$id};
    $bigaln->add_seq(Bio::LocatableSeq->new(-id  => $id,
					    -seq => $seq));
}

my $out = Bio::AlignIO->new(-format => $oformat,
			    -file   => ">$outfile");
$out->write_aln($bigaln);


