#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Text::NSP::Measures::2D::Fisher2::twotailed;
use List::Util qw(sum);

my $genomesdir = 'genomes';
my $rseg_genes = 'Neurospora_RSEGdomain_genes';


opendir(DIR, $genomesdir) || die $!;
my %signalp;
for my $file ( readdir(DIR) ) {
    if( $file =~ /(\S+)\.signalp/) {
	my $sp = $1;
	open(my $fh => "$genomesdir/$file") || die $!;
	my @header;

	while(<$fh>) {
	    if(/^\#\s+name/ ) {
		s/^\#//;
		@header = split;
	    } elsif( /^\#/ || /^Locus/ ) {
		next;
	    } else {
		my @row = split;
		$row[0] = (split(/\|/,$row[0]))[-1]; # last thing in | separate set is the gene name used
		$row[0] =~ s/(NCU\d+)T\d+$/$1/; # drop transcript info
		my $name = $row[0];		
		for my $h ( @header ) {		   
		    $signalp{$sp}->{$name}->{$h} = shift @row;
		}
	    }
	}
    }
}

my %K27;
my %count_cat;
my %cats;
opendir(DIR,$rseg_genes) || die $!;
for my $spname ( readdir(DIR) ) {
    next if ($spname =~ /^\./);
    opendir(SP,"$rseg_genes/$spname" ) || die $!;
    for my $file ( readdir(SP) ) {
	next unless $file =~ /^(\S+)\.txt/;
	my $base = $1;
	my ($sp,$type,$cat) = split(/_/,$base);
	open(my $fh => "$rseg_genes/$spname/$file" ) || die $!;
	while(<$fh>) {
	    chomp;
	    next if /^Gene|Locus/;
	    my ($gene,$chrom,$start,$end,$strand,$score) = split(/\t/,$_);	    
	    next if ! defined $gene;
	    if( exists $K27{$sp}->{$gene} ) {
		warn("already seen '$gene' ($sp) with value ",
		     $K27{$sp}->{$gene}, " will set it to ",$type."_".$cat,"\n");
	    }
	    $K27{$sp}->{$gene} = $type."_".$cat;
	    $cats{$type."_".$cat}++;
	    if( exists $signalp{$sp}->{$gene} ) {		 
		my $signalp_status = $signalp{$sp}->{$gene}->{'?'};
		$count_cat{$sp}->{$type."_".$cat}->{$signalp_status}++;
	    } else {
		warn("gene ($sp) $gene does not exist\n");
	    }
	}
    }
}

my @cats = sort keys %cats;
print join("\t", qw(SPECIES SIGNALP), @cats), "\n";
for my $sp ( keys %count_cat ) {
    print join("\t", $sp, 'Y', map { $count_cat{$sp}->{$_}->{'Y'} } @cats ), "\n";    
    print join("\t", $sp, 'N', map { $count_cat{$sp}->{$_}->{'N'} } @cats ), "\n";    
}


my @cat_compare = ( [qw(H3K27me3_genes nonH3K27me3_genes)],
		    [qw(H3K27me3_bordergenes nonH3K27me3_genes)],
		    [qw(H3K27me3_bordergenes H3K27me3_genes)]);
		    
#               CAT_COMPARE	
#           word2   ~word2
#  word1    n11      n12 | n1p   (YES)
# ~word1    n21      n22 | n2p   (NO)
#           --------------
#           np1      np2   npp	     

for my $sp ( keys %count_cat ) {
    print "\nSpecies: $sp\n";
    for my $cat ( @cat_compare ) {
	print join("\t", qw(SIGNALP), @$cat), "\n";
	my ($n11) = $count_cat{$sp}->{$cat->[0]}->{Y};
	my ($n21) = $count_cat{$sp}->{$cat->[0]}->{N};

	my ($n12) = $count_cat{$sp}->{$cat->[1]}->{Y};
	my ($n22) = $count_cat{$sp}->{$cat->[1]}->{N};


	my $n1p = $n11 + $n12;
	my $np1 = $n11 + $n21;

	my $npp = sum($n11,$n21,$n12,$n22);
	
	my $twotailed = calculateStatistic( n11=>$n11,
						  n1p=>$n1p,
						  np1=>$np1,
						  npp=>$npp);

	print join("\t", 'Y', $n11, $n21), "\n";
	print join("\t", 'N', $n12, $n22), "\n";
	print "--\n";	

	printf "n11=%d n1p=%d np1=%d npp=%d; \n\tP-value=%.4g\n\n",
	$n11,$n1p,$np1,$npp, $twotailed;
    }
}



