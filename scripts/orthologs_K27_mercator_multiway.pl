#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use List::Util qw(sum max min);
my $ref = 'Nc';
my $debug = 0;
# use to 10% of genes
my $top_cutoff = 0.10;
my %files = ('Nc' => 'data/RSEGdomain_genes/N_crassa/Nc_H3K27me3_genes.txt',
	     'Nd' => 'data/RSEGdomain_genes/N_discreta/Nd_H3K27me3_genes.txt',
	     'Nt' => 'data/RSEGdomain_genes/N_tetrasperma/Nt_H3K27me3_genes.txt');

my $fix_names = 'data/annotation_fix/v4junev5.csv';
my $orthologs = 'data/syntenic_orthologs/Nc-Nt-Nd.orthologs';

GetOptions(
	   'v|verbose!' => \$debug,
	   'f|fix:s'    => \$fix_names,
	   'orth:s'     => \$orthologs,
	   );

my %fixnames;
open(my $fh => $fix_names) || die $!;
while(<$fh>) {
    next if /^\#/;
    chomp;
    my ($old,$new) = split(/,/,$_);
    $old =~ s/\.\d+$//;
    $new =~ s/\.\d+$//;
    next if( $old eq $new);
    $fixnames{$old} = $new;
}

my %groups;
my %other_groups;
open($fh => $orthologs) || die $!;
my (undef,undef,$filebase) = File::Spec->splitpath($orthologs);
my ($name) = split(/\./,$filebase);
my ($from,@to) = split(/-/,$name);
for ( @to ) { $other_groups{$_}++ }
my $header = <$fh>;
chomp($header);
my @header = split(/\t/,$header);
my %genes;
while(<$fh>) {
    while(<$fh>) {
	chomp;
	my (@row) = split(/\t/,$_);
	my (undef,$gene) = @row;
	$gene =~ s/T\d+$//;
	my $dat = {};
	for( my $i = 0; $i < @header; $i++ ) {
	    last if ! defined $row[$i];
	    my ($t,$src) = split(/_/,$header[$i]);
	    if( $t =~ /MRNA|GENE/ ) {
		if( $src eq 'Nc') {
		    if( $t eq 'MRNA') {
			$row[$i] =~ s/T\d+//; 
		    }		    
		}
		$row[$i] = "$src|".$row[$i];
	    }
	    if( $src eq $ref ) {
		$genes{"$from|$gene"}++;
	    }
	    $groups{$from}->{$src}->{"$from|$gene"}->{$header[$i]} = $row[$i];
	}
    }
}
my %is_K27;
for my $species ( keys %files ) {
   open(my $fh => $files{$species}) || die $files{$species} . ": $!";
    my $header = <$fh>;
    while(<$fh>) {
	chomp;	
	my ($locus,$scaffold,$start,$end,$strand) = split(/\t/,$_);
	if( $fixnames{$locus} ) {
	    $locus = $fixnames{$locus};		    
	}
	$is_K27{$species}->{"$species|$locus"} = [$scaffold,$start,$end,$strand,];
    }
}

my $target = $ref;
my (@others) = keys %other_groups;
open(my $ofh_stat => ">K27_multigenestats.txt") || die $!;
my %counts = ( 'no_ortholog' => [ 0, 0] );
    
open(my $ofh => ">K27_multigeneinfo.tab") || die $!;
print $ofh join("\t", map { my $sp = $_; 
			    map { sprintf("%s_%s",$_,$sp) }  qw(MRNA K27) }
		$target,@others), "\n";
my %all;
for my $nm ( (keys %genes), (keys %{$is_K27{$target}}) ) {
    $all{$nm}++;
}
for my $gene (sort keys %all) {
#    $gene =~ s/T\d+$//;
    my $is_K27_marked = exists $is_K27{$target}->{$gene} ? 'yes' : 'no';;
    
    my @row = ( $gene,$is_K27_marked);
    
    my @K27_aggregate = ($is_K27_marked);
    my %noorth;
    for my $other (@others) {
	my $dat = $groups{$target}->{$other}->{$gene};
	if( ! $dat ) {
	    warn("No listing for $gene in $other\n");
	    push @row, 'NONE','NA';
	    next;	
	}
	if($dat->{"STRAND_$other"} =~ /^NO_/) {
	    warn "no $other ortholog for $gene K27 status: $is_K27_marked\n" if $debug;
	    $noorth{$other} = 1;
	    #$counts{'no_ortholog'}->[$is_K27_marked eq 'yes' ? 0 : 1]++;
	    # need to count these
	    push @row, 'NONE', 'NA';
	    next;
	}
	
	my @ortholog_K27_marked;
	my @to_nms;
	my $to_gene = $dat->{"MRNA_$other"};
	my %c;

	if( ! exists $is_K27{$other} ) {
	    warn("no K27 available for $other\n");
	}
	for my $to_gene_nm ( split(/,/,$to_gene) ) {
	    push @ortholog_K27_marked, exists $is_K27{$other}->{$to_gene_nm} ? 'yes' : 'no';
	    push @to_nms, $to_gene_nm;
	}
	push @row, join(",", @to_nms), join(",",@ortholog_K27_marked);
	for my $r (@ortholog_K27_marked ) { # figure out the unique pattern of K27 for this sp, for single-copy it is just 1 value 'yes' or 'no'
	    $c{$r}++;
	}
	push @K27_aggregate, join("-",sort keys %c); # summarize the patterns seen for this 1->n ortholog (most are 1->2 where gene is split in Nt or Nd
    }
    if( (scalar keys %noorth) == scalar @others ) {
	$counts{'no_ortholog'}->[$is_K27_marked eq 'yes' ? 0 : 1]++;
    }
    print $ofh join("\t",@row),"\n";
    $counts{'K27'}->{join(",",@K27_aggregate)}++;    
}
    
#1 are the no-ortholog ones more or less likely to be K27 marked?
printf $ofh_stat 
    " ($target) %d genes had no ortholog %d K27 marked, %d NOT K27 marked\n",
    sum (@{$counts{'no_ortholog'}}),@{$counts{'no_ortholog'}};
    
my @patterns = keys %{$counts{'K27'}};
printf $ofh_stat "\nThe breakdown of genes as K27 marks in:\n$target\t%s\tGENE COUNT\n",join("\t",@others);
# 2 what are the GO terms for the commonly K27 marked
#
# 3 What groups are the K27 commonly marked ones in -- which phylogenetic profile group?
for my $pat ( sort @patterns ) {
    print $ofh_stat join("\t",split(/,/,$pat),$counts{'K27'}->{$pat}), "\n";
}
