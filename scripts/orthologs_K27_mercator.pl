 #!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use List::Util qw(sum max min);

# use to 10% of genes
my $top_cutoff = 0.10;
my %files = ('Nc' => 'data/RSEGdomain_genes/N_crassa/Nc_H3K27me3_genes.txt',
	     'Nd' => 'data/RSEGdomain_genes/N_discreta/Nd_H3K27me3_genes.txt',
	     'Nt' => 'data/RSEGdomain_genes/N_tetrasperma/Nt_H3K27me3_genes.txt');

my @orthologs = ( 'data/syntenic_orthologs/Nc-Nd.orthologs.gz',
		  'data/syntenic_orthologs/Nc-Nt.orthologs.gz',
# 'data/syntenic_orthologs/Nt-Nd.orthologs',
		  );

GetOptions();

my %groups;
my %other_groups;
for my $file ( @orthologs ) {
    my $fh;
    if( $file =~ /\.gz$/ ) {
     open($fh => "zcat $file |") || die $!;
    } else {
     open($fh => $file) || die $!;
    }
    my (undef,undef,$filebase) = File::Spec->splitpath($file);
    my ($name) = split(/\./,$filebase);
    my ($from,$to) = split(/-/,$name);
    $other_groups{$to}++;
    my $header = <$fh>;
    chomp($header);
    my @header = split(/\t/,$header);
    while(<$fh>) {
	chomp;
	my (@row) = split(/\t/,$_);
	my (undef,$gene) = @row;
	my $dat = {};
	for( my $i = 0; $i < @header; $i++ ) {
	    last if ! defined $row[$i];
	    if( $header[$i] eq 'MRNA_TO' ||
		$header[$i] eq 'GENE_TO' ) {
		$row[$i] = "$to|".$row[$i];
	    }
	    $dat->{$header[$i]} = $row[$i];
	}
	push @{$groups{$from}->{$to}->{"$from|$gene"}}, $dat;	
    }
}
my %is_K27;
for my $species ( keys %files ) {
   open(my $fh => $files{$species}) || die $!;
    my $header = <$fh>;
    while(<$fh>) {
	chomp;
	my ($locus,$scaffold,$start,$end,$strand) = split(/\t/,$_);
	if( $species eq 'Nc' ) {
	    $locus =~ s/\.4$/T0/;
	    if( $locus !~ /T\d+$/ ) { # a real hack
		$locus .= "T0";
	    }
	}
	$is_K27{$species}->{"$species|$locus"} = [$scaffold,$start,$end,$strand,];
    }
}

my $target = 'Nc';
my (@others) = keys %other_groups;
for my $other ( @others ) {
    open(my $ofh_stat => ">$other\_K27_genestats.txt") || die $!;
    my %counts = ( 'no_ortholog' => [ 0, 0] );
    
    open(my $ofh => ">$other\_K27_geneinfo.tab") || die $!;
    print $ofh join("\t", qw(MRNA_FROM K27_MARKED_FROM K27_VAL_FROM SPECIES_TO MRNA_TO)), "\n";
    for my $gene (keys %{$is_K27{$target}}) {
	my $is_K27_marked = exists $is_K27{$target}->{$gene} ? 'yes' : 'no';;
	my @dat = @{$groups{$target}->{$other}->{$gene} || []};
	if( ! @dat ) {
	    warn("No listing for $gene\n");
	    #next;
	}
	for my $dat ( @dat ) {
	    if($dat->{STRAND_TO} =~ /^NO_/) {
		warn "no $other ortholog for $gene K27 status: $is_K27_marked\n";
		$counts{'no_ortholog'}->[$is_K27_marked eq 'yes' ? 0 : 1]++;
		# need to count these
		next;
	    } 
	    
	    my $ortholog_K27_marked = 'no';
	    my $ortholog_K27_val = "undef";
	    my $to_gene = $dat->{MRNA_TO};

	    if( exists $is_K27{$other} &&
		exists $is_K27{$other}->{$to_gene} ) {
		$ortholog_K27_marked = exists $is_K27{$other}->{$to_gene} ? 'yes' : 'no';
	    }
	    print $ofh join("\t", $dat->{MRNA_FROM}, 
			    $is_K27_marked,
			    $is_K27{$target}->{$gene}->[0],
			    $other, $dat->{MRNA_TO}, 
			    $ortholog_K27_marked,
			    ), "\n";
	    $counts{'K27_Ncra'}->{"$is_K27_marked,$ortholog_K27_marked"}++;
	}
    }

    # finish counting non-K27 marked target species genes to get the number of 
    # of the non-marked genes that also have orthologs
    for my $gene ( keys %{$groups{$target}->{$other}}) {
	next if $is_K27{$target}->{$gene};
	my @dat = @{$groups{$target}->{$other}->{$gene} || []};
	warn("here with $gene\n");
	for my $dat ( @dat ) {
	    $counts{'no_ortholog'}->[1]++ if $dat->{'STRAND_TO'} =~ /^NO_/;
	}
    }
    
    #1 are the no-ortholog ones more or less likely to be K27 marked?
    printf $ofh_stat 
	" ($target) %d genes had no ortholog %d K27 marked, %d NOT K27 marked\n",
	sum (@{$counts{'no_ortholog'}}),@{$counts{'no_ortholog'}};
    
    my @patterns = keys %{$counts{'K27_Ncra'}};
    printf $ofh_stat "\nThe breakdown of genes as K27 marks in:\n$target\t$other\tGENE COUNT\n";
    # 2 what are the GO terms for the commonly K27 marked
    # 3 What groups are the K27 commonly marked ones in -- which phylogenetic profile group?
    for my $pat ( sort @patterns ) {
	print $ofh_stat join("\t",split(/,/,$pat),$counts{'K27_Ncra'}->{$pat}), "\n";
    }

}

# merge the 3
for my $from ( keys %groups ) {
    for my $to ( keys %{$groups{$from}} ) {
	for my $gene (keys %{$groups{$from}->{$to}} ) {
	    
    }
    push @{$groups{$from}->{$to}->{"$from|$gene"}}, $dat;	
}
