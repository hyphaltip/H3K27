#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;
use Text::NSP::Measures::2D::Fisher2::twotailed;
use List::Util qw(sum);

my $genomesdir = 'genomes';
my $rseg_genes = 'data/RSEGdomain_genes';

my $profile_proteome = 'data/phylo_profile/neurospora_crassa_or74a_10.all_proteomes.FASTA.m9-profile.pergene.tab.gz';
my $proteins = 'data/orthomcl/input/neurospora_crassa_or74a__finished__10_proteins.fasta';
my $odir = 'results/phylo_by_clade';
mkdir($odir);
my $max_pfam = 25;
my $do_rename = 1;

my $pfamext = 'tab';
my $hmmerversion = 3;
my $evalue_cutoff = 0.01;
my $pfamdir = 'data/pfam';
GetOptions('r|rename!'  => \$do_rename,
	   'p|pep:s'    => \$proteins,	
	   'pprofile:s' => \$profile_proteome,
	   'd|pfam:s'  => \$pfamdir,
	   'ext:s'     => \$pfamext,
	   'hmmer|hv:s'=> \$hmmerversion,
	   'e|evalue:s' => \$evalue_cutoff,
	   );

my $in = Bio::SeqIO->new(-format => 'fasta',
			 -file   => $proteins);
my %len;
while(my $s = $in->next_seq ) {
    my ($sp,$id) = split(/\|/,$s->id);
    $id =~ s/\.t\d+$//;
    $id =~ s/T\d+$//;
    $len{$id} = $s->length;
}

my %pfam;
opendir(PFAM, $pfamdir) || die "cannot open $pfamdir: $!";
for my $file ( readdir(PFAM) ) {
    next unless ( $file =~ /\.\Q$pfamext\E$/);
    open(my $fh => "$pfamdir/$file" ) || die $!;
    if( $hmmerversion == 3 ) {	#parse HMMER3 domtblout output
	while(<$fh>) {
	    next if /^\#/;
	    my ($domain,$acesssion,$tlen, $qname, $qacc,
		$qlen, $seq_evalue, $seq_score,$seq_bias,
		$n, $totalhits, $dom_cvalue, $dom_ievalue,$dom_score,
		$dom_bias) = split(/\s+/,$_);
	    $pfam{$qname}->{$domain}++ if $dom_cvalue < $evalue_cutoff;
	}
    } elsif( $hmmerversion == 2 ) { # parse HMMER2 (hmmer_to_table output)
	while(<$fh>) {
	    next if /^\#/;
	    my ($qname, $qstart, $qend, $domain,
		$domstart, $domend, $score, $evalue) = split;	       
	    $pfam{$qname}->{$domain}++ if $evalue < $evalue_cutoff;
	}
    } else { 
	warn("unknown HMMER version (2 or 3 are expected)\n");
    }
}

my %rename = ('Neuro-Sord-anid' => 'Asco',
	      'Ncrassa' => 'Ncrassa',
	      'Asco-Basidio' => 'Fungi',
	      'Asco-Basidio-Chytrid-Metazoa-Oomycota-Planta-Zygo' => 'Euk',
	      'Asco-Basidio-Chytrid-Metazoa-Zygo' => 'Animal',
	      'Asco-Basidio-Chytrid-Metazoa-Planta-Zygo' => 'Euk',
	      'Asco-Metazoa' => 'Animal',
	      
	      'Asco-Basidio-Metazoa-Planta' => 'Euk',
	      'Asco-Basidio-Chytrid-Metazoa-Planta' => 'Euk',
	      'Asco-Basidio-Chytrid-Metazoa-Planta-Zygo' => 'Euk',
	      'Asco-Metazoa-Planta ' => 'Euk',
	      'Asco-Basidio-Metazoa-Planta-Zygo' => 'Euk',
	      
	      'Asco-Basidio-Metazoa' => 'Animal',
	      'Asco-Basidio-Chytrid-Metazoa' => 'Animal',
	      'Asco-Basidio-Metazoa-Zygo' => 'Animal',
	      'Asco-Metazoa-Zygo' => 'Animal',
	      'Asco-Chytrid-Metazoa-Zygo' => 'Animal',
	      'Asco-Chytrid-Metazoa' => 'Animal',
	      'Asco-Basidio-Chytrid-Metazoa-Planta' => 'Euk',	      
	      
	      'Basidio-Sord' => 'Fungi',
	      'Asco-Basidio-Chytrid-Planta-Zygo' => 'Plant',
	      'Asco-Basidio-Planta-Zygo' => 'Plant',
	      'Asco-Basidio-Planta' => 'Plant',
	      'Asco-Basidio-Chytrid-Planta' => 'Plant',
	      'Asco-Chytrid-Planta' => 'Plant',
	      'Asco-Planta-Zygo' => 'Plant',
	      'Asco-Planta' => 'Plant',

	      'Asco-Basidio-Chytrid-Oomycota-Zygo' => 'Euk',
	      'Asco-Basidio-Chytrid-Oomycota' => 'Euk',
	      'Asco-Basidio-Oomycota' => 'Euk',
	      'Asco-Basidio-Oomycota-Zygo' => 'Euk',
	      'Asco-Chytrid-Oomycota-Zygo' => 'Euk',

	      'Asco-Metazoa-Oomycota-Planta-Zygo' => 'Euk',
	      'Asco-Chytrid-Oomycota-Planta' => 'Euk',
	      'Asco-Metazoa-Oomycota-Zygo' => 'Euk',
	      'Asco-Metazoa-Oomycota-Planta' => 'Euk',
	      'Asco-Basidio-Chytrid-Oomycota-Planta-Zygo' => 'Euk',
	      'Asco-Basidio-Chytrid-Metazoa-Oomycota-Planta' => 'Euk',
	      'Asco-Basidio-Metazoa-Oomycota-Planta-Zygo' => 'Euk',
	      'Asco-Basidio-Chytrid-Metazoa-Oomycota-Zygo' => 'Euk',
	      'Asco-Basidio-Oomycota-Planta-Zygo' => 'Euk',
	      'Asco-Basidio-Chytrid-Oomycota-Planta' => 'Euk',
	      'Asco-Chytrid-Metazoa-Oomycota-Planta-Zygo' => 'Euk',
	      'Asco-Chytrid-Metazoa-Oomycota-Zygo' => 'Euk',
	      'Asco-Basidio-Chytrid-Oomycota-Planta' => 'Euk',
	      'Asco-Basidio-Metazoa-Oomycota' => 'Euk',
	      'Asco-Basidio-Metazoa-Oomycota-Zygo' => 'Euk',
	      'Asco-Basidio-Chytrid-Metazoa-Oomycota' => 'Euk',
	      'Asco-Chytrid-Metazoa-Planta' => 'Euk',
	      'Asco-Basidio-Oomycota-Planta' => 'Euk',
	      'Asco-Basidio-Metazoa-Oomycota-Planta' => 'Euk',
	      'Asco-Basidio-Oomycota-Planta' => 'Euk',
	      'Asco-Metazoa-Planta' => 'Euk',

	      'Chytrid-Sord' => 'Fungi',
	      'Asco-Chytrid-Zygo' => 'Fungi',
	      'Asco-Basidio-Chytrid-Zygo' => 'Fungi',
	      'Asco-Basidio-Zygo' => 'Fungi',
	      'Asco-Basidio-Chytrid' => 'Fungi',
	      'Asco-Chytrid' => 'Fungi',
	      'Asco-Zygo' => 'Fungi',
	      'Asco' => 'Asco',
	      'Sord' => 'Sord',
	      'Neuro' => 'Neuro',
	      );


my %pprof;
my $prf;
if( $profile_proteome =~ /\.gz/ ) {
 open($prf => "zcat $profile_proteome |") || die "$profile_proteome: $!";
} else {
 open($prf => $profile_proteome) || die "$profile_proteome: $!";
}
my $header1 = <$prf>;
while(<$prf>) {
    my ($gene,$clade) = split;
    my ($sp,$gene_name) = split(/\|/,$gene);
    $gene_name =~ s/\.t\d+$//;
    $gene_name =~ s/T\d+$//;
    $clade = $rename{$clade} || 'Other' if $do_rename;
    $pprof{$gene_name} = [$clade, join(",", keys %{$pfam{$gene} || {}}) ];
}
close($prf);

my %K27;
my %count_cat;
my %cats;
opendir(DIR,$rseg_genes) || die "$rseg_genes: $!";
for my $spname ( readdir(DIR) ) {
    next if ($spname =~ /^\./);
    next unless $spname eq 'N_crassa'; # ONLY KEEPING N.CRASSA FOR NOW
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
	    push @{$cats{$type."_".$cat}->{$sp}}, $gene;
	}
    }
}

my $i = 0;
my (%clades,%percentile_pfam, %percentile_name);
for my $category ( keys %cats ) {
    for my $gene_name ( @{$cats{$category}->{Nc}}  ) {	
	my ($clade,$pfam) = @{$pprof{$gene_name} || [qw(singleton)]};
	$clades{$category}->{$clade}++;    
    
	for my $p ( split(/,/,$pfam || '') ) {
	    $percentile_pfam{$category}->{$p}++;
	}
	push @{$percentile_name{$category}}, [$gene_name,
					      $len{$gene_name} || 
						warn"cannot find $gene_name\n"];
    
    }
}

open(my $R => ">$odir/lens.R") || die $!; 
print $R "pdf(\"K27_peplens.pdf\")\n";
{ 
    my (@l,@names);
    for my $p ( sort keys %percentile_name ) {
	open(my $lfh => ">$odir/$p\_lengths.dat") || die $!;
	print $lfh join("\t", qw(GENE LEN)),"\n";
	print $R "l$p <- read.table(\"$p\_lengths.dat\",header=T,sep=",'"\t"',");\n";
	for my $g ( @{$percentile_name{$p}} ) {
	    print $lfh join("\t", @$g),"\n";
	}       
	print $R "summary(l$p)\n";
	push @l, "subset(l$p\$LEN,l$p\$LEN < 4000)";
	push @names, "l$p";
    }
    print $R "boxplot(",join(",",@l),",xlab=\"K27 Status\", main=\"Protein size for K27 status\",ylab=\"Protein length\",range=2,names=c(",join(",",map { sprintf('"%s"',$_) } @names),"),col=rainbow(11, start=0, end=1));\n";
}
close($R);
open($R => ">$odir/clade.R") || die $!; 
print $R "pdf(\"K27_byclades.pdf\");\n";
for my $p ( sort keys %clades ) {
    open(my $pfh => ">$odir/$p\_clade.dat") || die $!;     
    print $R "c$p <- read.table(\"$p\_clade.dat\",header=F,sep=",'"\t"',");\n";
    for my $clade ( sort { $clades{$p}->{$b} <=> $clades{$p}->{$a}} keys %{$clades{$p}} ) {
	print $pfh join("\t", $clade, $clades{$p}->{$clade}),"\n";	
    }
    print $R "pct <- round(c$p\$V2/sum(c$p\$V2)*100);\n";
    print $R "pct <- round(c$p\$V2/sum(c$p\$V2)*100);\n";
    print $R "lbls <- paste(c$p\$V1, pct);\n"; # add percents to labels 
    print $R "lbls <- paste(lbls,\"%\",sep=\"\");\n"; # add % to labels 
    print $R "par(mar=c(10,0,2,2))\n";
    print $R "pie(c$p\$V2,labels=lbls,col=rainbow(length(lbls)+1),radius=0.8,main=\"K27 $p Conservation\");\n";
    print $R "par(mar=c(10,4,4,2))\n";
    print $R "barplot(c$p\$V2,names.arg=lbls,las=2,ylab=\"count\",main=\"K27 $p % (Conservation Counts)\");\n";

} 
close($R);

open($R => ">$odir/pfam.R") || die $!;
print $R "pdf(\"K27_bypfam.pdf\");\n";

for my $p ( sort keys %percentile_pfam) {
    open(my $ppfha => ">$odir/$p\_pfam_all.dat") || die $!;
    open(my $ppfh => ">$odir/$p\_pfam.dat") || die $!;
    my $i =0;
    for my $pfam ( sort { $percentile_pfam{$p}->{$b} <=> $percentile_pfam{$p}->{$a}} keys %{$percentile_pfam{$p}} ) {
	print $ppfha join("\t", $pfam, $percentile_pfam{$p}->{$pfam}),"\n";
	print $ppfh join("\t", $pfam, $percentile_pfam{$p}->{$pfam}),"\n" if $i++ < $max_pfam;

    }
    print $R "p$p <- read.table(\"$p\_pfam.dat\",header=F,sep=",'"\t"',");\n";    
    print $R "pct <- round(p$p\$V2/sum(p$p\$V2)*100);\n";
    print $R "pct <- round(p$p\$V2/sum(p$p\$V2)*100);\n";
    print $R "lbls <- paste(p$p\$V1, pct);\n"; # add percents to labels 
    print $R "lbls <- paste(lbls,\"%\",sep=\"\");\n"; # add % to labels 
    print $R "par(mar=c(10,4,4,2))\n";
    print $R "barplot(p$p\$V2,names.arg=lbls,las=2,ylab=\"count\",main=\"K27 $p% (Pfam counts)\");\n";
    print $R "par(mar=c(10,0,2,2))\n";
    print $R "pie(p$p\$V2,labels=lbls,col=rainbow(length(lbls)),main=\"K27 $p% Pfam\");\n";

}
