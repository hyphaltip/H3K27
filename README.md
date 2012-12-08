H3K27
=====

H3K27 neurospora project

scripts/
 - combine_fasaln_w_missing.pl for building a multiple alignment for phylogenetic tree, even with missing data
 - compare_K27_old_list.pl for simple tabulation of presence absence in the K27 list
 - count_genes_signalp_K27.pl for tabulating number of genes with K27 mark and signalP domains
 - extract_multiway_orthologs_mercator.pl for computing synteny based ortholog map
 - gatherFASTA_2orthogroups.pl for summarizing the ortholog group 
 - merge_OG_orthologs.pl combine the Neurospora-specific ortholog groups with Outgroup containing OG to build fasta filecontaining all
 - orthologs_K27_mercator_multiway.pl MAIN script computing the K27 status of orthologs and summarizing the output, doing it for the multi-way orthologs (3-way) to assess conservation and K27 status
 - orthologs_K27_mercator.pl - pairwise comparison of orthologs (e.g. how many Nc-Nd orthologs are K27 marked in both species, only one, etc)
 - orthomcl_phylogenetic_profile.pl for profiling species distribution for each ortholog group so that in total the phylogenetic profile of each Nc gene can be computed
 - outgroups_grabseqs.pl Generate fasta sequence files for ortholog groups 
 - phyloprofile.pl for a simpler phylogenetic profile based on one-direction FASTA hits
 - syntenic_ortholog_stats.pl to generate summary statistics about the number orthologs from the multiway ortholog report (perl scripts/syntenic_ortholog_stats.pl data/syntenic_orthologs/Nc-Nt-Nd.orthologs.gz)
