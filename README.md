# FRT
code from Male-Derived Transcripts Isolated From Mated Female Reproductive Tract in Drosophila melanogaster

This contains the set of custom perl scripts used to identify distinguishing SNPs between RAL genotypes and use them to identify full length transcripts from the male.

The order of the scripts is as follows:

1) Identify_SNPs.pl
2) Identify_genes_for_each_site.pl
3) label_SNPs_by_category.pl
4) run_bcf_for_ALLSNPs.sh
5) compare_dgrpSNPs_to_vcfl_ALL.pl
6) find_unexpected_hets.pl
7) run_masked.sh - runs copies of same_snps_in_masked_regions.pl for multiple input files
8) compare_transcripts_between_crosses.pl
9) add_label_to_combined_counts.pl
10) get_intron_exon_longest_CDS.pl
11) get_intron_exon_longest_transcript.pl
12) get_transcript_cov_from_SNPs.pl
13) get_FRT_high_confidence_expressed_transcripts.pl
14) make_reference_specific_transcripts.pl
15) 
