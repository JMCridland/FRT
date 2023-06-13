# FRT
code from Male-Derived Transcripts Isolated From Mated Female Reproductive Tract in Drosophila melanogaster

This contains the set of custom perl scripts used to identify distinguishing SNPs between RAL genotypes and use them to identify full length transcripts from the male.

The order of the scripts is as follows:

1) Identify_SNPs.pl - This script uses input from the DGRP SNPs and positions of genes from Flybase and produces a table with all SNPs between the focal lines.

2) Identify_genes_for_each_site.pl - This script takes the output from Identify_SNPs.pl and labels them by the gene they are associated with.

3) label_SNPs_by_category.pl - This script takes the output from Identify_genes_for_each_site.pl and labels them by context; i.e. UTR, exon, intron. 

4) run_bcf_for_ALLSNPs.sh - 

5) compare_dgrpSNPs_to_vcfl_ALL.pl

6) find_unexpected_hets.pl

7) same_snps_in_masked_regions.pl 

8) compare_transcripts_between_crosses.pl

9) add_label_to_combined_counts.pl

10) get_intron_exon_longest_transcript.pl

12) get_transcript_cov_from_SNPs.pl

13) get_FRT_high_confidence_expressed_transcripts.pl

###At this step parent specific transcripts are generated and used to improve the transcript coverage for genes.
14) make_reference_specific_transcripts.pl 

15) update_high_confidence_with_psr_reads.pl

16) make_FRT_abund.pl


