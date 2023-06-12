#!sh

bcftools mpileup -q 30 --regions-file ALL_SNPs.bed -f /data/FlyRef/Drosophila_melanogaster/6.41/fasta/dmel-all-chromosome-r6.41.fasta align/304_SR.bam | bcftools call -c -o /data/julie/FemaleRT/304_SR_ALL.vcf &

bcftools mpileup -q 30 --regions-file ALL_SNPs.bed -f /data/FlyRef/Drosophila_melanogaster/6.41/fasta/dmel-all-chromosome-r6.41.fasta align/304_ST.bam  | bcftools call -c -o /data/julie/FemaleRT/304_ST_ALL.vcf &

bcftools mpileup -q 30 --regions-file ALL_SNPs.bed -f /data/FlyRef/Drosophila_melanogaster/6.41/fasta/dmel-all-chromosome-r6.41.fasta align/304_PV.bam  | bcftools call -c -o /data/julie/FemaleRT/304_PV_ALL.vcf &

bcftools mpileup -q 30 --regions-file ALL_SNPs.bed -f /data/FlyRef/Drosophila_melanogaster/6.41/fasta/dmel-all-chromosome-r6.41.fasta align/307_SR.bam | bcftools call -c -o /data/julie/FemaleRT/307_SR_ALL.vcf &

bcftools mpileup -q 30 --regions-file ALL_SNPs.bed -f /data/FlyRef/Drosophila_melanogaster/6.41/fasta/dmel-all-chromosome-r6.41.fasta align/307_ST.bam | bcftools call -c -o /data/julie/FemaleRT/307_ST_ALL.vcf &

bcftools mpileup -q 30 --regions-file ALL_SNPs.bed -f /data/FlyRef/Drosophila_melanogaster/6.41/fasta/dmel-all-chromosome-r6.41.fasta align/307_PV.bam | bcftools call -c -o /data/julie/FemaleRT/307_PV_ALL.vcf &

bcftools mpileup -q 30 --regions-file ALL_SNPs.bed -f /data/FlyRef/Drosophila_melanogaster/6.41/fasta/dmel-all-chromosome-r6.41.fasta align/360_SR.bam | bcftools call -c -o /data/julie/FemaleRT/360_SR_ALL.vcf & 

bcftools mpileup -q 30 --regions-file ALL_SNPs.bed -f /data/FlyRef/Drosophila_melanogaster/6.41/fasta/dmel-all-chromosome-r6.41.fasta align/360_ST.bam | bcftools call -c -o /data/julie/FemaleRT/360_ST_ALL.vcf &

bcftools mpileup -q 30 --regions-file ALL_SNPs.bed -f /data/FlyRef/Drosophila_melanogaster/6.41/fasta/dmel-all-chromosome-r6.41.fasta align/399_ST.bam | bcftools call -c -o /data/julie/FemaleRT/399_ST_ALL.vcf &
