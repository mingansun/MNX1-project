#!/bin/sh
## Ming-an Sun, 2021-08-11

###############################################################################
# Pipelie for routine analysis and visualization of RNA-Seq data
###############################################################################

## set input file informatio
# RNA-Seq file
wt_rep1_r1=WT_rep1_1.fq
wt_rep1_r2=WT_rep1_2.fq
wt_rep2_r1=WT_rep2_1.fq
wt_rep2_r2=WT_rep2_2.fq
ko_rep1_r1=KO_rep1_1.fq
ko_rep1_r2=KO_rep1_2.fq
ko_rep2_r1=KO_rep2_1.fq
ko_rep2_r2=KO_rep2_2.fq
# Reference genome
star_index=/data/home/masun/db/StarIndex/GRCm38
rsem_index=/data/home/masun/db/RsemIndex/GRCm38
gene_gtf_file=/data/home/masun/db/Gene/GRCm38.gtf
chrom_size_file=/data/home/masun/db/Genome/GRCm38.chrom.sizes

## Trim fastq file using Trim Galore
wt_rep1_r1_trim=WT_rep1_1_val_1.fq
wt_rep1_r2_trim=WT_rep1_2_val_2.fq
wt_rep2_r1_trim=WT_rep2_1_val_1.fq
wt_rep2_r2_trim=WT_rep2_2_val_2.fq
ko_rep1_r1_trim=KO_rep1_1_val_1.fq
ko_rep1_r2_trim=KO_rep1_2_val_2.fq
ko_rep2_r1_trim=KO_rep2_1_val_1.fq
ko_rep2_r2_trim=KO_rep2_2_val_2.fq

trim_galore --illumina --paired wt_rep1_r1 wt_rep1_r2
trim_galore --illumina --paired wt_rep2_r1 wt_rep2_r2
trim_galore --illumina --paired ko_rep1_r1 ko_rep1_r2
trim_galore --illumina --paired ko_rep2_r1 ko_rep2_r2

## Perform alignment, filtering, sorting and PCR duplicate removal
wt_rep1_bam = wt_rep1.bam
# WT
STAR --runThreadN 12 --genomeDir $star_index --genomeLoad NoSharedMemory --readFilesIn $wt_rep1_r1_trim $wt_rep1_r2_trim \
	--outSAMtype BAM SortedByCoordinate --outStd SAM --outSAMattributes Standard --outSAMstrandField intronMotif --outSAMunmapped None \
	--outFilterType BySJout outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix wt_rep1
mv wt_rep1.Aligned.sortedByCoord.out.bam $wt_rep1.bam
samtools index ${wt_rep1}.bam
STAR --runThreadN 12 --genomeDir $star_index --genomeLoad NoSharedMemory --readFilesIn $wt_rep2_r1_trim $wt_rep2_r2_trim \
	--outSAMtype BAM SortedByCoordinate --outStd SAM --outSAMattributes Standard --outSAMstrandField intronMotif --outSAMunmapped None \
	--outFilterType BySJout outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix wt_rep1
mv wt_rep1.Aligned.sortedByCoord.out.bam $wt_rep2.bam
samtools index ${wt_rep2}.bam
# K
STAR --runThreadN 12 --genomeDir $star_index --genomeLoad NoSharedMemory --readFilesIn $ko_rep1_r1_trim $ko_rep1_r2_trim \
	--outSAMtype BAM SortedByCoordinate --outStd SAM --outSAMattributes Standard --outSAMstrandField intronMotif --outSAMunmapped None \
	--outFilterType BySJout outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix ko_rep1
mv ko_rep1.Aligned.sortedByCoord.out.bam $ko_rep1.bam
samtools index ${ko_rep1}.bam
STAR --runThreadN 12 --genomeDir $star_index --genomeLoad NoSharedMemory --readFilesIn $ko_rep2_r1_trim $ko_rep2_r2_trim \
	--outSAMtype BAM SortedByCoordinate --outStd SAM --outSAMattributes Standard --outSAMstrandField intronMotif --outSAMunmapped None \
	--outFilterType BySJout outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix ko_rep1
mv ko_rep1.Aligned.sortedByCoord.out.bam $ko_rep2.bam
samtools index ${ko_rep2}.bam

## Use IGVtools to convert BAM to TDF format, which can be visualized using IGV
wt_rep1_tdf =wt_rep1.tdf
wt_rep2_tdf =wt_rep2.tdf
ko_rep1_tdf =ko_rep1.tdf
ko_rep2_tdf =ko_rep2.tdf
igvtools count $wt_rep1_bam  $wt_rep1_tdf  $chrom_size_file
igvtools count $wt_rep2_bam  $wt_rep2_tdf  $chrom_size_file
igvtools count $ko_rep1_bam  $ko_rep1_tdf  $chrom_size_file
igvtools count $ko_rep2_bam  $ko_rep2_tdf  $chrom_size_file

## Use bamCoverage to convert BAM to bigWig format, which can be visualized with UCSC genome browser
wt_rep1_bw =wt_rep1.bw
wt_rep2_bw =wt_rep2.bw
ko_rep1_bw =ko_rep1.bw
ko_rep2_bw =ko_rep2.bw
bamCoverage -b $wt_rep1_bam  -o $wt_rep1_bw  --normalizeUsing RPKM
bamCoverage -b $wt_rep2_bam  -o $wt_rep2_bw  --normalizeUsing RPKM
bamCoverage -b $ko_rep1_bam  -o $ko_rep1_bw  --normalizeUsing RPKM
bamCoverage -b $ko_rep2_bam  -o $ko_rep2_bw  --normalizeUsing RPKM

## Get gene-level read counts using featureCount, which can be used for differential expression analysis and visualization
featureCounts -T 4 -p -a $gene_gtf_file -o MN_read_count.txt -s 2 *.bam

## Get TPM values for each gene using RSEM
rsem-calculate-expression --bowtie2 --paired-end -p 8 --strandedness reverse \
	$wt_rep1_r1_trim $wt_rep1_r2_trim $rsem_index
rsem-calculate-expression --bowtie2 --paired-end -p 8 --strandedness reverse \
	$wt_rep2_r1_trim $wt_rep2_r2_trim $rsem_index
rsem-calculate-expression --bowtie2 --paired-end -p 8 --strandedness reverse \
	$ko_rep1_r1_trim $ko_rep1_r2_trim $rsem_index
rsem-calculate-expression --bowtie2 --paired-end -p 8 --strandedness reverse \
	$ko_rep2_r1_trim $ko_rep2_r2_trim $rsem_index
	