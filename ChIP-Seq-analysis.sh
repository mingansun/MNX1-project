#!/bin/sh
## Ming-an Sun, 2021-08-11

###############################################################################
# Pipelie for routine analysis and visualization of ChIP-Seq data
###############################################################################

## set input file information
# ChIP file
chip_file=chip.fq
# Input file
input_file=input.fq
# Reference genome
bowtie2_index="~/db/Bowtie2Index/mm10"
chrom_size_file="~/db/Genome/mm10.chrom.sizes"
# ENCODE blacklist file v2
blacklist_file="~/db/Blacklist/mm10-blacklist.v2.bed"

## Trim fastq file using Trim Galore
chip_trim=chip_trimmed.fq
input_trim=input_trimmed.fq
trim_galore --illumina $chip_file
trim_galore --illumina $input_file

## Perform alignment, filtering, sorting and PCR duplicate removal
chip_bam=chip.bam
input_bam=input.bam
# ChIP
bowtie2 -p 4 -x $bowtie2_index $chip_trim | \
samtools view -b -F 4 - | \
samtools sort -@ 4 - | \
samtools rmdup -s - $chip_bam
samtools index $chip_bam
# Input
bowtie2 -p 4 -x $bowtie2_index $input_trim | \
samtools view -b -F 4 - | \
samtools sort -@ 4 - | \
samtools rmdup -s - $input_bam
samtools index $input_bam

## Use IGVtools to convert BAM to TDF format, which can be visualized using IGV
chip_tdf=chip.tdf
input_tdf=input.tdf
igvtools count $chip_bam  $chip_tdf  $chrom_size_file
igvtools count $input_bam $input_tdf $chrom_size_file

## Use bamCoverage to convert BAM to bigWig format, which can be visualized with DeepTools or UCSC genome browser
chip_bw=chip.bw
input_bw=input.bw
bamCoverage -b $chip_bam  -o $chip_bw  --normalizeUsing RPKM
bamCoverage -b $input_bam -o $input_bw --normalizeUsing RPKM

## Perform peak calling with MACS2 and filter out blacklist regions
peak_prefix=chip
peak_file=${peak_prefix}_peaks.narrowPeak
peak_file_filtered=${peak_prefix}.bed
macs2 callpeak -c $input_bam -t $chip_bam -n chip -g mm --keep-dup all -q 0.05
windowBed -a $peak_file -b $blacklist_file -w 0  -v >$peak_file_filt

## Visualize ChIP signal flanking given regions as heatmap using DeepTools
matrix_file=heatmap.matrix
heatmap_file=heatmap.png
computeMatrix reference-point --referencePoint center 
	-b 5000 -a 5000 \
    -R $peak_file_filtered \
    -S $chip_bam $input_bam \
    --outFileName $matrix_file
plotHeatmap -m $matrix_file -out $heatmap_file \
	--plotFileFormat png \
	--heatmapWidth 3 --heatmapHeight 18 \
	--refPointLabel 0 --xAxisLabel 'Center' 
	--missingDataColor 0.9 --colorMap Purples Greys
