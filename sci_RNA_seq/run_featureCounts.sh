#! /usr/bin/bash
gtf="/root/wangje/Reference/Homo_sapiens/UCSC/hg38/Annotation/Genes/hg38.gtf"
thread=20
featureCounts -a $gtf \
-o gene_assigned \
-R BAM test.Aligned.sortedByCoord.out.bam \
-T $thread
