thread=20
genomeDir="/root/wangje/Reference/Star_index/UCSC"
input_file="/root/wangje/Project/Yin/Data/extract_umi_barcode/I230303_I230303-11_2.barcode.clearn.fastq.gz"
output_file_prefix="I230303_I230303-11_2.barcode.clearn.star"

STAR --runThreadN $thread \
--genomeDir $genomeDir \
--readFilesIn $input_file \
--readFilesCommand zcat \
#比对时允许的最大错配数
--outFilterMultimapNmax 1 \ 
--outFileNamePrefix $output_file_prefix \
--outSAMtype BAM SortedByCoordinate
