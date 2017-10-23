#!/bin/bash
# Runs a single subject/window pair in CNVnator

id=$1
window=$2
out_dir='/data/NCR_SBRB/simplex/cnvnator'
bam_file="/data/NCR_SBRB/simplex/BAM/${id}/${id}_sorted_RG_markduplicate_recalibrated.bam"
ref_dir='/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes'

cd /lscratch/$SLURM_JOBID
for i in {1..22} X Y; do
   ln -s ${ref_dir}/chr${i}.fa chr${i}.fa;
done

# make chromossomal string
my_chr='' ; for c in {1..22} X Y; do my_chr=${my_chr}' chr'$c; done

module load cnvnator

cnvnator -unique -root ${id}.root -chrom $my_chr -tree $bam_file;
cnvnator -root ${id}.root -chrom $my_chr -his $window;
cnvnator -root ${id}.root -chrom $my_chr -stat $window;
cnvnator -root ${id}.root -chrom $my_chr -partition $window;
cnvnator -root ${id}.root -chrom $my_chr -call $window > ${id}_calls_w${window}.txt;
cnvnator2VCF.pl ${id}_calls_w${window}.txt > ${id}_calls_w${window}.vcf;

mv *_calls_* ${out_dir}/
