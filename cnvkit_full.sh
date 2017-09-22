#~/bin/bash
# Runs the full CNVkit pipeline
id=$1
w=$2
exome_targets='/data/NCR_SBRB/fake_trios/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed' 
ref_fa='/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
out_dir='/data/NCR_SBRB/big_fake_simplex/cnvkit'
bam_dir='/data/NCR_SBRB/big_fake_simplex/BAM'
bam_name="${id}_sorted_RG_markduplicate_recalibrated"
ncpu=6

cd $out_dir
module load cnvkit

cnvkit.py batch ${bam_dir}/${id}/${bam_name}.bam -n -t ${exome_targets} \
	-f ${ref_fa} --access access-${w}kb.hg19.bed --output-reference my_flat_reference.cnn -d output_${w}kb/ -p $ncpu;

cnvkit.py call output_${w}kb/${bam_name}.cns -y -m threshold -t=-1.1,-0.4,0.3,0.7 -o output_${w}kb/${bam_name}.call.cns;

cnvkit.py export vcf output_${w}kb/${bam_name}.call.cns -i "${id}" -o output_${w}kb/${bam_name}.call.cnv.vcf;
