#~/bin/bash
# Runs the full CNVkit pipeline
id=$1
w=$2
exome_targets='/data/NCR_SBRB/simplex/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed' 
ref_fa='/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
out_dir='/data/NCR_SBRB/simplex/cnvkit'
tmp_dir="/lscratch/$SLURM_JOBID"
bam_dir='/data/NCR_SBRB/simplex/BAM'
bam_name="${id}_sorted_RG_markduplicate_recalibrated"
ncpu=6

# let's use a temp dir so we don't keep overwriting the target and antitarget beds
cd $tmp_dir
module load cnvkit

cnvkit.py batch ${bam_dir}/${id}/${bam_name}.bam -n -t ${exome_targets} \
	-f ${ref_fa} --access ${out_dir}/access-${w}kb.hg19.bed --output-reference my_flat_reference.cnn -d ./ -p $ncpu;

cnvkit.py call ./${bam_name}.cns -y -m threshold -t=-1.1,-0.4,0.3,0.7 -o ./${bam_name}.call.cns;

cnvkit.py export vcf ./${bam_name}.call.cns -i "${id}" -o ./${bam_name}.call.cnv.vcf;

mkdir $out_dir/output_${w}kb
mv ${id}* $out_dir/output_${w}kb/
