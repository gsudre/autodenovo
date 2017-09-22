#!/bin/bash
#SBATCH --job-name=gatk1
###SBATCH --mail-type=ALL
#BROAD Best Practice

#####################
#Parameters to check#
#####################
id=$1
bwa_fa_index="/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
samtools_sort_memory="3G" # this is per thread!
java_memory="-Xmx50g"
gatk_memory="50g"
my_gc_threads=$(($SLURM_CPUS_PER_TASK-1))
ref_fa="/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
hg19_known_snps="/fdb/GATK_resource_bundle/hg19/dbsnp_138.hg19.vcf.gz"
hg19_known_indels="/fdb/GATK_resource_bundle/hg19/1000G_phase1.indels.hg19.vcf.gz"

home_directory="/data/NCR_SBRB/big_fake_simplex"
fastq_directory="${home_directory}/FASTQ"
nisc_bam_directory="/data/NCR_SBRB/NISC/SampleBams"
qc_top_directory="${home_directory}/QC"
qc_this_sample_directory="${qc_top_directory}/${id}"
qc_this_sample_R1_directory="${qc_this_sample_directory}/R1"
qc_this_sample_R2_directory="${qc_this_sample_directory}/R2"
bam_directory="${home_directory}/BAM"
bam_this_sample_directory="${bam_directory}/${id}"
vcf_directory="${home_directory}/VCF"
vcf_this_sample_directory="${vcf_directory}/${id}"

nisc_bam="${nisc_bam_directory}/${id}.bam"
fastq_R1="${fastq_directory}/${id}_R1.fastq"
fastq_R2="${fastq_directory}/${id}_R2.fastq"

################
#Upload modules#
################
module load picard/2.9.2
module load fastqc/0.11.5
module load bwa/0.7.15
module load samtools/1.5
module load GATK/3.6


###################
###Getting FASTQ###
###################
mkdir -p ${fastq_directory}
cd ${fastq_directory}
if [ ! -e ${fastq_R1} ]; then
	java ${java_memory} -jar $PICARDJARPATH/picard.jar SamToFastq I=${nisc_bam} FASTQ=${fastq_R1} SECOND_END_FASTQ=${fastq_R2}
fi
cd ..


############
###FASTQC###
############
cd ${home_directory}
mkdir -p ${qc_top_directory}
if [ ! -d ${qc_this_sample_directory} ]; then
	mkdir -p ${qc_this_sample_directory}
	mkdir -p ${qc_this_sample_R1_directory}
	mkdir -p ${qc_this_sample_R2_directory}
	fastqc -o ${qc_this_sample_R1_directory} -f fastq ${fastq_R1}
	fastqc -o ${qc_this_sample_R2_directory} -f fastq ${fastq_R2}
fi


#########
###BWA###
#########
mkdir -p ${bam_directory}
mkdir -p ${bam_this_sample_directory}
cd ${bam_this_sample_directory}
if [ ! -e ${id}.sam ]; then
	bwa mem -M -t $SLURM_CPUS_PER_TASK  ${bwa_fa_index} ${fastq_R1} ${fastq_R2} > ${id}.sam
fi


################
###Sorted BAM###
################
cd ${bam_this_sample_directory}
if [ ! -e ${id}_sorted.bam ]; then
	samtools sort --threads $SLURM_CPUS_PER_TASK -m ${samtools_sort_memory} -T /lscratch/$SLURM_JOB_ID -o ${id}_sorted.bam ${id}.sam
	samtools index ${id}_sorted.bam
fi


#####################
###Add Read Groups###
#####################
cd ${bam_this_sample_directory}
if [ ! -e ${id}_sorted_RG.bam ]; then
	java ${java_memory} -XX:ParallelGCThreads=${my_gc_threads} -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups I=${id}_sorted.bam O=${id}_sorted_RG.bam RGID=${id} RGPL=illumina RGSM=${id} RGPU=unit1 RGLB=${id}
fi

#####################
###Mark Duplicates###
#####################
cd ${bam_this_sample_directory}
if [ ! -e ${id}_sorted_RG_markduplicate.bai ]; then
	java ${java_memory} -XX:ParallelGCThreads=${my_gc_threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 I=${id}_sorted_RG.bam O=${id}_sorted_RG_markduplicate.bam METRICS_FILE=${id}_sorted_RG_markduplicate.bam.metrics.txt
	samtools index ${id}_sorted_RG_markduplicate.bam
fi


########################
###Base Recalibration###
########################
cd ${bam_this_sample_directory}
if [ ! -e ${id}_sorted_RG_markduplicate_recalibrated.bam ]; then 
	#1st pass: Analyze patterns of covariation in the sequence dataset
	GATK -m ${gatk_memory} BaseRecalibrator -R ${ref_fa} -I ${id}_sorted_RG_markduplicate.bam -o ${id}_sorted_RG_markduplicate.bam.base_recalibration_table_1st -knownSites ${hg19_known_snps} -knownSites ${hg19_known_indels}
	#2nd pass: Analyze covariation remaining after recalibration
	GATK -m ${gatk_memory} BaseRecalibrator -R ${ref_fa} -I ${id}_sorted_RG_markduplicate.bam -o ${id}_sorted_RG_markduplicate.bam.base_recalibration_table_2nd -knownSites ${hg19_known_snps} -knownSites ${hg19_known_indels} -BQSR ${id}_sorted_RG_markduplicate.bam.base_recalibration_table_1st 
	#3rd step: Generate before & after plots
	GATK -m ${gatk_memory} AnalyzeCovariates -R ${ref_fa} -before ${id}_sorted_RG_markduplicate.bam.base_recalibration_table_1st -after ${id}_sorted_RG_markduplicate.bam.base_recalibration_table_2nd -plots ${id}_sorted_RG_markduplicate.bam.recalibration_plots.pdf
	#4th step: Apply the recalibration to the sequence data
	GATK -m ${gatk_memory} PrintReads -R ${ref_fa} -I ${id}_sorted_RG_markduplicate.bam -o ${id}_sorted_RG_markduplicate_recalibrated.bam -BQSR ${id}_sorted_RG_markduplicate.bam.base_recalibration_table_1st  #-BQSR will do the on-the-fly recalibration, which is equal to 2nd pass
fi

#####################################
###Call Variants to generate g.vcf###
#####################################
mkdir -p ${vcf_directory}
mkdir -p ${vcf_this_sample_directory}
cd ${vcf_this_sample_directory}
if [ ! -e ${id}.g.vcf.idx ]; then
	GATK -m ${gatk_memory} HaplotypeCaller -R ${ref_fa} -I ${bam_this_sample_directory}/${id}_sorted_RG_markduplicate_recalibrated.bam -o ${id}.g.vcf -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000
fi

