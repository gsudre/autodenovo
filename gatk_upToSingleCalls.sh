#!/bin/bash
#SBATCH --job-name=gatkv5
###SBATCH --mail-type=ALL
#BROAD Best Practice

#####################
#Parameters to check#
#####################
# sbatch --cpus-per-task=32 --mem=120g --time=10-00:00:00 /data/NCR_SBRB/ADHDQTL/SCRIPTS/BASH/gatk_3.8_best_practice_step1_v6.7_CLIA_CCGO.script
# id="XXX"
id=$1
#bwa_fa_index="/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"
bwa_fa_index="/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"
samtools_sort_memory="118G"
java_memory="-Xmx118g"
gatk_memory="118g"
my_gc_threads=$(($SLURM_CPUS_PER_TASK-1))
#ref_fa="/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
ref_fa="/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
hg38_known_snps="/fdb/GATK_resource_bundle/hg38/dbsnp_146.hg38.vcf.gz"
hg38_known_indels="/fdb/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
hg19_known_snps="/fdb/GATK_resource_bundle/hg19-2.8/dbsnp_138.hg19.vcf.gz"
hg19_known_indels="/fdb/GATK_resource_bundle/hg19-2.8/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz"
na12878_known_sites="/data/yuns/common/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf"
GIAB_SNPs="/data/yuns/common/GIAB_GRCh38_SNPs.vcf"
GIAB_INDELs="/data/yuns/common/GIAB_GRCh38_INDELs.vcf"
my_1000G_hg38_vcf="/data/yuns/common/1000G_phase1.snps.high_confidence.hg38.vcf"
my_dbsnp_146_hg38_vcf="/data/yuns/common/dbsnp_146.hg38.vcf"
my_TruSeq3_PE_2_fa="/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa"
my_TruSeq2_PE_fa="/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq2-PE.fa"

hg19_hapmap_vcf="/fdb/GATK_resource_bundle/hg19-2.8/hapmap_3.3.hg19.vcf.gz"
hg19_1000G_omni_vcf="/fdb/GATK_resource_bundle/hg19-2.8/1000G_omni2.5.hg19.vcf.gz"
hg19_1000G_phase1_snps="/fdb/GATK_resource_bundle/hg19-2.8/1000G_phase1.snps.high_confidence.hg19.vcf.gz"
hg19_dbsnp_138="/fdb/GATK_resource_bundle/hg19-2.8/dbsnp_138.hg19.vcf.gz"
hg19_mills_indels_vcf="/fdb/GATK_resource_bundle/hg19-2.8/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz"

home_directory="/data/NCR_SBRB/ADHDQTL/GATK"
fastq_directory="${home_directory}/FASTQ"
nisc_bam_directory="/data/NCR_SBRB/ADHDQTL/NISC/SampleBams"
qc_top_directory="${home_directory}/QC"
qc_this_sample_directory="${qc_top_directory}/${id}"
qc_this_sample_R1_directory="${qc_this_sample_directory}/R1"
qc_this_sample_R2_directory="${qc_this_sample_directory}/R2"
bam_directory="${home_directory}/BAM"
bam_this_sample_directory="${bam_directory}/${id}"
vcf_directory="${home_directory}/VCF"
vcf_this_sample_directory="${vcf_directory}/${id}"
annovar_directory="${home_directory}/ANNOVAR"

java_memory="-Xmx118g"
nisc_bam="${nisc_bam_directory}/${id}.bam"
fastq_R1="${fastq_directory}/${id}_R1.fastq.gz"
fastq_R2="${fastq_directory}/${id}_R2.fastq.gz"

################
#Upload modules#
################
module load picard/2.9.2
module load fastqc/0.11.5
module load bwa/0.7.15
module load samtools/1.6
module load GATK/3.8-0
module load annovar/2017-07-16
module load vcftools/0.1.15
module load bbtools/37.36

#:<<COMMENT
###################
###Getting FASTQ###
###################
cd ${home_directory}
mkdir -p ${fastq_directory}
cd ${fastq_directory}
java ${java_memory} -jar $PICARDJARPATH/picard.jar SamToFastq I=${nisc_bam} FASTQ=${fastq_R1} SECOND_END_FASTQ=${fastq_R2}
cd ..

###########
###bbduk###
###########
# bbtools bbduk threads=1 in=${fastq_R1} in2=${fastq_R2} out=${fastq_directory}/${id}_bbduk_R1.fastq.gz out2=${fastq_directory}/${id}_bbduk_R2.fastq.gz ftl=16 ftr=90 overwrite=true

bbtools bbduk threads=1 in=${fastq_R1} in2=${fastq_R2} out=${fastq_directory}/${id}_tbo_R1.fastq.gz out2=${fastq_directory}/${id}_tbo_R2.fastq.gz tbo overwrite=true

bbtools bbduk threads=1 in=${fastq_directory}/${id}_tbo_R1.fastq.gz in2=${fastq_directory}/${id}_tbo_R2.fastq.gz out=${fastq_directory}/${id}_bbduk_R1.fastq.gz out2=${fastq_directory}/${id}_bbduk_R2.fastq.gz ref=/data/NCR_SBRB/ADHDQTL/SCRIPTS/BASH/adapters.u.i.fa k=21 hdist=1 qtrim=rl trimq=19 minlen=30 overwrite=true

############
###FASTQC###
############
cd ${home_directory}
mkdir -p ${qc_top_directory}
mkdir -p ${qc_this_sample_directory}
mkdir -p ${qc_this_sample_R1_directory}
mkdir -p ${qc_this_sample_R2_directory}
fastqc -o ${qc_this_sample_R1_directory} -f fastq ${fastq_R1}
fastqc -o ${qc_this_sample_R2_directory} -f fastq ${fastq_R2}
fastqc -o ${qc_this_sample_R1_directory} -f fastq ${fastq_directory}/${id}_bbduk_R1.fastq.gz
fastqc -o ${qc_this_sample_R2_directory} -f fastq ${fastq_directory}/${id}_bbduk_R2.fastq.gz

#########
###BWA###
#########
mkdir -p ${bam_directory}
mkdir -p ${bam_this_sample_directory}
cd ${bam_this_sample_directory}
bwa mem -M -t $SLURM_CPUS_PER_TASK  ${bwa_fa_index} ${fastq_directory}/${id}_bbduk_R1.fastq.gz ${fastq_directory}/${id}_bbduk_R2.fastq.gz > ${id}_trimmed.sam

################
###Sorted BAM###
################
cd ${bam_this_sample_directory}
samtools sort --threads $SLURM_CPUS_PER_TASK -m ${samtools_sort_memory} -o ${id}_trimmed_sorted.bam ${id}_trimmed.sam
samtools index ${id}_trimmed_sorted.bam

#####################
###Add Read Groups###
#####################
cd ${bam_this_sample_directory}
java ${java_memory} -XX:ParallelGCThreads=${my_gc_threads} -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups I=${id}_trimmed_sorted.bam O=${id}_trimmed_sorted_RG.bam RGID=${id} RGPL=illumina RGSM=${id} RGPU=unit1 RGLB=${id}

#####################
###Mark Duplicates###
#####################
cd ${bam_this_sample_directory}
java ${java_memory} -XX:ParallelGCThreads=${my_gc_threads} -jar $PICARDJARPATH/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 I=${id}_trimmed_sorted_RG.bam O=${id}_trimmed_sorted_RG_markduplicate.bam METRICS_FILE=${id}_trimmed_sorted_RG_markduplicate.bam.metrics.txt
samtools index ${id}_trimmed_sorted_RG_markduplicate.bam

########################
###Base Recalibration###
########################
cd ${bam_this_sample_directory}
#1st pass: Analyze patterns of covariation in the sequence dataset
GATK -m ${gatk_memory} BaseRecalibrator -R ${ref_fa} -I ${id}_trimmed_sorted_RG_markduplicate.bam -o ${id}_trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_1st -knownSites ${hg19_known_snps} -knownSites ${hg19_known_indels}
#GATK ==> /usr/local/apps/GATK/4.0.0.0/gatk BaseRecalibrator
# /usr/local/apps/GATK/4.0.0.0/gatk --java-options ${java_memory} BaseRecalibrator -R ${ref_fa} -I ${id}_trimmed_sorted_RG_markduplicate.bam -O ${id}_trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_1st --known-sites ${hg19_known_snps} --known-sites ${hg19_known_indels}
#2nd pass: Analyze covariation remaining after recalibration
GATK -m ${gatk_memory} BaseRecalibrator -R ${ref_fa} -I ${id}_trimmed_sorted_RG_markduplicate.bam -o ${id}_trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_2nd -knownSites ${hg19_known_snps} -knownSites ${hg19_known_indels} -BQSR ${id}_trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_1st
#/usr/local/apps/GATK/4.0.0.0/gatk --java-options ${java_memory} BaseRecalibrator -R ${ref_fa} -I ${id}_trimmed_sorted_RG_markduplicate.bam -O ${id}_trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_2nd --known-sites ${hg19_known_snps} --known-sites ${hg19_known_indels} --bqsr-recal-file ${id}_trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_1st
# /usr/local/apps/GATK/4.0.0.0/gatk --java-options ${java_memory} ApplyBQSR -R ${ref_fa} -I ${id}_trimmed_sorted_RG_markduplicate.bam -O ${id}_trimmed_sorted_RG_markduplicate_recalibrated.bam --bqsr-recal-file ${id}_trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_1st
# _trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_2nd --known-sites ${hg19_known_snps} --known-sites ${hg19_known_indels}
#3rd step: Generate before & after plots
GATK -m ${gatk_memory} AnalyzeCovariates -R ${ref_fa} -before ${id}_trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_1st -after ${id}_trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_2nd -plots ${id}_trimmed_sorted_RG_markduplicate.bam.recalibration_plots.pdf
#/usr/local/apps/GATK/4.0.0.0/gatk --java-options ${java_memory} AnalyzeCovariates -R ${ref_fa} -before ${id}_trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_1st -after ${id}_trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_2nd -plots ${id}_trimmed_sorted_RG_markduplicate.bam.recalibration_plots.pdf
#4th step: Apply the recalibration to the sequence data
GATK -m ${gatk_memory} PrintReads -R ${ref_fa} -I ${id}_trimmed_sorted_RG_markduplicate.bam -o ${id}_trimmed_sorted_RG_markduplicate_recalibrated.bam -BQSR ${id}_trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_1st  #-BQSR will do the on-the-fly recalibration, which is equal to 2nd pass
#/usr/local/apps/GATK/4.0.0.0/gatk --java-options ${java_memory} PrintReads -R ${ref_fa} -I ${id}_trimmed_sorted_RG_markduplicate.bam -o ${id}_trimmed_sorted_RG_markduplicate_recalibrated.bam -BQSR ${id}_trimmed_sorted_RG_markduplicate.bam.base_recalibration_table_1st  #-BQSR will do the on-the-fly recalibration, which is equal to 2nd pass

# COMMENT
#####################################
###Call Variants to generate g.vcf###
#####################################
mkdir -p ${vcf_directory}
mkdir -p ${vcf_this_sample_directory}
cd ${vcf_this_sample_directory}
GATK -m ${gatk_memory} HaplotypeCaller -R ${ref_fa} --dbsnp ${hg19_dbsnp_138} -I ${bam_this_sample_directory}/${id}_trimmed_sorted_RG_markduplicate_recalibrated.bam -o ${id}_trimmed.rsid4.g.vcf -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --doNotRunPhysicalPhasing
#/usr/local/apps/GATK/4.0.0.0/gatk --java-options ${java_memory} HaplotypeCaller -R ${ref_fa} --dbsnp ${hg19_dbsnp_138} -I ${bam_this_sample_directory}/${id}_trimmed_sorted_RG_markduplicate_recalibrated.bam -O ${id}_trimmed.rsid4.g.vcf -ERC GVCF --do-not-run-physical-phasing true --dont-use-soft-clipped-bases true
# TMP_DATA="/data/N...." --variant_index_type LINEAR --variant_index_parameter 128000 --doNotRunPhysicalPhasing





