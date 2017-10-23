#!/bin/bash
#SBATCH --job-name=gatk2
###SBATCH --mail-type=ALL
#BROAD Best Practice

#####################
#Parameters to check#
#####################
bam_list=$1
gatk_memory="118g"
my_gc_threads=$(($SLURM_CPUS_PER_TASK-1))
ref_fa="/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
gatk_files="/fdb/GATK_resource_bundle/hg19"

home_directory="/data/NCR_SBRB/simplex"
vcf_directory="${home_directory}/VCF"
annovar_directory="${home_directory}/ANNOVAR"


module load GATK
module load annovar

######################
###Joint Genotyping###
######################
cd ${vcf_directory}
v_str=''
while read id; do
	v_str=$v_str" -V ${id}/${id}.g.vcf"
done < $bam_list
GATK -m ${gatk_memory} GenotypeGVCFs -R ${ref_fa} ${v_str} -o joint.vcf -nt $SLURM_CPUS_PER_TASK

################################################
###VQSR (Variant Quality Score Recalibration)###
################################################
cd ${vcf_directory}
#Step1: Build SNP recalibration model
GATK -m ${gatk_memory} VariantRecalibrator -R ${ref_fa} -input joint.vcf \
	-recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches \
	-rscriptFile recalibrate_SNP_plots.R -nt $SLURM_CPUS_PER_TASK \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${gatk_files}/hapmap_3.3.hg19.vcf.gz \
   	-resource:omni,known=false,training=true,truth=true,prior=12.0 ${gatk_files}/1000G_omni2.5.hg19.vcf.gz \
   	-resource:1000G,known=false,training=true,truth=false,prior=10.0 ${gatk_files}/1000G_phase1.snps.high_confidence.hg19.vcf.gz \
   	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${gatk_files}/dbsnp_138.hg19.vcf.gz \
   	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
   	-mode SNP \
	

#Step2: Apply the desired level of recalibration to the SNPs in the call set
GATK -m ${gatk_memory} ApplyRecalibration -R ${ref_fa} -input joint.vcf -mode SNP --ts_filter_level 99.5 \
	-recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -o recalibrated_snps_raw_indels.vcf

#Step3: Build INDEL recalibration model
GATK -m ${gatk_memory} VariantRecalibrator -R ${ref_fa} -input recalibrated_snps_raw_indels.vcf \
	-recalFile recalibrate_INDEL.recal -tranchesFile recalibrate_INDEL.tranches -rscriptFile recalibrate_INDEL_plots.R \
	--maxGaussians 4 \
   	-resource:mills,known=false,training=true,truth=true,prior=12.0 ${gatk_files}/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz \
   	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${gatk_files}/dbsnp_138.hg19.vcf.gz \
   	-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
   	-mode INDEL \

#Step4: Apply the desired level of recalibration to the SNPs in the call set
GATK -m ${gatk_memory} ApplyRecalibration -R ${ref_fa} -input recalibrated_snps_raw_indels.vcf \
	-mode INDEL --ts_filter_level 99.0 -recalFile recalibrate_INDEL.recal \
	-tranchesFile recalibrate_INDEL.tranches -o recalibrated_variants.vcf
cd ..


#############
###ANNOVAR###
#############
cd ${annovar_directory}
table_annovar.pl ${vcf_directory}/recalibrated_variants.vcf  $ANNOVAR_DATA/hg19 --tempdir ${annovar_directory}/temp --thread $SLURM_CPUS_PER_TASK --buildver hg19 --out big_fake_simplex --remove --protocol refGene,exac03,avsnp147,clinvar_20170130,esp6500siv2_all -operation g,f,f,f,f --nastring . --vcfinput
