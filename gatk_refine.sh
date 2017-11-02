#!/bin/bash
#SBATCH --job-name=gatk1
###SBATCH --mail-type=ALL
#BROAD Best Practice

#####################
#Parameters to check#
#####################
id=$1
gatk_memory="50g"
ref_fa="/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
hg19_known_snps="/fdb/GATK_resource_bundle/b37/1000G_phase3_v4_20130502.sites.vcf.gz"
home_directory="/data/NCR_SBRB/simplex"
vcf_directory="${home_directory}/VCF"
my_gc_threads=$(($SLURM_CPUS_PER_TASK-1))

################
#Upload modules#
################
module load GATK/3.6

trio_fname=${home_directory}/${1}.ped
trio_name=$1

cd $vcf_directory
# Step 1: Derive posterior probabilities of genotypes
GATK -m ${gatk_memory} CalculateGenotypePosteriors -R ${ref_fa} --supporting ${hg19_known_snps} -ped ${trio_fname} -pedValidationType SILENT -V recalibrated_variants.vcf -o recalibrated_variants.${trio_name}.vcf 

# Step 2: Filter low quality genotypes
GATK -m ${gatk_memory} VariantFiltration -R ${ref_fa} -V recalibrated_variants.${trio_name}.vcf -G_filter "GQ < 20.0" -G_filterName lowGQ -o recalibrated_variants.${trio_name}.Gfiltered.vcf

# Step 3: Annotate possible de novo mutations
GATK -m ${gatk_memory} VariantAnnotator -R ${ref_fa} -V recalibrated_variants.${trio_name}.Gfiltered.vcf -A PossibleDeNovo -ped ${trio_fname} -pedValidationType SILENT -o recalibrated_variants.${trio_name}.Gfiltered.deNovos.vcf

# Finally filter the high confidence calls
mkdir ../gatk_refine
grep '#' recalibrated_variants.${trio_name}.Gfiltered.deNovos.vcf > ../gatk_refine/${trio_name}_hiConfDeNovo.vcf
grep hiConfDeNovo recalibrated_variants.${trio_name}.Gfiltered.deNovos.vcf >> ../gatk_refine/${trio_name}_hiConfDeNovo.vcf
grep '#' recalibrated_variants.${trio_name}.Gfiltered.deNovos.vcf > ../gatk_refine/${trio_name}_allDeNovo.vcf
grep DeNovo recalibrated_variants.${trio_name}.Gfiltered.deNovos.vcf >> ../gatk_refine/${trio_name}_allDeNovo.vcf
