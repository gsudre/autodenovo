#~/bin/bash
# Runs GATK Depth of coverage, first step in XHMM
id=$1
exome_targets='/data/NCR_SBRB/simplex/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed' 
gatk_memory="50g"
ref_fa='/fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta'
bam_dir='/data/NCR_SBRB/simplex/BAM'
out_dir='/data/NCR_SBRB/simplex/xhmm'

cd $out_dir
module load GATK
GATK -m ${gatk_memory} DepthOfCoverage \
-I  ${bam_dir}/${id}/${id}_sorted_RG_markduplicate_recalibrated.bam \
-L ${exome_targets} \
-R ${ref_fa} \
-dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable \
--minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 \
--includeRefNSites \
--countType COUNT_FRAGMENTS \
-o ${id}.DOC
