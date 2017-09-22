#~/bin/bash
# Evaluate the rest of XHMM< after running DepthOfCoverage
exome_targets='/data/NCR_SBRB/fake_trios/SeqCapEZ_Exome_v3.0_Design_Annotation_files/SeqCap_EZ_Exome_v3_hg19_capture_targets.bed' 
gatk_memory="50g"
ref_fa='/fdb/GATK_resource_bundle/hg19-2.8/ucsc.hg19.fasta'
out_dir='/data/NCR_SBRB/big_fake_simplex/xhmm'

cd $out_dir
module load GATK
module load XHMM

GATK -m ${gatk_memory} GCContentByInterval -L ${exome_targets} -R ${ref_fa} -o ./DATA.locus_GC.txt

cat ./DATA.locus_GC.txt | awk '{if ($2 < 0.1 || $2 > 0.9) print $1}' > ./extreme_gc_targets.txt

# merging all subjects in the directory
ls -1 *.sample_interval_summary > depth_list.txt;
xhmm --mergeGATKdepths --GATKdepthsList=depth_list.txt -o ./DATA.RD.txt; 

/usr/local/apps/XHMM/2016-01-04/sources/scripts/interval_list_to_pseq_reg ${exome_targets} > ./EXOME.targets.reg

xhmm --matrix -r ./DATA.RD.txt --centerData --centerType target \
-o ./DATA.filtered_centered.RD.txt \
--outputExcludedTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt \
--outputExcludedSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt \
--minTargetSize 10 --maxTargetSize 10000 \
--minMeanTargetRD 10 --maxMeanTargetRD 500 \
--minMeanSampleRD 25 --maxMeanSampleRD 200 \
--maxSdSampleRD 150

xhmm --PCA -r ./DATA.filtered_centered.RD.txt --PCAfiles ./DATA.RD_PCA

xhmm --normalize -r ./DATA.filtered_centered.RD.txt --PCAfiles ./DATA.RD_PCA \
--normalizeOutput ./DATA.PCA_normalized.txt \
--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7

xhmm --matrix -r ./DATA.PCA_normalized.txt --centerData --centerType sample --zScoreData \
-o ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
--outputExcludedTargets ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--outputExcludedSamples ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
--maxSdTargetRD 30

xhmm --matrix -r ./DATA.RD.txt \
--excludeTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt \
--excludeTargets ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--excludeSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt \
--excludeSamples ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
-o ./DATA.same_filtered.RD.txt

xhmm --discover -p params.txt -r ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R ./DATA.same_filtered.RD.txt -c ./DATA.xcnv -a ./DATA.aux_xcnv -s ./DATA

