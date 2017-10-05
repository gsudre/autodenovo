# Creates EGT file for PennCNV wave adjustment
# Takes in the PFB file, a genome reference, and the output file name as arguments

pfb=$1
ref=$2
out_fname=$3

module load R
module load bedtools

echo Padding PFB positions
Rscript ~/autodenovo/penncnv_pad500kb.R $pfb

echo Calculating GC content
bedtools nuc -fi $ref -bed for_gc.bed > gc_content.txt

echo Reformating GC model file
Rscript ~/autodenovo/penncnv_format_gc.R $out_fname

rm gc_content.txt for_gc.bed
