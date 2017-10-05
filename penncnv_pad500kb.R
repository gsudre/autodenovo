# Creates a .bed file for later GC region extraction.
# Takes a pfb file name as its only argument

args <- commandArgs(TRUE)
pos = read.table(args[1])
start = pos[,3]-500000
end = pos[,3]+500000
start[start<0] = 0
chrs = vector()
for (i in pos[2]) {
    chrs = c(chrs, sprintf('chr%s', i))
}
df=data.frame(chr=chrs,start=start,end=end, rsid=pos[,1], pos=pos[,3])  # might need rsid and pos later for merging

# remove non-autossomal chromosomes
keep_me = df$chr!='chrMT' & df$chr!='chrXY' & df$chr!='chrX' & df$chr!='chrY' & df$chr!='chr0'

df = df[keep_me,]
# compromising based on chromosome lengths based on error outputs from bedtools
df[df$chr=='chr19' & df$end > 59128983,]$end = 59128983
df[df$chr=='chr16' & df$end > 90354753,]$end = 90354753
df[df$chr=='chr5' & df$end > 180915260,]$end = 180915260
df[df$chr=='chr22' & df$end > 51304566,]$end = 51304566
df[df$chr=='chr18' & df$end > 78077248,]$end = 78077248
df[df$chr=='chr8' & df$end > 146364022,]$end = 146364022
df[df$chr=='chr12' & df$end > 133851895,]$end = 133851895
df[df$chr=='chr13' & df$end > 115169878,]$end = 115169878
df[df$chr=='chr15' & df$end > 102531392,]$end = 102531392
df[df$chr=='chr17' & df$end > 81195210,]$end = 81195210
df[df$chr=='chr20' & df$end > 63025520,]$end = 63025520
df[df$chr=='chr21' & df$end > 48129895,]$end = 48129895
df[df$chr=='chr10' & df$end > 135534747,]$end = 135534747
df[df$chr=='chr1' & df$end > 249250621,]$end = 249250621
df[df$chr=='chr2' & df$end > 243199373,]$end = 243199373
df[df$chr=='chr3' & df$end > 198022430,]$end = 198022430
df[df$chr=='chr7' & df$end > 159138663,]$end = 159138663
df[df$chr=='chr4' & df$end > 191154276,]$end = 191154276
df[df$chr=='chr6' & df$end > 171115067,]$end = 171115067
df[df$chr=='chr9' & df$end > 141213431,]$end = 141213431
df[df$chr=='chr14' & df$end > 107349540,]$end = 107349540
df[df$chr=='chr11' & df$end > 135006516,]$end = 135006516
write.table(df, file='for_gc.bed',row.names=F,col.names=F,sep='\t',quote=F)
