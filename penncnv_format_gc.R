# Reads in a gc_content file from bedtools and spits out a GCModel file for PennCNV

args <- commandArgs(TRUE)
gc = read.table('gc_content.txt')
df = gc[,c(4,1,5,7)]
df[,4] = df[,4]*100
write.table(df, file=args[1],row.names=F,col.names=F,sep='\t',quote=F)
