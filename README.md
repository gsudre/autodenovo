# AutoDenovo

The goal in this project is to construct a pipeline for finding denovo mutations related to neurodevelopmental disorders using family study designs, similar to pipelines previously developed for cancer research. It combines different types of variants by multiple callers, and uses both chip and exome sequencing data. 

The end goal is to compare annotated variants to public databases (dbVar, ClinVar), possibly indicating further sequencing or mechanistic avenues.

This is an ongoing project, and I've been journaling about my experiences in this project's GitHub [wiki](https://github.com/gsudre/autodenovo/wiki/AutoDenovo-Wiki).

## Graphical Overview
![Flowchart](main_flow.png "MainFlow")

## Description

The current implementation needs the following files:

* Whole Exome Sequencing (WES) data for each participant
* Trio pedigree files (.ped), one trio per file
* SNP chip data for each participant

Note that some tools (e.g. XHMM, GATK joint calling) make use of all samples at the same time, regardless of their family structure. Also note that SNP chip data is not required, but you'll only be able to run 2 of the 3 arms in the pipeline. This approach of consensus calling and multiple sources of information tries to reduce false positives and increase a clinician's confidence on the results provided. 

The end result of each analysis arm is a file containing de-novo mutations for each trio. The options for analysis using each result are many. For example:

* If some trios contain affected individuals and other don't, one option for filtering variants is to include only the ones that appear in most (or all) affected trios. Then, one can check how often those de-novo mutations appeared in the unnafected trios.
* If affected and unnafected trios come from the same family, this variant comparison can be performed within family first, before checking across families.
* The parameters specific to each caller can be tweaked to increase/decrease sensitivity and specificity, and hence obtain more/less variants in the final consensus call of each arm.
* Compare annotated variants to public databases, possibly indicating further sequencing or mechanistic avenues.

## Usage

The SNP Chip pipeline can be run in paralle with the other two, but SNV and CNV need the files from GATK Best Practices pipeline. So, we'll start by running that.

### CNV and SNV arms

#### 1. Run GATK pipeline script to do variant calling for single samples

First, modify the script gatk_upToSingleCalls.sh to reflect your environment. This script runs the GATK pipeline up to calls for the single sample. It takes the sample name, assuming all variables are properly set. I highly recommend running it in a computer cluster. In my case, we have a cluster running SLURM, so I swarm the jobs this way:

```bash
while read s; do echo "bash ~/autodenovo/gatk_upToSingleCalls.sh $s" >> swarm.single; done < sample_ids.txt
swarm -f swarm.single -t 16 -g 55 --job-name gatk_single --logdir trash --time=48:00:00 --gres=lscratch:100
```

where sample_ids.txt is a file with one sample name per line. This should take about 1-2 days to run per participant, so ther'es the beauty of running them all in parallel.

### SNP Chip arm

As the SNP chip arm doesn't depend on the exome data, let's run it while we wait for the GATK pipeline results. 

#### 1. Prepare files

I __strongly__ advise you to look at the [PennCNV documentation](http://penncnv.openbioinformatics.org/en/latest/user-guide/input/) to figure out how to prepare your input files. 

In my case, we used Illumina's InfiniumExome-24_v1.0 chips. So, I used GenomeStudio (free to [download](https://support.illumina.com/array/array_software/genomestudio/downloads.html), but it does require Windows) to export a Final Report with the columns we need `(SNP Name, Sample ID, Allele1 - Top, Allele2 - Top, GC Score, Log R Ratio, B Allele Freq.)`. 

We also have a few samples that were genotyped on Illumina's HumanExome-12v1-2. Louckily they are all in the same family, otherwise you'd have to come up with files representing the intersection of the array. Contact me if you'd like instructions on how to do that, or check the Wiki, because I had to do it for the sample data I was using there.

For now, let's go ahead and split the data:

```bash
mkdir penncnv
cd penncnv
module load penncnv
mkdir HumanExome
split_illumina_report.pl -prefix HumanExome/ fs_ccgo_box4_FinalReport.txt 
mkdir InfiniumExome
split_illumina_report.pl -prefix InfiniumExome/ newBox225_FinalReport.txt 
split_illumina_report.pl -prefix InfiniumExome/ newBox226_FinalReport.txt 
```

We then have to download PFB files from [Illumina's website](https://support.illumina.com/array/downloads.html). They don't have those exact files there, but all we're looking for is the Population Frequency of B Allele, which is calculated in their MAF files. We just need to extract the column that corresponds to our population. In my case, I'll get [this](ftp://webdata:webdata@ussd-ftp.illumina.com/Downloads/ProductFiles/HumanExome/ProductSupportFiles/HumanExome-12v1-2_A_MAF.txt) and [this](https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/infinium-exome-24/infinium-exome-24-v1-0-a1-population-reports-maf-copy-numbers.zip). So, to make our PFB files, we need:

```bash
# get the first 4 columns, and ignore the first line
cut -f 1,2,3,4 HumanExome-12v1-2_A_MAF.txt | tail -n +2 > HumanExome.pfb
cut -f 1,2,3,4 Population\ Reports\ \(MAF\,\ Copy\ Numbers\)/InfiniumExome-24v1-0_A1_PopulationReport_MAF.txt | tail -n +2 > InfiniumExome.pfb
```

I'm going to use the default HMM parameters, but feel free to play with them if you want to increase/decrease sensitivity and specificity. Finally, I also created a GC content file to test the corrections implemented in PennCNV. There's a simple script to do it, but you'll likely need to change it to fit your needs. Also, __this takes a few hours__ to run, so use it wisely.

```bash
bash ~/autodenovo/penncnv_create_gcmodel.sh HumanExome.pfb /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa HumanExome.h19.gcmodel
bash ~/autodenovo/penncnv_create_gcmodel.sh InfiniumExome.pfb /fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa InfiniumExome.h19.gcmodel
```
