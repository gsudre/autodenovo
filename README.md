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

### 1. Modify the script gatk_upToSingleCalls.sh to reflect your environment. 

This script runs the GATK pipeline up to calls for the single sample. It takes the sample name, assuming all variables are properly set. I highly recommend running it in a computer cluster. In my case, we have a cluster running SLURM, so I swarm the jobs this way:

```bash
while read s; do echo "bash ~/autodenovo/gatk_upToSingleCalls.sh $s" >> swarm.single; done < sample_ids.txt
swarm -f swarm.single -t 16 -g 55 --job-name gatk_single --logdir trash --time=48:00:00 --gres=lscratch:100
```

where sample_ids.txt is a file with one sample name per line. This should take about 1-2 days to run per participant, so ther'es the beauty of running them all in parallel.

