**09/07/2017**

As I write this, the GATK best practices pipeline still doesn't have recommendations for CNV callers. So, while I'm adopting that for data processing and calling SNPs/Indels (mostly through biobc-nextgen), I'll use a group of other callers to get CNVs. Because the data I'll primarily be analyzing is WES, I'll focus on the tools that are specialized for that:

* XHMM
* EXCAVATORtool
* Conifer (I'm having issues with this code, looks too old for its dependencies)
* cn.mops?
* CNVnator
* ExomeCNV
* ExomeCopy
* CNVkit
* ExomeDepth: calls CNV from WES BAMs using RD. Reading a bit more about this (http://www.cureffi.org/2014/01/17/comparison-of-tools-for-calling-cnvs-from-sequence-data/), it looks like it needs to create its reference from the samples, which should be unrelated and not contain the “target” CNV. So, we could potentially do it in all parents of the trios? Just not sure if we’d have enough people. Worth a try.

There was also CONTRA, but it looked like it was more focuse don the normal/tumor cases.

This review of tools that find CNV from exomes (http://onlinelibrary.wiley.com/doi/10.1002/humu.22537/full) looked at XHMM, CoNIFER, ExomeDepth, and CONTRA. See above for my impressions on them. If I find many variants using just the parents in the family as the template (or all parents together), I can use the SNP data as well. It can be used either to add more candidate variants to the pool, or to filter the results from the WES analysis.

Even though most of them use read depth to infer CNVs, they have different algorithms. As expected, each of their papers show how much better one tool is when compared to the other. XHMM compares it to CONIFER, and finds XHMM better. The XHMM paper is actually quite interesttring, and provides an application example using SCZ trios (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3484655/, and the supplem,ental material is here http://ars.els-cdn.com/content/image/1-s2.0-S000292971200417X-mmc1.pdf). 

Note that this paper (https://molecularcytogenetics.biomedcentral.com/articles/10.1186/s13039-017-0333-5) that compares 3 read-depth tools to find CNV in WES looked at XHMM, CONIFER and CNVNator (modified). XHMM and CONIFER were somewhat similar, but different than CNVNator. Maybe worth taking a look at CNVnator as well? Need to look at the help page, or the paper itself, to check how CNVnator needs to be modified to run WES. I also found http://cnvkit.readthedocs.io/en/stable/, which seems to be a quite new tool based on its publication date (2016, http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873). They also include copywriteR and CONTRA in their comparisons. Note that bcbio also implements cnvkit, so maybe worth doing it through there?

The excavator2 paper (https://academic.oup.com/nar/article/44/20/e154/2607979/Enhanced-copy-number-variants-detection-from-whole?keytype=ref&ijkey=O8r64Qj81gfMzLo) compares it to CODEX, CONIFER, COpywriteR, and XHMM. 

So, because of the different mentions, let's stick with these for now:

* XHMM
* Excavator
* CNVKit
* CNVnator
* CODEX
* CopywriteR

So, the goal for CNVs from WES is to use SURVIVOR to combine the calls made by the different tools.

On the SNP chip front, I plan to stick with PennCNV, which has its own algorithm for denovo calling. I might try to apply it to the variants I find in the other arms as well. Scan2CNV looks like a wrapper to it, but I'm not sure if it implements all functionalities I need.

For SNPs/INDELS, I'll use multiple callers implemented in bcbio-gen, but I'm not sure yet whether it's a good idea to call them individually, or do some sort of joint/population call. If the latter, should we do it within families, or across all samples?

We can use finally trioDenovo or polymutt in the bcbio-nextgen results.

The final goal is to annotate the variants obtained in the 3 arms, and compare them to public (SV?) databases (dbVar, ClinVar).
