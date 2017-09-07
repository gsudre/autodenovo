**09/07/2017**

As I write this, the GATK best practices pipeline still doesn't have recommendations for CNV callers. So, while I'm adopting that for data processing and calling SNPs/Indels (mostly through biobc-nextgen), I'll use a group of other callers to get CNVs. Because the data I'll primarily be analyzing is WES, I'll focus on the tools that are specialized for that:

* XHMM
* EXCAVATORtool

Even though most of them use read depth to infer CNVs, they have different algorithms. As expected, each of their papers show how much better one tool is when compared to the other. So, the goal for CNVs from WES is to use SURVIVOR to combine the calls made by the different tools.

On the SNP chip front, I plan to stick with PennCNV, which has its own algorithm for denovo calling. I might try to apply it to the variants I find in the other arms as well. Scan2CNV looks like a wrapper to it, but I'm not sure if it implements all functionalities I need.

For SNPs/INDELS, I'll use multiple callers implemented in bcbio-gen, but I'm not sure yet whether it's a good idea to call them individually, or do some sort of joint/population call. If the latter, should we do it within families, or across all samples?

The final goal is to annotate the variants obtained in the 3 arms, and compare them to public (SV?) databases.
