Processing CDR3 clonality from RNA-seq data.
============================================

July 22, 2016
Scott Brown

Citation: [Brown et al. Bioinformatics 2016](https://academic.oup.com/bioinformatics/article/33/8/1111/2730232/Defining-the-clonality-of-peripheral-T-cell)
doi: 10.1093/bioinformatics/btw810
https://academic.oup.com/bioinformatics/article/33/8/1111/2730232/Defining-the-clonality-of-peripheral-T-cell


--------------------------------------------------


This script assumes you have already run MiTCR on your RNA-seq data. If not, please 
see [Extracting TCR CDR3 sequences from RNAseq data](https://github.com/scottdbrown/TCR-from-RNAseq2015) for details.

##Requirements:##
This script was developed using R v3.1.1 on 64-bit Linux.
Required R packages: ggplot2, entropy.

##Description:##
*Note:* This script aims to keep things generalized, but some nuanced lines may still exist, and editing may be required to make it fit your data.

MiTCR output files need to be parsed and combined into a single tsv file with the following headers:

|library|sample|chain|aaSeq|nucSeq|abundance|
|---|---|---|---|---|---|

During analysis, this dataframe will be merged with a "sample" dataframe. This is mainly to ensure that samples which have no TCR yield are maintained through the analysis, but other sample variables can also be introduced at this point.

The script will remove non-productive CDR3 sequences, calculate the relative abundance of each CDR3, resolves ambiguity in CDR3 chain assignment, and classifies CDR3 clonotypes as being dominant or background (using control samples to determine the background level). It then calculates Shannon Entropy and estimates Tumor Purity for each sample.

The script has an example plot showing the clonal landscape of each sample, highlighting dominant CDR3s in each.

Finally, it includes code for analyzing TCR C gene expression. For our analysis, RSEM v1.2.29 was used to calculate gene expression levels for all transcripts, and transcripts from TCR C genes were extracted.

