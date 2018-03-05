# small-RNA Expression Analysis Pipeline
Expression analysis of cellular and plasma micro RNA
The following pipeline allows you to trim small RNA reads generated from Qia-Seq protocol.
#Reads are then aligned using very sensitive settings in Bowtie2 and counted using GenomicRanges
#Differential gene expression analysis using edgeR

## Step1: Run small_RNA_trimming.sh [Input file]
## Step2: Generate the targets file, data, and results directory as described in systemPipeR
## Step3: Run small_RNA_alignment.sh line by line
## Step4: Run small_RNA_readcounting.sh line by line
## Step5: Run small_RNA_edgeR.R line by line
