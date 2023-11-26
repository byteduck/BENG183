# Gene Expression Analysis with RNA-seq

## Introduction

Understanding gene expression is crucial for unraveling the molecular mechanisms underlying various biological processes, such as development, disease, and response to environmental stimuli. One powerful tool for studying gene expression is RNA sequencing (RNA-seq), which allows for a comprehensive quantitative assessment of the transcriptome, or in other words, the complete set of RNA molecules in a cell or tissue.

## RNA-seq Overview

### Principle

RNA-seq is a high-throughput sequencing technique that involves the following key steps:

1. **RNA Extraction**: Total RNA is isolated from a sample of interest, which can be cells, tissues, or even single cells.

2. **Library Preparation**: The extracted RNA is converted into a library of cDNA (complementary DNA) molecules. This step typically involves the fragmentation of RNA and the addition of adapter sequences.

3. **Sequencing**: The cDNA library is sequenced using next-generation sequencing (NGS) technology, generating millions of short reads.

5. **Data Analysis**: The sequenced reads are processed to infer gene expression levels and identify differentially expressed genes.

### RNA-seq Pipeline

Here is a basic overview of the tools used in the pipeline for RNA-Seq differential expression analysis:

| Tool          | Input                                         | Output                                | Purpose                                          |
|---------------|-----------------------------------------------|---------------------------------------|--------------------------------------------------|
| `fastqc`        | Fastq Reads                                   | Quality Report HTML                   | Check read quality                               |
| `Trim Galore`   | Fastq Reads                                   | Trimmed fastq reads                   | Trim adapter sequences                           |
| `STAR`          | Trimmed fastq reads + Transcriptome reference | SAM / BAM aligned reads               | Align reads to reference genome or transcriptome |
| `featureCounts` | SAM/BAM + Reference                           | Count matrix                          | Quantify expression                              |
| `DESeq2`        | Count Matrices for Samples                    | Differential expression for each gene | Figure out difference between samples            |

Figure 1: Overview of the RNA-seq analysis pipeline. The process includes RNA extraction, library preparation, sequencing, and data analysis.

### File Formats

- **FASTQ**: Raw sequencing data, containing sequence reads and their associated quality scores.

- **SAM/BAM**: Sequence Alignment/Map files for mapping reads to a reference genome.

- **GTF/GFF**: Annotation files defining gene and transcript structures.

- **Counts Table**: A table that quantifies the number of reads or fragments aligned to each gene.

### Altogether...
Essentially, the bioinformatic RNA-Seq pipeline looks something like this:
- First, we use fastqc to determine the quality of our reads.
- If everything checks out, we use a tool, such as trim galore, to trim the adapter sequences (if using a technology like Illumina sequencing, for example), and check the quality again.
- Next, we use an alignment tool such as STAR along with our trimmed reads and reference transcriptome to align the reads to the reference transcriptome. This will output aligned reads.
- Next, we use a tool like featureCounts to quantify the expression of our mapped reads, to get a count matrix.
- Finally, we use a tool such as `DESeq2` to compare counts between samples, and figure out the fold change.

## Quantification of Gene Expression

To quantify gene expression, RNA-seq data is typically analyzed in terms of read counts or normalized expression values. Three commonly used metrics for expression quantification are RPKM, FPKM, and TPM.

### RPKM (Reads Per Kilobase Million)

RPKM measures gene expression as the number of reads mapping to a gene per kilobase of its coding sequence per million mapped reads. It normalizes for gene length and sequencing depth.

![RPKM Formula](rpkm.png)

Figure 2: RPKM calculation formula.

### FPKM (Fragments Per Kilobase Million)

FPKM is similar to RPKM but considers fragments (paired-end reads) rather than individual reads. It is especially useful when working with paired-end sequencing data.

![FPKM Formula](fpkm.png)

Figure 3: FPKM calculation formula.

### TPM (Transcripts Per Million)

TPM quantifies gene expression as the proportion of transcripts attributed to a specific gene, normalized to the total number of transcripts in the sample.

![TPM Formula](tpm.png)

Figure 4: TPM calculation formula.

TPM is typically the preferred measure of gene expression for differential expression analysis for a handful of good reasons:

1. **Normalization for Transcriptome Size**: TPM is normalized to the total RNA reads in each sample. This means that the sum of all TPM values in a sample is the same, which makes it a lot easier to compare expression levels across different samples. In contrast, RPKM and FPKM instead normalize for sequencing depth and gene length, but not for the total number of reads, which can vary a lot between samples - making comparing them between samples difficult.
2. **Better Handling of RNA Composition Effects**: TPM is more effective in accounting for differences in RNA make-up between samples. In RPKM and FPKM, a disproportionately expressed gene can take up a huge portion of the reads sampled, which can skew the perceived expression levels of other genes. TPM, on the other hand, fixes this issue by normalizing each gene’s expression to the total transcriptome size, making it less susceptible to variations caused by a handful of highly expressed genes.
3. **Look-ability**: TPM provides an instantly intuitive understanding of the data. A TPM value of 100, for instance, means that for every million transcripts in your sample, 100 are from the gene of interest. This direct interpretation is not as straightforward with RPKM or FPKM; and this way, it is a lot more obvious just by looking at the numbers what sort of difference there may be between samples. Not that we don't still need to run statistics, though!

## RNA-Seq Differential Expression Analysis

### So, why differential expression?

In certain situations, it can be useful to compare quantified gene expression using RNA-Seq between samples; for instance:
- **Identifying Disease Mechanisms**: By analyzing gene expression with respect to diseases and infections, we can figure out what pathways they use and possible treatments.
- **Drug Development and Testing**: By figuring out what genes drugs affect the expression of, we can better understand them and develop new ones.
- **Cancer Gene Identification**: By identifying genes that are expressed differently in cancer cells, we can create better treatments.
- **Gene Regulation Research**: Knowing which genes regulate others can be helpful in developing treatments and technologies that make use of the natural mechanisms of the genome.

### DESeq2 for differential expression

There are many different tools that can be used to perform differential expression analysis, but the one we'll be explaining is [**DESeq2**](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

### DESeq2 pipeline
<Steps for using DESEq2>

### DESeq2 example output
<example output and link to dataset, maybe from a paper?>

### The magic behind DESeq2
<Explanation of statistics behind DESeq2>

### DESeq2 vs Cuffdiff

