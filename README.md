
# LoopExtrusion

Loop extrusion as a mechanism for DNA Double-Strand Breaks repair foci
formation

## Overview

This is the github repository of our paper [Loop extrusion as a
mechanism for DNA Double-Strand Breaks repair foci
formation](https://www.biorxiv.org/content/10.1101/2020.02.12.945311v1).

## Data Availability

High throughput sequencing data have been deposited to Array Express
under accession number
[E-MTAB-8851](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8851/).
ChIP-chip data have been deposited to Array Express under accession
number
[E-MTAB-8793](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8793/).
Other data and source codes are available upon request.

## System requirements

### Alignment

The workflow was written with snakemake (python), and need some genomic
tools :

  - `bwa 0.7.12-r1039`.
  - `samtools 1.9`.
  - `FastQC v0.11.5`.
  - `macs2 2.1.2`.
  - `bedtools v2.26.0`.
  - `deeptools 3.4.3`.
  - `R version 3.6.3`.

### Figure generation

the scripts were written with `R`, and need some packages
:

  - `library(rtracklayer)`.
  - `library(BSgenome.Hsapiens.UCSC.hg19)`.
  - `library(ggplot2)`.
  - `library(reshape2)`.
  - `library(dplyr)`.
  - `library(Homo.sapiens)`.

| Scripts                                              | Description                                                                                      | Figures                             |
| :--------------------------------------------------- | :----------------------------------------------------------------------------------------------- | :---------------------------------- |
| analysis\_Virtual4C\_DSB.R                           |                                                                                                  |                                     |
| insulation\_score\_profil.R                          |                                                                                                  |                                     |
| make\_APA\_heatmap.R                                 | Get the APA from juicer tools and produce graphics with ggplot2.                                 | Figs 2c; ext2f; ext2g; ext6f        |
| make\_4Cseq\_barplot.R                               | Compute differential 4C-seq signal (log2 +DSB/-DSB) on 1mb around DSBs viewpoints.               | Figs 4c; ext3d; ext3h               |
| make\_chipseq\_barplot\_TADs.R                       | Compute ChIP-seq quantification within the damaged TAD and neighboring TADs.                     | Figs ext1f                          |
| make\_chipseq\_boxplot\_cohesine.R                   | Compute signal over cohesine peaks.                                                              | Figs ext6b                          |
| make\_chipseq\_boxplot\_profile.R                    | Compute signal over DSBs and produce some profiles/boxplots.                                     | Figs ext1d; ext1k; 2e; ext2b; ext7b |
| make\_HiC\_heatmap.R                                 | Compute average heatmaps from HiC dumped matrixes (juicer) and HiTC.                             | Figs 2b; 2d; 2g; ext2d; ext2e       |
| make\_HiTC.R                                         | Produce HiTC format from dumped matrix to be loaded in R.                                        |                                     |
| make\_loopanchor\_chipseq\_boxplot\_over\_distance.R | Compute the quantification of SCC1 recruitment on loop anchors at different distances from DSBs. | Figs ext6d                          |
| make\_loopanchor\_chipseq\_boxplot.R                 | Compute the ChIP-seq signal on loop anchors.                                                     | Figs ext7c                          |
| make\_loopanchor\_HiC\_boxplot.R                     | Compute the differential loop strength in undamaged or damaged TADs                              | Figs ext6g                          |
| plot\_gH2AX\_TADborder.R                             | Compute average gammaH2AX profile centered to the closest TAD border to the DSBs.                | Figs 1e                             |
