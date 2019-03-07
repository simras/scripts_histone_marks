# Scripts used to produce metagene profiles of ChIP-seq histone marks in Arabidopsis Thaliana.

## 1. Retrieval of data
Will download datasets from SRA specified in a file SRRID.txt and unpack the SRA files to fastq.gz.

    download_data.sh

## 2. processing and mapping of ChIP-seq data

### Pipeline for mapping ChIP-seq histone Mark data (and RNA-seq with a few alterations)

    map_RNAseq_data.sh

### Script to estimate the most frequent 3'end adapter sequence in a dataset

    pick_adapter.sh

## 3. Making metagene profiles

### Script to produce binnned metagene profiles

    mk_binned_metagene_ex_shortlist_50.sh

### Script to plot a metagene profile of a single dataset for different gneesets selected by rate of transcript

    print_profile_ex_50pct.R

### Script to plot many metagene profiles in one plot

    plot_histone_marks.R

## 4. Analysis of pNET-seq and estimation of rate of transcription

### Pipeline for mapping and processing of pNET-seq data.

    map_NETseq.sh

### Script that calculates the rate of transcription from pNET-seq data

    calc_nascent_expression.sh

### Script that can estimate new gene boundries based on TSS-seq, TIF-seq and PAS-seq data in bedgraph format

    estimate_representative_transcript.py

Details described below
