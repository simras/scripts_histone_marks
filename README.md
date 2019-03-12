# Scripts used to produce metagene profiles of ChIP-seq histone marks in <i>Arabidopsis Thaliana<i>.

## 1. Retrieval of data
Will download datasets from SRA specified in a file SRRID.txt and unpack the SRA files to fastq.gz.

    download_data.sh

## 2. processing and mapping of ChIP-seq data

### Pipeline for mapping ChIP-seq histone Mark data (and RNA-seq with a few alterations)

    map_RNAseq_data.sh

### Script to estimate the most frequent 3'end adapter sequence in a dataset

    pick_adapter.sh

## 3. Making metagene profiles

### Scripts to produce binnned metagene profiles

    mk_binned_metagene_ex_shortlist_50.sh
    
It depends on this scripts, first profiles are calculated across each genomic region and plotted by

    mk_feature_profile_ex.py
    
Then the following script calculates averages across bins

    binnned_profile.py

### Script to plot a metagene profile of a single dataset for different gneesets selected by rate of transcript

    print_profile_ex_50pct.R

### Plotting many metagene profiles in one plot

    plot_histone_marks.R

## 4. Analysis of pNET-seq and estimation of rate of transcription
Commandline for merging pNET-seq dataset can be found in 

    bashcommands.sh

Merging is done after initial overlap-analysis of each dataset in calc_nascent_expression.sh 
    
### Pipeline for mapping and processing of pNET-seq data.
The pipeline for analysing pNET.seq depends on the below scripts for processing the pNET-seq data.

    map_NETseq.sh

### Processing.
To be able to map pNET-seq reads barcodes have to be trimmed from first 4 bases in both mates with.

    trimBases.py
    
these 8 base barcodes are used to remove duplicated reads with
    
    rem_dup.py

### Script that calculates the rate of transcription from pNET-seq data

    calc_nascent_expression.sh

### Script that can estimate new gene boundries based on TSS-seq, TIF-seq and PAS-seq data in bedgraph format

    estimate_representative_transcript.py

## Details described below

All custom Python, Bash and R scripts used in the computational analyses below are shared here.

## Retrieval of data
All histone mark ChIP-seq datasets are wild type Col-0 Arabidopsis Thaliana seedlings, 5 days to 3 weeks old (see supplementary table S1 for more information). These data sets were identified through queries to the SRA (1) and DNA Data Bank of Japan (DDBJ) (2). Upon identification, SRA-files were retrieved from the SRA FTP-server and uncompressed from SRA format using fastq-dump.

## Mapping and calculation of genomic coverage
Before mapping, 3’adapters were removed by a custom script that removes the most frequent of 4 different commonly used adapter sequences from the 3’end of single-end reads or 3’end of both mates in paired-end reads with cutadapt (3). Reads from ChIP-seq libraries were aligned to the Arabidopsis Thaliana genome TAIR10 using the STAR Ver 2.60c aligner (4) (options: --outSAMmultNmax 1, --seedSearchStartLmax 30, --alignEndsType EndToEnd, --alignIntronMax 1). Upon mapping, aligned reads were sorted by samtools (5) and clustering, normalization and peak calling was performed with MACS (options: -w -S -g 1.35+08 -m 3,50) to produce bedgraph files of genomic coverage.

## Binned metagene profiles
To calculate mean coverage across a metagene, genomic coverage was overlapped with protein coding genes from the above custom gene annotation with bedtools intersect (6). The genomic coverage profiles along the gene body were then divided into 100 equal sized parts, averages were calculated across values given by MACS in each bin and then across all genes by a custom script and this was plotted in R (7).

## Estimation of gene borders
To make boundaries of protein coding genes more tissue specific to Col-0 seedling the Araport11 gene annotation (8) was retrieved andin-house and published datasets of TSS-seq, TIF-seq and PAS-seq were merged and normalized into a single data track for each strand separately for transcription start sites (TSS) and termination sites (TTS) along the genome. Normalization was done by first calculating the average coverage e for each dataset as the fraction of basecalls R over the genome length l.

e=R/l

The basecall count in each position was normalized to e for coverage to be comparable between datasets. After this, the coverage in position i, ci was calculated as ni and it was added across all datasets to produce a single track separated by strands and whether the data is indicative of TSS or TTS.

n_i=c_i/e

## Quantification of rate of transcription
In order to quantify rate of transcription pNET-seq datasets from (9) in all conditions (total RNA, Unphosphorylated, Serine 2P and Serine 5P) corresponding to run IDs in supplementary table S1 were retrieved from the SRA FTP server. These libraries are paired-end with UMIs in the first 4 positions in each mate, where the first mate in shorter than the 2nd mate. Thus, the UMIs were trimmed and used for removal of reads that map to the same position and strand, allowing 1 mismatch in the 8 base UMI. The 2nd mate was mapped with STAR Ver 2.6.0c (options: --outSAMmultNmax 1, --seedSearchStartLmax 30, --alignEndsType EndToEnd, --alignIntronMax 1) and all uniquely mapping reads were selected, strand inverted and overlapped with the above gene annotation. Then overlapping reads from all libraries were merged and Transcripts Per Million (TPM) was calculated for each protein coding gene. 

## References
1. 	Leinonen R, Sugawara H, Shumway M. The sequence read archive. Nucleic Acids Res. 2011; 
2. 	Kaminuma E, Mashima J, Kodama Y, Gojobori T, Ogasawara O, Okubo K, et al. DDBJ launches a new archive database with analytical tools for next-generation sequence data. Nucleic Acids Res. 2009; 
3. 	Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal. 2011; 
4. 	Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, et al. STAR: Ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15–21. 
5. 	Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009; 
6. 	Quinlan AR. BEDTools: The Swiss-Army tool for genome feature analysis. Curr Protoc Bioinforma. 2014;2014:11.12.1-11.12.34. 
7. 	R Development Core Team R. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing. 2011. 
8. 	Cheng CY, Krishnakumar V, Chan AP, Thibaud-Nissen F, Schobel S, Town CD. Araport11: a complete reannotation of the Arabidopsis thaliana reference genome. Plant J. 2017; 
9. 	Zhu J, Liu M, Liu X, Dong Z. RNA polymerase II activity revealed by GRO-seq and pNET-seq in Arabidopsis. Nat Plants [Internet]. 2018;4(12):1112–23. Available from: http://dx.doi.org/10.1038/s41477-018-0280-0

