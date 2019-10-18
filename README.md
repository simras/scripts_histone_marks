# Scripts used to produce metagene profiles of ChIP-seq histone marks in <i>Arabidopsis Thaliana<i>.
The following set of scripts define a pipeline for retrieval of sequencing libraries in Sequence Read Archive (SRA), an analysis pipeline and scripts for data preparation and visualization. To run they require installation of a number of freely available Bioinformatics tools and paths to these may need to be changed in the scripts for proper operation.
    
## 1. Retrieval of data
To retrieve and uncompress the sequencing data files that will be analyzed, their run ID which starts with SRR (if submitted in the US server) has to be specified in a file (SRRID.txt) and this file should be present at a path specified in the below script. The script will produce a fastq.gz data format as output and organize the data in the specified folder structure. The name of the script is as follows.

    download_data.sh

## 2. processing and mapping of ChIP-seq data

### Pipeline for mapping ChIP-seq histone Mark data (and RNA-seq with a few alterations)
This pipeline has been used to analyze paired-end and single end ChIP-seq data and RNA-seq data. The STAR aligner is used to map against the genome for the ChIP-seq data.

    map_RNAseq_data.sh

### Estimating which is the most frequent 3'end adapter sequence in a dataset
In order to map the ChIP-seq data 3'end adapter sequences had to be removed in most datasets. We do not know with certainty which is the sequence of this adapter, but we know the most commonly used adapter sequences. The following script tests which of the know adapters is present in most reads and selects this for later trimming.

    pick_adapter.sh

## 3. Making metagene profiles
In order to visualize coverage along transcripts mapped reads are quantified at each genomic position and mapped onto most relevant gene annotated and its flanks. To be able to visualize data coverage on a standard metagene, all gene were binned where the size of sequence associated to the bin is proportionel to the gene length. Up and down stream flanks were binned in fixed sized bins (10 nt).

### Producing binnned metagene profiles
Binning of genomic positions happens in two steps, first coverage is calculated along each annotated gene selected and its flanks with the below script

    mk_feature_profile_ex.py

Then these profiles are binned

    mk_binned_metagene_ex_shortlist_50.sh  
    
Then the following script calculates averages across bins

    binnned_profile.py

### Script to plot a metagene profile of a single dataset for different datesets selected by rate of transcription
To explore how ChIP-se occupancy changes with rate of transcription (measured by pNET-seq).

    print_profile_ex_50pct.R

### Plotting many metagene profiles in one plot

    plot_histone_marks.R

## 4. Analysis of pNET-seq and estimation of rate of transcription
Commandline for merging pNET-seq dataset can be found in 

    bash_commands.sh

Merging is done after initial overlap-analysis of each dataset in calc_nascent_expression.sh 
    
### Pipeline for mapping and processing of pNET-seq data.
The pipeline for analysing pNET.seq is the following and it depends on the below scripts for processing the pNET-seq data.

    map_NETseq.sh

### Processing.
To be able to map pNET-seq reads barcodes have to be trimmed from first 4 bases in both mates with.

    trimBases.py
    
these 8 base barcodes are used to remove duplicated reads with
    
    rem_dup.py

### Script that calculates the rate of transcription from pNET-seq data

    calc_nascent_expression.sh
    
### Script to calculate TPM
    GRO_seq_expression.py

### Script that can estimate new gene boundries based on TSS-seq, TIF-seq and PAS-seq data in bedgraph format

    estimate_representative_transcript.py

## More details described below

All important custom Python, Bash and R scripts used in the computational analyses below are shared here.

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

