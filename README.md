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
In order to visualize coverage along transcripts mapped reads are quantified at each genomic position and mapped onto most relevant gene annotated and its flanks. To be able to visualize data coverage on a standard metagene, all gene were binned where the size of sequence associated to the bin is proportionel to the gene length. Up and down stream flanks were binned in fixed sized bins (5 nt).

### Producing binnned metagene profiles
Binning of genomic positions happens in two steps, first coverage is calculated along each annotated gene selected and its flanks with the below script

    mk_feature_profile_ex.py

First appropriate data sets are produced

    mk_binned_metagene_ex_shortlist_50.sh
    
Then this script does the binning and calculations

    binnned_profile.py

### Plotting a metagene profile of a single dataset for different datesets selected by rate of transcription
To explore how ChIP-seq occupancy changes with rate of transcription (measured by pNET-seq) scripts have been developed to plot ChIP-seq occupancy on metagenes

### Plotting many metagene profiles in one plot
In one visualization genes are divided by transcrion rate in 10th, 20th,... quantiles, then line representing average genomic coverage is plotted for each group across the metagene.

    plot_histone_marks.R

### Plotting many metagene profiles in one plot
First an intermediate data format is created to easy processing of the big dataset. Then the data sets are loaded, z-score normalized and this script will plot data, disregarding genes with rate of transcription in the lower and higher 25th quantile.

    print_profile_ex_50pct.R

## 4. Analysis of pNET-seq and estimation of rate of transcription
Commandline for merging pNET-seq dataset can be found in 

    bash_commands.sh

Merging is done after initial overlap-analysis of each dataset in calc_nascent_expression.sh 
    
### Pipeline for mapping and processing of pNET-seq data.
To estimate rate of transcription from pNET-seq data it takes a number of steps of processing before mapping and quantification. This is all done in the following script. Only 2nd mate is used and mapped by STAR and quantification is done using custom scripts.

    map_NETseq.sh

### Processing.
To be able to map pNET-seq reads barcodes have to be trimmed from first 4 bases in both mates with.

    trimBases.py
    
these 8 base barcodes are used to remove duplicated reads with
    
    rem_dup.py

3'end adapters are also removed with cutadapt.

### Calculating the rate of transcription from pNET-seq data
Following scripts are needed for calculating rate of transcription. 

    calc_nascent_expression.sh
    
    GRO_seq_expression_readdict.py

### Estimating new gene boundries based on TSS-seq, TIF-seq and direct RNA-seq data
To improve the annotation and make it more tissue specific, boundaries of genes were corrected by data indicative of transcription start site and termination site. The data sets were normalized and aggregated and then genes were corrected by the maximal data point in a plausible region around annotated TSS and TTS. 

    estimate_representative_transcript.py

## More details described below

All important custom Python, Bash and R scripts used in the computational analyses of (<paper ref>) below are shared here.
## Retrieval of data
All histone mark ChIP-seq datasets are wild type Col-0 Arabidopsis Thaliana seedlings, 5 days to 3 weeks old (see supplementary table S1 for more information). These data sets were identified through queries to the SRA (1) and DNA Data Bank of Japan (DDBJ) (2). Upon identification, SRA-files were retrieved from the SRA FTP-server and uncompressed from SRA format using fastq-dump part of SRA tools.
## Mapping and calculation of genomic coverage
Before mapping, 3’adapters were removed by a custom script that removes the most frequent of 4 different commonly used adapter sequences from the 3’end of single-end reads or 3’end of both mates in paired-end reads with cutadapt (3). Reads from ChIP-seq libraries were aligned to the Arabidopsis Thaliana genome TAIR10 using the STAR Ver 2.60c aligner (4) (options: --outSAMmultNmax 1, --seedSearchStartLmax 30, --alignEndsType EndToEnd, --alignIntronMax 1). Upon mapping, aligned reads were sorted by samtools (5) and clustering, normalization and peak calling was performed with MACS2 (options: -w -S -g 1.35+08 -m 3,50) to produce bedgraph files of genomic coverage. With these options MACS2 uses reads on both strands to estimate the average fragment size and extends reads accordingly to correct shifted signal of single-end protocols.   
## Binned metagene profiles
To calculate mean coverage across a metagene, genomic coverage was overlapped with the above annotation of protein coding genes using bedtools intersect (6). The genomic coverage profiles along the gene body were then divided into 200 equal sized parts for each gene, then average coverage was calculated across all genes for each bin separately by a custom script. Upstream and downstream 500 nt flanks of the gene were divided into 100 bins each 5 nt and average genomic coverage was calculated like for the gene body. To quantify relative enrichment along the metagene and flank, z-scores were finally calculated by withdrawing the mean enrichment and dividing by the standard deviation and these Z-scores were plotted in R (7).
## Estimation of gene borders
To make boundaries of protein coding genes more tissue specific to Col-0 seedling the Araport11 gene annotation (8) was retrieved and in-house and published datasets of TSS-seq, TIF-seq and Direct RNA-Seq (DR-Seq) were merged and normalized into a four data tracks along the genome, one for transcription start site track for each strand (TSS) and one termination site track (TTS) for each strand. Normalization was done by first calculating the average coverage e for each dataset as the fraction of basecalls R over the genome length l.
    
    e=R/l
    
The basecall count in each position c_i was normalized to n_i as follows. 
    
    n_i=c_i/e
    
Assignment of a new gene TSS or TTS was done if there was a non-zero datapoint in the corresponding data track which was inside the boundary window around respectively TSS or TTS. If there were more than one data non-zero point the coordinate of the maximal datapoint was selected and assigned to the gene feature in question. These boundary windows were estimated by including all transcript boundaries and protein boundaries annotated for a gene, choosing the minimal and maximal boundaries as window borders around TSS and TTS separately. The idea is that we can only refine boundaries of UTRs and not shorten the protein or define a new UTR longer than ever observed.

## Quantification of rate of transcription
In order to quantify rate of transcription, pNET-seq datasets from (9) in all conditions (total RNA, Unphosphorylated, Serine 2P and Serine 5P) corresponding to run IDs in supplementary table S1 were retrieved from the SRA FTP server. These libraries are paired-end with UMIs in the first 4 positions in each mate. 3’end adapters were removed with cutadapt (3), the UMIs were trimmed and used for removal of reads with the same sequence, allowing 1 mismatch in the 8 base UMI.  The 2nd mate was mapped with STAR ver 2.6.0c (options: --outSAMmultNmax 1, --seedSearchStartLmax 30, --alignEndsType EndToEnd, --alignIntronMax 1) and all uniquely mapping reads were selected, strand inverted and overlapped with the above described gene annotation. Then overlapping reads from all libraries were merged, associated to genes by overlapping with gene annotation with bedtools intersect (6) and Transcripts Per Million (TPM) was calculated for each protein coding gene with a custom python script. Finally, each gene was associated the corresponding TPM value and the genes with highest (25th higher quantile) and lowest (25th lower quantile) rates of transcription exhibited length bias. So, these data in the 25th higher and lower quantiles were discarded from the final dataset, which was used to plot the metagene profiles.

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


