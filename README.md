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

First appropriate data sets are produced

    mk_binned_metagene_ex_shortlist_50.sh
    
Then this script does the binning and calculations

    binnned_profile.py

### Plotting a metagene profile of a single dataset for different datesets selected by rate of transcription
To explore how ChIP-seq occupancy changes with rate of transcription (measured by pNET-seq) scripts have been developed to plot ChIP-seq occupancy on metagenes

### Plotting many metagene profiles in one plot
In one visualization genes are divided by transcrion rate in 10th, 20th,... quantiles, the an occupancy line is plotted for each group across the metagene.

    plot_histone_marks.R

### Plotting many metagene profiles in one plot
First an intermediate data format is created to easy processing of the big dataset. Then the data sets are loaded z-score normalized and this script will plot data on genes with rate of transcription not in the lower and higher 25th quantiles.

    print_profile_ex_50pct.R

## 4. Analysis of pNET-seq and estimation of rate of transcription
Commandline for merging pNET-seq dataset can be found in 

    bash_commands.sh

Merging is done after initial overlap-analysis of each dataset in calc_nascent_expression.sh 
    
### Pipeline for mapping and processing of pNET-seq data.
To estimate rate of transcription from pNET-seq data it takes a number of step of processing. The following pipelione does

    map_NETseq.sh

### Processing.
To be able to map pNET-seq reads barcodes have to be trimmed from first 4 bases in both mates with.

    trimBases.py
    
these 8 base barcodes are used to remove duplicated reads with
    
    rem_dup.py

### Calculating the rate of transcription from pNET-seq data

    calc_nascent_expression.sh
    
### Script to calculate TPM

    GRO_seq_expression.py

### Estimating new gene boundries based on TSS-seq, TIF-seq and direct RNA-seq data
To improve the annotation and make it more tissue specific, boundaries of genes were corrected by data indicative of transcription start site and termination site. The data sets were normalized and aggregated and then genes were corrected by the maximal data point in a plausible region around annotated TSS and TTS. 

    estimate_representative_transcript.py

## More details described below

All important custom Python, Bash and R scripts used in the computational analyses of (<paper ref>) below are shared here.
    
