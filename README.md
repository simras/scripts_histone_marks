Scripts used to produce metagene profiles of ChIP-seq histone marks in Arabidopsis Thaliana.

Retrieval of data
Will download datasets from SRA specified in a file SRRID.txt and unpack the SRA files to fastq.gz.

download_data.sh


calc_nascent_expression.sh
estimate_representative_transcript.py

Pipeline for mapping and processing of pNET-seq data.
map_NETseq.sh

Pipeline for mapping ChIP-seq histone Mark data (and RNA-seq with a few alterations)
map_RNAseq_data.sh

Script to estimate the most frequent 3'end adapter sequence in a dataset
pick_adapter.sh

Script to produce binnned metagene profiles
mk_binned_metagene_ex_shortlist_50.sh

Script to plot a metagene profile of a single dataset for different gneesets selected by rate of transcript
print_profile_ex_50pct.R

Script to plot many metagene profiles in one plot
plot_histone_marks.R
