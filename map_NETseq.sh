#!/bin/bash
#
# Automatic ChIP pipeline
#

#set -e
set -x


# Bedgraph folder
bedgraphs="."
mapped="."
logs="."

# data folder
data="."

#TAIR10
baseGDIR="/home/simras/Work/histone_mark_data/genome"
genome="TAIR10.fa"

#Download_Folder
downloads="."

STAR="/home/simras/Work/STAR/bin/Linux_x86_64/STAR"
MACS=/home/simras/Work/MACS-1.4.2/bin/macs14
bam=Aligned.out.bam

adapt="AGATCGGAAGAGC"

i=0
j=0
mapped_b=0

# put run IDs
for SRR in $(cat SRRIDs.txt)
do
    tmp=${SRR%.fastq.gz}
    tmp1=${tmp%_*}
    SRR_ID=${tmp1##*/}
    
    # map data to assembly 1
    dat1=$downloads/$SRR_ID".fastq.gz"
    dat2=$downloads/$SRR_ID"_2.fastq.gz"
    
    if [ -f $dat2 ]
    then
	PAIRED="T"
	dat1=$downloads/$SRR_ID"_1.fastq.gz"
    else
	adapt=$(./pick_adapter.sh $dat1)   
	PAIRED="F"
	dat2=""
    fi
    
    new_bam=$mapped/$SRR_ID\_$genome.bam
    sorted_bam=$mapped/$SRR_ID\_$genome"_"sorted.bam
    #rm  $sorted_bam
    if [ ! -f $sorted_bam ]
    then
	mapped_b=1
	if [ $PAIRED == "T" ]
	then
	    datno1=$downloads/$SRR_ID"_nodups_1.fastq.gz"
	    datno2=$downloads/$SRR_ID"_nodups_2.fastq.gz"
	    
	    ad1=$(./pick_adapter.sh $dat1)
	    ad2=$(./pick_adapter.sh $dat2)
	    datno1=$downloads/$SRR_ID"_noadt_1.fastq.gz"
	    datno2=$downloads/$SRR_ID"_noadt_2.fastq.gz"
	    	    
	    dat_1_trim=$downloads/$SRR_ID"_trimmed_1.fastq.gz"
	    dat_2_trim=$downloads/$SRR_ID"_trimmed.fastq.gz" 
	    cutadapt -a $ad1 -j 4 -o $datno1 $dat1
	    cutadapt -a $ad2 -j 4 -o $datno2 $dat2
	    
	    ./trimBases.py -f $datno1 -5 4 -3 4 -o $SRR_ID"_barcodes_1.txt" 1> /dev/null 
	    ./trimBases.py -f $datno2 -5 4 -3 4 -o $SRR_ID"_barcodes_2.txt"|gzip -9 > $dat_2_trim

	    $STAR --genomeDir $baseGDIR  --readFilesIn ../$dat_2_trim --runThreadN 8 --outSAMtype BAM Unsorted --readFilesCommand "zcat" --seedSearchStartLmax 30 --alignEndsType EndToEnd --alignIntronMax 1 --outSAMmultNmax 1
	    mv Aligned.out.bam "Aligned_"$SRR_ID".bam"
	    
	    join $SRR_ID"_barcodes_1.txt" $SRR_ID"_barcodes_2.txt"|awk '{con=$2$3;print $1"\t"con}' > $SRR_ID"_barcodes.txt"
	    cat $SRR_ID"_barcodes.txt"|sort -k1,1 > $SRR_ID"_tmp.txt"
	    mv $SRR_ID"_tmp.txt" $SRR_ID"_barcodes.txt" 
	else
	    $STAR --genomeDir $baseGDIR  --readFilesIn $dat1 --runThreadN 8 --outSAMtype BAM Unsorted --clip3pAdapterSeq $adapt --readFilesCommand "zcat" --seedSearchStartLmax 30 --alignEndsType EndToEnd --alignIntronMax 1 --outSAMmultNmax 1
	fi
	new_bam=$mapped/$SRR_ID\_$genome.bam
	
	if [ ! -f $new_bam ]
	then
	    new_sam=$mapped/$SRR_ID\_$genome.sam
	    samtools view -H "Aligned_"$SRR_ID".bam" > $new_sam
	    samtools view "Aligned_"$SRR_ID".bam"|sort -k1,1|./rem_dup.py -d $SRR_ID"_barcodes.txt" >> $new_sam 
	    samtools view -Sb $new_sam > $new_bam
	    rm $new_sam
	    mv Log.final.out $logs/$SRR_ID\_Log.final_2.out
	    mv Log.out $logs/$SRR_ID\_Log_2.out
	    mv Log.progress.out $logs/$SRR_ID\_Log.progress_2.out
	    samtools sort $new_bam > $sorted_bam
	    rm $new_bam
	fi
    fi
done
