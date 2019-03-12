#!/bin/bash
#
# Automatic ChIP-seq mapping pipeline
#

#set -e
set -x

#TAIR10
baseGDIR="/home/simras/Work/histone_mark_data/genome"
genome="TAIR10.fa"

#Download_Folder
downloads=data/
STAR="/home/simras/Work/STAR/bin/Linux_x86_64/STAR"
MACS=/home/simras/Work/MACS-1.4.2/bin/macs14
bam=Aligned.out.bam

adapt="AGATCGGAAGAGC"

i=0
j=0
mapped=0

for SRR in SRR5833577 SRR5833580 SRR5833583
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
    
    new_bam=/home/simras/data/mapped/$SRR_ID\_$genome.bam
    sorted_bam=/home/simras/data/mapped/$SRR_ID\_$genome"_"sorted.bam
    if [ ! -f $sorted_bam ]
    then
	mapped=1
	if [ $PAIRED == "T" ]
	then
	    ad1=$(./pick_adapter.sh data/$SRR_ID\_1.fastq.gz)
	    ad2=$(./pick_adapter.sh data/$SRR_ID\_2.fastq.gz)
	    datno1=$downloads/$SRR_ID"_noadt_1.fastq.gz"
	    datno2=$downloads/$SRR_ID"_noadt_2.fastq.gz"
	    #cutadapt -a $ad1 -A $ad2 -o $datno1 -p $datno2 $dat1 $dat2
	    if [ ! -f $datno2 ]
	    then
		cutadapt -a $ad2 -j 4 -o $datno2 $dat2
	    fi
	    if [ ! -f $datno1 ]
	    then
		cutadapt -a $ad1 -j 4 -o $datno1 $dat1
	    fi
	    $STAR --genomeDir $baseGDIR  --readFilesIn $datno1 $datno2 --runThreadN 8 --outSAMtype BAM Unsorted --readFilesCommand "zcat" --seedSearchStartLmax 30 --alignEndsType EndToEnd --alignIntronMax 1 --outSAMmultNmax 1
	   
	else
	    $STAR --genomeDir $baseGDIR  --readFilesIn $dat1 --runThreadN 8 --outSAMtype BAM Unsorted --clip3pAdapterSeq $adapt --readFilesCommand "zcat" --seedSearchStartLmax 30 --alignEndsType EndToEnd --alignIntronMax 1 --outSAMmultNmax 1
	fi
	
	if [ $? -eq 0 ]
	then	    
	    new_bam=mapped/$SRR_ID\_$genome.bam
	    mv $bam $new_bam
	    mv Log.final.out logs/$SRR_ID\_Log.final_2.out
	    mv Log.out logs/$SRR_ID\_Log_2.out
	    mv Log.progress.out logs/$SRR_ID\_Log.progress_2.out
	    
	    samtools sort $new_bam > $sorted_bam
	    rm $new_bam
	    if [ $? -ne 0 ]
	    then
		echo "samtools sort error" $SRR_ID > errors/$SRR_ID\_error.report
	    fi
	else
	    echo "Mapping error " $SRR_ID > errors/$SRR_ID\_error.report    
	fi
    fi
    rm  bedgraphs/$SRR_ID"_treat_afterfiting_all.bdg.gz"
    if [ -f $sorted_bam ] && [ ! -f bedgraphs/$SRR_ID"_treat_afterfiting_all.bdg.gz" ]
    then
	#cd mapped
	if [ $PAIRED == "T" ]
	then
	    touch bedgraphs/$SRR_ID"_treat_afterfiting_all.bdg.gz"
	    sorted_bam_PE=/home/simras/data/mapped/$SRR_ID\_$genome"_"sorted_PE.bam
	    # histone mark data sets
	    #samtools view -b -f 0x2 $sorted_bam > $sorted_bam_PE
	    #$MACS -t $sorted_bam -f BAM -w -n $SRR_ID -S -g 1.35e+08 -m 3,50 --bdg --nomodel
	    bedtools genomecov -trackline -bg -ibam $sorted_bam -strand +|gzip -9 > bedgraphs/$SRR_ID"_+_treat_afterfiting_all.bdg.gz"
	    bedtools genomecov -trackline -bg -ibam $sorted_bam -strand -|gzip -9 > bedgraphs/$SRR_ID"_-_treat_afterfiting_all.bdg.gz"
	    
	    bedtools genomecov -trackline -bg -ibam $sorted_bam -strand +|gzip -9 > bedgraphs/$SRR_ID"_-+_treat_afterfiting_all.bdg.gz"
	    bedtools genomecov -trackline -bg -ibam $sorted_bam -strand -|gzip -9 >> bedgraphs/$SRR_ID"_-+_treat_afterfiting_all.bdg.gz"
	else
	    touch bedgraphs/$SRR_ID"_treat_afterfiting_all.bdg.gz"
	    #$MACS -t $sorted_bam -f BAM -w -n $SRR_ID -S -g 1.35e+08 -m 3,50 --bdg
	    bedtools genomecov -trackline -bg -ibam $sorted_bam -strand +|gzip -9 > bedgraphs/$SRR_ID"_+_treat_afterfiting_all.bdg.gz"
	    bedtools genomecov -trackline -bg -ibam $sorted_bam -strand -|gzip -9 > bedgraphs/$SRR_ID"_-_treat_afterfiting_all.bdg.gz"

	fi
	# cleanup
	mv $SRR_ID\_MACS_bedGraph/treat/$SRR_ID"_treat_afterfiting_all.bdg.gz" bedgraphs/
	rm -r $SRR_ID"_peaks"* $SRR_ID"_summits.bed" $SRR_ID"_model.r" $SRR_ID"_MACS_bedGraph"
    fi
done
