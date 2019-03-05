#!/bin/bash
#
# Automatic ChIP pipeline
#

#TAIR10
baseGDIR="/home/simras/Work/histone_mark_data/genome"
genome="TAIR10.fa"
#genome2="LFQT01.2.fsa_nt.fa"
#Download_Folder
downloads=$(pwd)/$1
dat1="data/SRR2452455.fastq"
#dat2="RNAseq/SRR2452456.fastq"
#dat3="RNAseq/SRR2452457.fastq"
STAR="../STAR/bin/Linux_x86_64/STAR"
bam=Aligned.sortedByCoord.out.bam
data=data
set -x
#set -e
#
#$STAR --runMode genomeGenerate --genomeFastaFiles $baseGDIR/$genome --genomeDir $baseGDIR/ --runThreadN 8


# process URL-List
i=0
j=0
k=0
space_consumption=0
while read ID
do

    if [ -z $(cat blacklist.txt |grep $ID) ]
    then
	#    echo $link
	# Download raw data wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR001/SRR001115/SRR001115.sra
	#echo ${ID:0:3}
	if [ ${ID:0:3} == "SRR" ]
	then
	    fname=$ID.sra
	    fname2=$ID.fastq.gz	
	    if [ -f $data/$ID\_1.fastq.gz ]
	    then
		PAIRED=1
		fname2=$ID\_1.fastq.gz
	    else
		fname2=$ID.fastq.gz
		PAIRED=0
	    fi
	    link=ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${ID:0:6}/$ID/$fname
	else
	    if [ ${ID:0:3} == "DRR" ]
	    then
		echo "JAPANESE DATABANK DOWNLOAD MANUALLY", $fname
	    fi
	    fname=$ID.fastq.gz
	    link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${ID:0:6}/$ID/$fname
	fi
	if [ -f $data/$fname2 ]
	then
	    echo "Already processed" $fname
	else
	    if [ ! -f $data/$fname ]
	    then
		cd $data/
		wget -t 10 --timeout=10 $link 1> /dev/null #2> /dev/null
		cd ..
	    else
		echo "Already downloaded" $fname	    
	    fi
	fi
	if [ ! -f $data/$fname2 ]
	then
	    echo "Processing file" $fname
	    cd $data/
	    fastq-dump -split-3 --gzip $fname
	    rm $fname
	    cd ..
	else
	    rm $data/$fname
	fi
    else
	echo $ID "is on blacklist"
    fi
done < SRRIDs.txt
