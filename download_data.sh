#!/bin/bash
#
# Automatic ChIP pipeline
#

# TAIR10 genome location
baseGDIR="/home/simras/Work/histone_mark_data/genome"
genome="TAIR10.fa"

# Download_Folder
downloads=$(pwd)/$1

# Data folder
data=$downloads

set -x

# SRA collects data from the Japanese Databank, EBI and NCBI. Origin of data is coded in the first three letters of the run ID
# SRR, NCBI
# ERR, EBI
# DRR, Japanese Databank
# These datacollections hav different download links which ccan be used for automated download.

# process URL-List
i=0
j=0
k=0
space_consumption=0
while read ID
do

    if [ -z $(cat blacklist.txt |grep $ID) ]
    then
	#
	# Download raw data wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR001/SRR001115/SRR001115.sra
	
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
		echo "JAPANESE DATABANK", $fname
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
