#!/bin/bash
#
# Automatic adapter detection in ChIP pipeline
#
# Use:
# for f in data/*.fastq.gz;do ./pick_adapter.sh $f;done

#set -x

# Standard illumina adapter
ad1=AGATCGGAAGAGC

# 
ad2=GCGGGTTGAGCCAAGCCTGTGCAGTA

#
ad4=GATCGTCGGACTGTAGAAC

# 
ad3=TGGAATTCTCGGGTGCCAA

ll=$(zcat $1|head -2|tail -1|wc -c)

c1=$(zcat $1|head -400000 |awk 'BEGIN{i = 0}{if(i % 4 == 1){print $0};i = i + 1}'|grep -c $ad1)
c2=$(zcat $1|head -400000 |awk 'BEGIN{i = 0}{if(i % 4 == 1){print $0};i = i + 1}'|grep -c $ad2)
c3=$(zcat $1|head -400000 |awk 'BEGIN{i = 0}{if(i % 4 == 1){print $0};i = i + 1}'|grep -c $ad3)
c4=$(zcat $1|head -400000 |awk 'BEGIN{i = 0}{if(i % 4 == 1){print $0};i = i + 1}'|grep -c $ad4)
call=$((c1 + c2 + c3 + c4))

if [ $((c4 + c4)) -gt $((call + 10)) ]
then
    echo $ad4
    exit 0
fi

if [ $((c3 + c3)) -gt $((call + 10)) ]
then
    echo $ad3
    exit 0
fi

if [ $((c2 + c2)) -gt $((call + 10)) ]
then
    echo $ad2
    exit 0
else
    # default on standard illumina adapter
    echo $ad1
    exit 0
fi

