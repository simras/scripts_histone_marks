#!/bin/bash
#
# Automatic ChIP pipeline
#
# Use:
# for f in data/*.fastq.gz;do ./pick_adapter.sh $f;don

#set -x

ad1=AGATCGGAAGAGC
#ad2=CTGTCTCTTATACACATCT
ad2=GCGGGTTGAGCCAAGCCTGTGCAGTA
#ad3=TGGAATTCTCGG
ad4=GATCGTCGGACTGTAGAAC
ad3=TGGAATTCTCGGGTGCCAA
#ad5=GATCGTCGGACTGTAGAAC
#pa=AAAAAAAAAA
#pt=TTTTTTTTTT

#echo "Adapters"
ll=$(zcat $1|head -2|tail -1|wc -c)
#c4=0
#IFS=""
#if (( "$ll" < 50 ))
#then
#    echo "short reads" $ll #$(zcat $1|head -2|tail -1)
#    ad4=CTGTCTCTTATA

c1=$(zcat $1|head -400000 |awk 'BEGIN{i = 0}{if(i % 4 == 1){print $0};i = i + 1}'|grep -c $ad1)
c2=$(zcat $1|head -400000 |awk 'BEGIN{i = 0}{if(i % 4 == 1){print $0};i = i + 1}'|grep -c $ad2)
c3=$(zcat $1|head -400000 |awk 'BEGIN{i = 0}{if(i % 4 == 1){print $0};i = i + 1}'|grep -c $ad3)
c4=$(zcat $1|head -400000 |awk 'BEGIN{i = 0}{if(i % 4 == 1){print $0};i = i + 1}'|grep -c $ad4)
#c5=$(zcat $1|head -400000 |awk 'BEGIN{i = 0}{if(i % 4 == 1){print $0};i = i + 1}'|grep -c $ad5)
call=$((c1 + c2 + c3 + c4))
#echo "Counts" $c1 $c2 $c3 $c4 
#if [ $((c5 + c5)) -gt $((call + 10)) ]
#then
#    echo $ad5
#        exit 0
#fi

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
    echo $ad1
    exit 0
fi

