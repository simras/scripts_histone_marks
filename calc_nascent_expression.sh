#/bin/bash
#
#
#
# Simon H. Rasmussen
#

set -e
set -x

win=500
region="gene"
#feature_ann="genome/araport_"$region"_ex.bed"
feature_ann="genome/araport_annotation_gene_reptrans.bed"
#feature_ann=genome/araport_annotation_gene.bed
#./extract_TSS_term_mid.py genome/araport_annotation_gene.bed $region 0 0|sort|uniq > genome/tmp.bed
#python add_expression.py -1 RNA-seq/norm_col0_RNA-seq_counts_mean.txt -2 genome/tmp.bed > $feature_ann

#rm genome/tmp.bed #genome/tmp2.bed
#bins=100
# donor sites allpNETseq 
for SRR in allpNETseq SRR6661081 SRR6661082 SRR6661085 SRR6661086 SRR6661087 SRR6661088 TKR181200217 TKR181200218 SRR6661083 SRR6661084
do
    
    
    #  
    d=mapped/$SRR"_TAIR10.fa_sorted.bam"
    title=$(cat ID_to_title.txt|awk -v id=$SRR '{if (id == $1){print $2} }'|head -1)
    bam_dataset=$d
    
    overlap_file="GRO_overlap_"$SRR"_"$region"_ex.bed"

    #cat GRO_overlap_SRR6661081_gene_ex.bed GRO_overlap_SRR6661082_gene_ex.bed GRO_overlap_SRR6661083_gene_ex.bed GRO_overlap_SRR6661084_gene_ex.bed GRO_overlap_SRR6661085_gene_ex.bed GRO_overlap_SRR6661086_gene_ex.bed GRO_overlap_SRR6661087_gene_ex.bed GRO_overlap_SRR6661088_gene_ex.bed > GRO_overlap_allpNETseq_gene_ex.bed
    
    if [ ! -f $overlap_file ]
    then
	bedtools intersect -S -abam $bam_dataset -b $feature_ann -wb -bed > $overlap_file
#	bedtools intersect -s -abam $bam_dataset -b $feature_ann -wb -bed > $overlap_file
 #   else
#	if [ ! $(wc -l $overlap_file|cut -f 1 -d " ") -gt 0 ]
#	then
#	    bedtools intersect -abam $bam_dataset -b $feature_ann -wb -bed > $overlap_file
#	fi
    fi    
    if [ ! -f "gene_expression_"$SRR"_TPM.txt" ]
    then
	./GRO_seq_expression_readdict.py -1 $feature_ann -2 $overlap_file|sort -k1,1 > "gene_expression_"$SRR\_TPM.txt
 #   else
#	if [ ! $(wc -l "gene_expression_"$SRR"_TPM.txt"|cut -f 1 -d " ") -gt 0 ]
#	then
#	    ./GRO_seq_expression_readdict.py -1 $feature_ann -2 $overlap_file|sort -k1,1 > "gene_expression_"$SRR\_TPM.txt
#	fi
    fi
    #rm $overlap_file
done

#SRR6821567 SRR6821571 SRR6821575 SRR5681049 SRR5681050 SRR5681055 SRR5681056 SRR6661079 SRR6661080
	    #SRR6661081 SRR6661082 SRR6661083 SRR6661084 SRR6661085 SRR6661086 SRR6661087 SRR6661088 TKR181200217 TKR181200218
echo -n "gene_ID" > GRO_all.TPM.txt
for SRR in  allpNETseq SRR6821567 SRR6821571 SRR6821575 SRR5681049 SRR5681050 SRR5681055 SRR5681056 SRR6661079 SRR6661080 SRR6661081 SRR6661082 SRR6661083 SRR6661084 SRR6661085 SRR6661086 SRR6661087 SRR6661088 TKR181200217 TKR181200218
do
    awk -v sid=$SRR 'BEGIN{printf "\tTPM_%s", sid}' >> GRO_all.TPM.txt
   # echo -n "TMP_" >> GRO_all.TMP.txt

done
#SRR6821571 SRR6821575 SRR5681049 SRR5681050 SRR5681055 SRR5681056 
echo "" >> GRO_all.TPM.txt
cat "gene_expression_"allpNETseq\_TPM.txt > tmptmp.here
#for SRR in  SRR6661082 SRR6661083 SRR6661084 SRR6661085 SRR6661086 SRR6661087 SRR6661088 TKR181200217 TKR181200218
for SRR in SRR6821567 SRR6821571 SRR6821575 SRR5681049 SRR5681050 SRR5681055 SRR5681056 SRR6661079 SRR6661080 SRR6661081 SRR6661082 SRR6661083 SRR6661084 SRR6661085 SRR6661086 SRR6661087 SRR6661088 TKR181200217 TKR181200218
do
    cat tmptmp.here > gene_expression_tmp_TPM.txt
    join -j 1 gene_expression_tmp_TPM.txt "gene_expression_"$SRR"_TPM.txt" > tmptmp.here
    #|awk '{print $1"\t"$2"\t"$3"\t"($2+$3)/2}' >> GRO_all.TMP.txt
#    head tmptmp.here
done
cat tmptmp.here > gene_expression_tmp_TPM.txt

#awk 'BEGIN {print "gene_ID\tTPM_SRR5681049\tTMP_SRR5681050"}' > GRO_all.TMP.txt
#join -j 1 "gene_expression_"SRR5681049\_TPM.txt "gene_expression_"SRR5681050\_TPM.txt|awk '{print $1"\t"$2"\t"$3"\t"($2+$3)/2}' >> GRO_all.TMP.txt

cat gene_expression_tmp_TPM.txt|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21}' >> GRO_all.TPM.txt
exit 0
#echo -n "gene_ID" > GRO_all.TMP.txt
#for SRR in  SRR7474123 SRR5681049 SRR5681050 SRR5681055 SRR7474124 SRR7474125
#do
#    awk -v sid=$SRR 'BEGIN{printf "\tTPM_%s", sid}' >> GRO_all.TMP.txt
#done

#cat "gene_expression_"SRR7474123\_TPM.txt > tmptmp.here
#for SRR in  SRR5681049 SRR5681050 SRR5681055 SRR7474124 SRR7474125
#do
#    cat tmptmp.here > gene_expression_tmp_TPM.txt
#    join -j 1 gene_expression_tmp_TPM.txt "gene_expression_"$SRR"_TPM.txt" > tmptmp.here
#done
#cat gene_expression_tmp_TPM.txt|awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' >> GRO_all.TPM.txt
