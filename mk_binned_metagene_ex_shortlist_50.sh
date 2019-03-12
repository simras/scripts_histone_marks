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
feature_ann="genome/araport_"$region"_ex.bed"
feature_ann_p="genome/araport_"$region"_ex_p.bed"
feature_ann_m="genome/araport_"$region"_ex_m.bed"

./extract_TSS_term_mid_geneID.py genome/araport_annotation_gene_reptrans.bed $region 0 0|sort|uniq > genome/tmp.bed
python add_expression.py -1 GRO_all.TPM.txt -2 genome/tmp.bed -c "1" > $feature_ann

# Divide on strands
cat $feature_ann|awk '{if($6 == "-"){print $0}}' > $feature_ann_m
cat $feature_ann|awk '{if($6 == "+"){print $0}}' > $feature_ann_p
    
expr="NET"
rm genome/tmp.bed
bins=100
stranded="1"

# shortlist_plot.txt contains run IDs in first colunmn
#for SRR in $(cut -f 1 shortlist_plot.txt)
for SRR in TKR181200217 TKR181200218 SRR6661083 SRR6661084 SRR6661082 SRR6661081  SRR6661085 SRR6661086 SRR6661087 SRR6661088
do

    if [ $stranded ]
    then
	dataset=$SRR"_-+_treat_afterfiting_all.bdg.gz"
	dataset_m=$SRR"_-_treat_afterfiting_all.bdg.gz"
	dataset_p=$SRR"_+_treat_afterfiting_all.bdg.gz"
    else
	dataset=$SRR"_treat_afterfiting_all.bdg.gz"
    fi
    d=bedgraphs/$dataset
    
    y=${d%%_*}
    name=${y##*/}
    
    echo $name
    
    nocompress=bedgraphs/${dataset%.gz}
    nocompress_m=bedgraphs/${dataset_m%.gz}
    nocompress_p=bedgraphs/${dataset_p%.gz}

    # run ID to title translation files
    #title=$(cat shortlist_WT.txt|awk -v id=$name '{if (id == $1){print $2} }'|head -1)
    title=$(cat ID_to_title.txt|awk -v id=$name '{if (id == $1){print $2} }'|head -1)
    
    echo $dataset
    # Clean annotation file
    ann=genome/araport_annotation_exon.bed
    
    region="gene"
    feature_ann="genome/araport_"$region"_ex.bed"
    overlap_file_p="overlap_"$name"_"$region"_ex_p.txt"
    overlap_file_m="overlap_"$name"_"$region"_ex_m.txt"
    overlap_file="overlap_"$name"_"$region"_ex.txt"

    if [ ! -f  profiles/"profile_"$bins"_reptrans"$expr"_bins_"$name"_"$region"_ex.txt" ]
    then
	touch  profiles/"profile_"$bins"_reptrans"$expr"_bins_"$name"_"$region"_ex.txt"

	sorted_bam=mapped/$SRR\_TAIR10.fa_sorted.bam
	if [ $stranded ]
	then
	    bedtools genomecov -trackline -bg -ibam $sorted_bam -strand +|gzip -9 > bedgraphs/$dataset_p
	    bedtools genomecov -trackline -bg -ibam $sorted_bam -strand -|gzip -9 > bedgraphs/$dataset_m
	    bedtools genomecov -trackline -bg -ibam $sorted_bam -strand +|gzip -9 > bedgraphs/$dataset
	    bedtools genomecov -bg -ibam $sorted_bam -strand -|awk 'BEGIN {OFS = "\t"}{$4=-$4;print $0}'|gzip -9 >> bedgraphs/$dataset
	    
	fi
	#	bedtools genomecov -trackline -bg -ibam $sorted_bam -strand +             |gzip -9 > bedgraphs/$SRR"_-+_treat_afterfiting_all2.bdg.gz"
	#	bedtools genomecov -bg -ibam $sorted_bam -strand -|awk '{$4=-$4;print $0}'|gzip -9 >> bedgraphs/$SRR"_-+_treat_afterfiting_all2.bdg.gz"
	echo $nocompress
	if [ ! -f  $nocompress ]
	then
	    gunzip bedgraphs/$dataset
	    if [ $stranded ]
	    then
		gunzip bedgraphs/$dataset_m
		gunzip bedgraphs/$dataset_p
	    fi
	fi
	
	if [ $stranded ]
	then
	    bedtools intersect -b $feature_ann_p -a $nocompress_m -wb|awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$4"\t"$10}' > $overlap_file_m
	    bedtools intersect -b $feature_ann_m -a $nocompress_p -wb|awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$4"\t"$10}' > $overlap_file_p
	    cat $overlap_file_m $overlap_file_p > $overlap_file
	else
	    bedtools intersect -b $feature_ann -a $nocompress -wb|awk '{print $1"\t"$2"\t"$3"\t"$8"\t"$4"\t"$10}' > $overlap_file
	fi

	./mk_feature_profile_ex.py $feature_ann $overlap_file 1 >  profiles/"tmpprofile_"$bins"_bins_"$name"_"$region"_ex.txt"
	./binned_profile.py -1 profiles/"tmpprofile_"$bins"_bins_"$name"_"$region"_ex.txt" -b $bins > profiles/"profile_"$bins"_reptrans"$expr"_bins_"$name"_"$region"_ex.txt"

	if [ $stranded ]
	then
	    rm $overlap_file $overlap_file_p $overlap_file_m profiles/"tmpprofile_"$bins"_bins_"$name"_"$region"_ex.txt" 
	else
	    rm $overlap_file profiles/"tmpprofile_"$bins"_bins_"$name"_"$region"_ex.txt"
	fi
    fi
    
    if [ ! -f bedgraphs/$dataset ]
    then
	gzip $nocompress
	if [ $stranded ]
	then
	    gzip $nocompress_p
	    gzip $nocompress_m
	fi
    fi
    
done

if [ $stranded ]
then
    mat_fc="NET_heatmap_matrix_reptrans50_"$expr"_bins_fc.txt"
    mat="NET_heatmap_matrix_reptrans50_"$expr"_bins.txt"
else
    mat_fc="heatmap_matrix_reptrans50_"$expr"_bins_fc.txt"
    mat="heatmap_matrix_reptrans50_"$expr"_bins.txt"
fi

echo -n "" > $mat_fc
echo -n "" > $mat

# plot
#for SRR in  $(cut -f 1 shortlist_plot.txt)
for SRR in  TKR181200217 TKR181200218 SRR6661081 SRR6661082 SRR6661085 SRR6661086 SRR6661087 SRR6661088 SRR6661083 SRR6661084
do
    name=$SRR
    if [ $stranded ]
    then
	title=$(cat ID_to_title.txt|awk -v id=$SRR '{if (id == $1){print $2} }'|head -1)
    else
	title=$(cat shortlist_plot.txt|awk -v id=$SRR '{if (id == $1){print $2} }'|head -1)	
    fi
    echo $dataset 
    # Clean annotation file
    ann=genome/araport_annotation_exon.bed
    
    region="gene"
    if [ ! -f  $title"_"$name"_reptrans"$expr"_bins_"$region"_50_50.pdf" ]
    then
    
	if [ $stranded ]
	then
	    Rscript print_profile_ex_50pct.R 50 50 profiles/"profile_"$bins"_reptrans"$expr"_bins_"$name"_"$region"_ex.txt" $title"_"$name"_reptrans50"$expr"_bins_"$region $title "NET_heatmap_matrix_reptrans50_"$expr
	else
	    Rscript print_profile_ex_50pct.R 50 50 profiles/"profile_"$bins"_reptrans"$expr"_bins_"$name"_"$region"_ex.txt" $title"_"$name"_reptrans50"$expr"_bins_"$region $title "heatmap_matrix_reptrans50_"$expr
	fi
    fi
    
done
