#/bin/bash
#
#
# Make pooled dataset 

# 
rm gene_expression_allpNETseq_TPM.txt

# 
cat GRO_overlap_SRR6661081_gene_ex.bed GRO_overlap_SRR6661082_gene_ex.bed GRO_overlap_SRR6661083_gene_ex.bed GRO_overlap_SRR6661084_gene_ex.bed GRO_overlap_SRR6661085_gene_ex.bed GRO_overlap_SRR6661086_gene_ex.bed GRO_overlap_SRR6661087_gene_ex.bed GRO_overlap_SRR6661088_gene_ex.bed >  GRO_overlap_allpNETseq_gene_ex.bed

