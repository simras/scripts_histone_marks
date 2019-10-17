#!/usr/bin/python
# -*- coding: latin-1 -*-
#
#
# 197160 CDS
#      7 chromosome
# 215909 exon
#  34621 five_prime_UTR
#  28775 gene
#    180 miRNA
#  35386 mRNA
#   3911 mRNA_TE_gene
#    480 ncRNA
#  35386 protein
#    924 pseudogene
#   1274 pseudogenic_exon
#    926 pseudogenic_transcript
#     15 rRNA
#     71 snoRNA
#     13 snRNA
#  30634 three_prime_UTR
#   3903 transposable_element_gene
#    689 tRNA
import sys
def extract_gene_ID(meta_info):
    for l in meta_info.split(";"):
        ll = l.split("=")
        if ll[0] == "ID":
            this_id = ll[1]
            return this_id.split(":")[0]
        elif ll[0] == "Name":
            this_id = ll[1]
    return this_id.split(":")[0]
            
def tpm_const(len_dic,read_count):
    a = 0
    for ID,l in len_dic.items():
        a = a + (read_count[ID] * 1000)/l
    return a

def mk_anno(ann_file, f2_name):
    '''
    Extracts annotations from a gff3 file.
    '''
    f1 = open(ann_file)
    f2 = open(f2_name)
    current_transcript = ""
    current_gene = ""
    exon_nr = 0
#    boundary_site_dic = {}
    len_dic = {}
    ex_dic = {}
    read_dic = {}
    i = 0
   # gID="AT1G28450.1"
    for l in f1:
        fields = l.split()
        #if i == 0:
        #    i = i + 1
        #    continue
        
        # assignment
        chrom = fields[0]
        st = float(fields[1])
        end = float(fields[2])
        strand = (fields[5])
        gene_ID = (fields[3]).split("_")[0]
    #    print gene_ID
    #    if gID.split(".")[0] in gene_ID:
    #        print l
        len_dic[gene_ID] = end - st
        ex_dic[gene_ID] = 0
        read_dic[gene_ID] = {}

        # print "1",gene_ID
        # print gene_ID.split("_")[0]
    ii = 0
    #print "Done"
    for l in f2:
        read_fields = (l.split())[0:6]
        ann_fields = (l.split())[12:18]
        
        #if ii % 4000000 == 0:
        #    print ii/183940931.0
 #           break      
        # annotation
        ann_chrom = ann_fields[0]
        ann_start = (ann_fields[1])
        ann_end = (ann_fields[2])
        ann_ID = ann_fields[3].split("_")[0]
        ann_strand = ann_fields[5]

        # Read
        read_chrom = read_fields[0]
        read_start = (read_fields[1])
        read_end = (read_fields[2])
        read_ID = read_fields[3].split("/")[0]
        read_strand = read_fields[5]
  #      print "2",ann_ID
        #print ann_chrom,read_chrom
        # Count
        if not read_dic[ann_ID].has_key(read_ID):
            ex_dic[ann_ID] = ex_dic[ann_ID] + 1
            (read_dic[ann_ID])[read_ID] = ""
        ii = ii + 1
        
 #   sys.exit("Error message")
    del read_dic
    
    Asum = 1000000/tpm_const(len_dic,ex_dic)
    
    for ID,c in ex_dic.items():
        l = float(len_dic[ID])
        A=(c * 1000)/l
        TPM = A*Asum
        print "\t".join([ID,str(TPM)])
    
        #try:
        #    expression = "{0:.2f}".format(ex_dic[ID])
        #    print "\t".join([chrom, start, end, fields[3], expression, strand])
        #except:
        #    print "\t".join([chrom, start, end, fields[3], "NA", strand])
# Regions
# gene, 3utr, 5utr, cds, all, ncrna, lncrna


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option("-1", action="store", type="string", dest="f1",default="", help="Annotation file")
    parser.add_option("-2", action="store", type="string", dest="f2",default="", help="overlap file")

    (options, args) = parser.parse_args()

    mk_anno(options.f1, options.f2)
