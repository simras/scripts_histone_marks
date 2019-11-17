#!/usr/bin/python
# -*- coding: latin-1 -*-
#
#

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
    len_dic = {}
    ex_dic = {}
    read_dic = {}
    i = 0   
    for l in f1:
        fields = l.split()
        
        # assignment
        chrom = fields[0]
        st = float(fields[1])
        end = float(fields[2])
        strand = (fields[5])
        gene_ID = (fields[3]).split("_")[0]
        len_dic[gene_ID] = end - st
        ex_dic[gene_ID] = 0
        read_dic[gene_ID] = {}
    ii = 0
    
    for l in f2:
        read_fields = (l.split())[0:6]
        ann_fields = (l.split())[12:18]
         
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

        # Count
        if not read_dic[ann_ID].has_key(read_ID):
            ex_dic[ann_ID] = ex_dic[ann_ID] + 1
            (read_dic[ann_ID])[read_ID] = ""
        ii = ii + 1
        
    del read_dic
    
    Asum = 1000000/tpm_const(len_dic,ex_dic)
    
    for ID,c in ex_dic.items():
        l = float(len_dic[ID])
        A=(c * 1000)/l
        TPM = A*Asum
        print "\t".join([ID,str(TPM)])
    


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option("-1", action="store", type="string", dest="f1",default="", help="Annotation file")
    parser.add_option("-2", action="store", type="string", dest="f2",default="", help="overlap file")

    (options, args) = parser.parse_args()

    mk_anno(options.f1, options.f2)
