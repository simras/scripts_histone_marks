#!/usr/bin/python
# -*- coding: latin-1 -*-
#
#

def extract_gene_ID(meta_info):
    for l in meta_info.split(";"):
        ll = l.split("=")
        if ll[0] == "ID":
            this_id = ll[1]
            return this_id.split(":")[0]
        elif ll[0] == "Name":
            this_id = ll[1]
    return this_id.split(":")[0]
            

def mk_anno(expr_file, f2_name,cols):
    '''
    Extracts annotations from a gff3 file
    '''
    f1 = open(expr_file)
    f2 = open(f2_name)
    current_transcript = ""
    current_gene = ""
    exon_nr = 0
    boundary_site_dic = {}
    ex_dic = {}
    i = 0
    for l in f1:
        fields = l.split()
        if "gene_ID" == fields[0]:
            i = i + 1
            continue

        # assignment
        ex_ID = fields[0]
        vals = []
        val_mean = 0 
        for c in cols.split(","):
            val_mean = val_mean + float(fields[int(c)])
            vals.append(float(fields[int(c)]))

        # Possibly calc mean
        val_mean = val_mean/len(vals)

        ex_dic[ex_ID] = val_mean

    for l in f2:
        fields = l.split()

        # assignment
        chrom = fields[0]
        start = (fields[1])
        end = (fields[2])

        # gene ID
        ID = fields[3]

        # transcript ID
        ID = fields[3].split("_")[0].split(".")[0]
        
        strand = fields[5]
        try:
            expression = "{0:.2f}".format(ex_dic[ID])
            print "\t".join([chrom, start, end, fields[3], expression, strand])
        except:
            print "\t".join([chrom, start, end, fields[3], "NA", strand])
# Regions
# gene, 3utr, 5utr, cds, all, ncrna, lncrna


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option("-1", action="store", type="string", dest="f1",default="", help="expression file")
    parser.add_option("-2", action="store", type="string", dest="f2",default="", help="Annotation file")
    parser.add_option("-c", action="store", type="string", dest="cols",default="", help="Columns to consider (comma-separated list)")

    (options, args) = parser.parse_args()

    mk_anno(options.f1, options.f2,options.cols)
