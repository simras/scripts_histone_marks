#!/usr/bin/python
#
#
# Simon H. Rasmussen

# TIF-seq
# 79842
#
# All
# ./estimate_representative_transcript.py -1 ../Araport11_GFF3_genes_transposons.201606.gff -2 TSS_seq_+.bedgraph -3 TSS_seq_-.bedgraph -4 TIF_seq.bed -5 Wu2011_wt_leaf_merged_fw_rev_corrected.bedgraph -6 Thomas2012_wt_merged_fw_rev_corrected.bedgraph|wc -l
# 31556
#
# unique start
# ./estimate_representative_transcript.py -1 ../Araport11_GFF3_genes_transposons.201606.gff -2 TSS_seq_+.bedgraph -3 TSS_seq_-.bedgraph -4 TIF_seq.bed -5 Wu2011_wt_leaf_merged_fw_rev_corrected.bedgraph -6 Thomas2012_wt_merged_fw_rev_corrected.bedgraph|awk '{if($7 == 0){print $7,$8}}'|wc -l
#
# 22889
# Multiple starts
# ./estimate_representative_transcript.py -1 ../Araport11_GFF3_genes_transposons.201606.gff -2 TSS_seq_+.bedgraph -3 TSS_seq_-.bedgraph -4 TIF_seq.bed -5 Wu2011_wt_leaf_merged_fw_rev_corrected.bedgraph -6 Thomas2012_wt_merged_fw_rev_corrected.bedgraph|awk '{if($7 > 0){print $7,$8}}'|wc -l
# 8667
#
# Unique end
# ./estimate_representative_transcript.py -1 ../Araport11_GFF3_genes_transposons.201606.gff -2 TSS_seq_+.bedgraph -3 TSS_seq_-.bedgraph -4 TIF_seq.bed -5 Wu2011_wt_leaf_merged_fw_rev_corrected.bedgraph -6 Thomas2012_wt_merged_fw_rev_corrected.bedgraph|awk '{if($8 == 0){print $7,$8}}'|wc -l
# 22851
#
# Multiple ends
# ./estimate_representative_transcript.py -1 ../Araport11_GFF3_genes_transposons.201606.gff -2 TSS_seq_+.bedgraph -3 TSS_seq_-.bedgraph -4 TIF_seq.bed -5 Wu2011_wt_leaf_merged_fw_rev_corrected.bedgraph -6 Thomas2012_wt_merged_fw_rev_corrected.bedgraph|awk '{if($8 > 0){print $7,$8}}'|wc -l
# 8705

import numpy as np
import copy

def mk_chrom_arrays():
    chrom_hash = {}
    l = [34964571 + 10000, 22037565 + 10000, 25499034 + 10000 ,20862711 + 10000, 31270811 + 10000,371576,156441]
    chrom_hash["Chr1"] = (np.zeros(l[0])).astype(float)
    chrom_hash["Chr2"] = (np.zeros(l[1])).astype(float)
    chrom_hash["Chr3"] = (np.zeros(l[2])).astype(float)
    chrom_hash["Chr4"] = (np.zeros(l[3])).astype(float)
    chrom_hash["Chr5"] = (np.zeros(l[4])).astype(float)
    chrom_hash["ChrM"] = (np.zeros(l[5])).astype(float)
    chrom_hash["ChrC"] = (np.zeros(l[6])).astype(float)
    return chrom_hash

def write_bedgraph(dicp,header,region):
    if header != "":
        print header
    for chrom in ["Chr1","Chr2","Chr3","Chr4","Chr5","ChrM","ChrC"]:
        i = 0
        for dp in dicp[chrom]:
            if not dp == 0:
                if region == "start":
                    print "\t".join([chrom, str(i), str(i+1), str(int(dp))])
                else:
                    print "\t".join([chrom, str(i), str(i+1), str(int(-1*dp))])
            i = i + 1
    
def norm_data(dat):
    # Number of reads and loci
    loci = 0
    R_all = 0
    for datapoint in dat:
        if datapoint[0] == "t":
            continue
        dplist = datapoint.split()
        chrom = dplist[0]
        start = int(dplist[1])
        end = int(dplist[2])
        try:
            strand = dplist[5]
        except:
            strand = ""
            None
        if strand != "":
            # bedfile
            val = abs(float(dplist[4]))                
        else:
            val = float(dplist[3])
        R_all = R_all + abs(val)
        loci = loci + 1
    return (float(R_all), loci)

def read_data(dat, start_dicp, start_dicm, end_dicp, end_dicm):
    # TSS-seq plus
    R, loci = norm_data(dat)
    ep = R/loci
    em = R/loci
    for datapoint in dat:
        if datapoint[0] == "t":
            continue
        dplist = datapoint.split()
        chrom = dplist[0]
        start = int(dplist[1])
        end = int(dplist[2])
        try:
            strand = dplist[5]
        except:
            strand = ""
            None
        if strand != "":
            # TIF-seq
            val = abs(float(dplist[4]))
            if strand == "+":
                (start_dicp[correct_chrom(chrom)])[start] = (start_dicp[correct_chrom(chrom)])[start] + val/ep
                (end_dicp[correct_chrom(chrom)])[end] = (end_dicp[correct_chrom(chrom)])[end] + val/ep
            else:
                (start_dicm[correct_chrom(chrom)])[end] = (start_dicm[correct_chrom(chrom)])[end] + val/em
                (end_dicm[correct_chrom(chrom)])[start] = (end_dicm[correct_chrom(chrom)])[start] + val/em
        else:
            val = (float(dplist[3]))
            # not TIF-seq
        if not end_dicp and not end_dicm:
            # TSS
            if val > 0:
                if not start_dicp:
                    (start_dicm[correct_chrom(chrom)])[start] = (start_dicm[correct_chrom(chrom)])[start] + abs(val)/em
                else:
                    (start_dicp[correct_chrom(chrom)])[start] = (start_dicp[correct_chrom(chrom)])[start] + abs(val)/ep
            else:
                (start_dicm[correct_chrom(chrom)])[start] = (start_dicm[correct_chrom(chrom)])[start] + abs(val)/em
        elif not start_dicp  and not start_dicm :
            # TTS
            if val > 0:
                (end_dicp[correct_chrom(chrom)])[start] = (end_dicp[correct_chrom(chrom)])[start] + abs(val)/ep
            else:
                (end_dicm[correct_chrom(chrom)])[start] = (end_dicm[correct_chrom(chrom)])[start] + abs(val)/em
    return (start_dicp, start_dicm, end_dicp, end_dicm)

def correct_chrom(chrom):
    trans_table=["Chr1","Chr2","Chr3","Chr4","Chr5","ChrM","ChrC"]
    #trans_table=["chr1","chr2","chr3","chr4","chr5","chrM","chrC"]
    if chrom == "Chr1" or chrom == "chr1" or chrom == "1":
        return trans_table[0]
    elif chrom == "Chr2" or chrom == "chr2" or chrom == "2" :
        return trans_table[1]
    elif chrom == "Chr3" or chrom == "chr3" or chrom == "3" :
        return trans_table[2]
    elif chrom == "Chr4" or chrom == "chr4" or chrom == "4" :
        return trans_table[3]
    elif chrom == "Chr5" or chrom == "chr5" or chrom == "5" :
        return trans_table[4]
    elif chrom == "ChrM" or chrom == "chrMt" or chrom == "Mt" :
        return trans_table[5]
    elif chrom == "ChrC" or chrom == "chrPt" or chrom == "Pt" :
        return trans_table[6]
    else:
        exit
        
def output_chrom(chrom):
    trans_table=["chr1","chr2","chr3","chr4","chr5","chrMT","chrCP"]
    #trans_table=["chr1","chr2","chr3","chr4","chr5","chrM","chrC"]
    if chrom == "Chr1" or chrom == "chr1" or chrom == "1":
        return trans_table[0]
    elif chrom == "Chr2" or chrom == "chr2" or chrom == "2" :
        return trans_table[1]
    elif chrom == "Chr3" or chrom == "chr3" or chrom == "3" :
        return trans_table[2]
    elif chrom == "Chr4" or chrom == "chr4" or chrom == "4" :
        return trans_table[3]
    elif chrom == "Chr5" or chrom == "chr5" or chrom == "5" :
        return trans_table[4]
    elif chrom == "ChrM" or chrom == "chrMt" or chrom == "Mt" :
        return trans_table[5]
    elif chrom == "ChrC" or chrom == "chrPt" or chrom == "Pt" :
        return trans_table[6]
    else:
        exit

def find_extremes(transcript_hash):
    min_st = 999999999
    max_st = -999999999
    min_end = 9999999999
    max_end = -9999999999
    extreme_list_init = [min_st, max_st, min_end, max_end, "", ""]
    extreme_hash = {} 
    for ID,genes in transcript_hash.items():
        extreme_list = copy.copy(extreme_list_init)
        for transcript in genes:
            start = transcript[1]
            end = transcript[2]
            strand = transcript[3]
            chrom = transcript[0]
            if start < extreme_list[0]:
                extreme_list[0] = start
            if start >= extreme_list[1]:
                extreme_list[1] = start
            if end < extreme_list[2]:
                extreme_list[2] = end
            if end >= extreme_list[3]:
                extreme_list[3] = end
            extreme_list[4] = strand
            extreme_list[5] = chrom
        extreme_hash[ID] = extreme_list
    return extreme_hash

def find_max_point(data,start,end,w):
    #
    
    if start == end:
        return (start,0)
    max_dp = -999
    max_i = -1
    i = 0
    ties = []
    for dp in data[(start - w):(end + w)]:
        
        if dp > max_dp:
            max_dp = dp
            max_i = i
        elif dp == max_dp:
            ties.append(max_i)
            max_i = i
        i = i + 1
    return (start + max_i + 1 - w,max_dp)
            
def overlay_data(ID,chrom, start_region, end_region, strand, data_hashes, w):
    #
    start_array_p = data_hashes[0]
    start_array_m = data_hashes[1]

    # 
    end_array_p = data_hashes[2]
    end_array_m = data_hashes[3]

    # find TSS max in start region
    if strand == "+":
        new_start, score_start = find_max_point((start_array_p[chrom]), start_region[0], start_region[1],w)
        new_end, score_end     = find_max_point((end_array_p[chrom]), end_region[0], end_region[1],w)
    else:
        new_start, score_start = find_max_point((end_array_m[chrom]), start_region[0], start_region[1],w)
        new_end, score_end = find_max_point((start_array_m[chrom]), end_region[0], end_region[1],w)
        
    if score_start == 0:
        new_start = start_region[0]
    else:
        new_start = new_start - 1

    if score_end == 0:
        new_end = end_region[1]
  
    if score_start > 0 and score_end > 0:
        return (new_start,new_end,(score_start + score_end)/2)
    else:
        return (new_start,new_end,max(score_start, score_end))
    
    
    
def find_new_boundaries(genes_hash, data_hashes,n,w):
    
    for ID, transcript in genes_hash.items():
        # Expand by w

        min_start = transcript[0]# - w
        max_start = transcript[1] #+ w # may break the protein
        min_end = transcript[2] #- w # may break the protein
        max_end = transcript[3]# + w
        strand = transcript[4]
        chrom = transcript[5]
        
        start_region =  max_start - min_start
        end_region = max_end - min_end
        max_length = max_end - min_start
        min_length = min_end - max_start
        
        new_start, new_end,score = overlay_data(ID,chrom, [min_start, max_start], [min_end, max_end], strand, data_hashes,w)
        # compute statistics
        # 
        if strand == "-":
            print "HERE",new_start - min_start, new_end - max_end
        if not n:
            if new_end-new_start > 10:
                print "\t".join([output_chrom(chrom), str(new_start), str(new_end), ID, str(score), strand])
            

def extract_tID(mInfo):
    ilist = mInfo.split(";")
    for it in ilist:
        if it[0:3] == "ID=":
            ID = it[3:]
    return ID

import sys
def estimate(ann_file,TSSp,TSSm,TIF,TTS1,TTS2,n,full,w):
    # headers
    h_TSS_m = "track type=bedGraph name=S07_Col-0-rep2_rev description=S07_Col-0-rep2_rev visibility=full color=200,100,0 autoScale=off viewLimits=0:60 maxHeightPixels=60"
    h_TSS_p = "track type=bedGraph name=S07_Col-0-rep2_fw description=S07_Col-0-rep2_fw visibility=full color=0,100,200 autoScale=off viewLimits=0:60 maxHeightPixels=60"

    h_TIF = "track name=TIF_seq type=bedGraph color=0,100,200 altColor=200,100,0" 

    h_TTS_1 = "track name=TTS_WU2011 type=bedGraph color=0,100,200 altColor=200,100,0"
    h_TTS_2 = "track name=TTS_Thomas2012 type=bedGraph color=0,100,200 altColor=200,100,0"

    header_all = "track name=All_data type=bedGraph description=TSS_TIF_TTS visibility=full color=200,100,0 autoScale=off viewLimits=0:60 maxHeightPixels=60"
    
    start_array_p = mk_chrom_arrays()
    start_array_m = mk_chrom_arrays()
    end_array_p = mk_chrom_arrays()
    end_array_m = mk_chrom_arrays()

    f_TSS_p = open(TSSp)
    TSS_seq_p = f_TSS_p.readlines()
    f_TSS_p.close()
    
    start_array_p, tmp, tmp, tmp = read_data(TSS_seq_p, start_array_p, {}, {}, {})

    f_TSS_m = open(TSSm)
    TSS_seq_m = f_TSS_m.readlines()
    f_TSS_m.close()
    
    tmp, start_array_m, tmp, tmp = read_data(TSS_seq_m, {}, start_array_m, {}, {})

    TIF_seq = open(TIF)
    TIF_seq_data = TIF_seq.readlines()
    
    start_array_p, start_array_m, end_array_p, end_array_m = read_data(TIF_seq_data, start_array_p, start_array_m, end_array_p, end_array_m)

    f_TTS_1 = open(TTS1)
    TTS_seq_1 = f_TTS_1.readlines()
    
    tmp, tmp, end_array_p, end_array_m = read_data(TTS_seq_1, {}, {}, end_array_p, end_array_m)

    f_TTS_2 = open(TTS2)
    TTS_seq_2 = f_TTS_2.readlines()
    
    tmp, tmp, end_array_p, end_array_m = read_data(TTS_seq_2, {}, {}, end_array_p, end_array_m)

    dat_list = [start_array_p, start_array_m, end_array_p, end_array_m]
    gene_annotation = open(ann_file)
    transcript_hash = {}
    extreme_hash = {}
    
    f_TTS_1.close()
    f_TTS_2.close()
    TIF_seq.close()
    
    # Analyze annotation
    for l in gene_annotation:
        if l[0] == "#":
            continue
        line_list = l.split()
        ann_chr = line_list[0]
        ann_feat_name = line_list[2]
        ann_st = int(line_list[3])
        ann_end = int(line_list[4])
        ann_str = line_list[6]
        mInfo = line_list[8]
        
        ann_transcript_ID = extract_tID(mInfo)
        ann_gene_ID = ann_transcript_ID.split(".")[0]
        
        if ann_feat_name == "mRNA" or ann_feat_name == "protein":
            if transcript_hash.has_key(ann_gene_ID):
                transcript_hash[ann_gene_ID].append([ann_chr,ann_st,ann_end,ann_str])
            else:
                transcript_hash[ann_gene_ID] = [[ann_chr,ann_st,ann_end,ann_str]]
            if full:
                print "\t".join([ann_chr,str(ann_st),str(ann_end),ann_transcript_ID,".",ann_str])
    transcript_extremes = find_extremes(transcript_hash)
    find_new_boundaries(transcript_extremes,dat_list,n,w)

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option("-1", action="store", type="string", dest="f1",default="", help="Annotation file"       )
    parser.add_option("-2", action="store", type="string", dest="f2",default="", help="TSS+"                  )
    parser.add_option("-3", action="store", type="string", dest="f3",default="", help="TSS-"                  )
    parser.add_option("-4", action="store", type="string", dest="f4",default="", help="TIF"                   )
    parser.add_option("-5", action="store", type="string", dest="f5",default="", help="TTS"                   )
    parser.add_option("-6", action="store", type="string", dest="f6",default="", help="TTS"                   )
    parser.add_option("-w", action="store", type="int", dest="w",default=0, help="Window to consider around min-max boundaries.")
    parser.add_option("-n" , action="store_true", dest="nn" ,default=False, help="Don't print annotation")
    parser.add_option("-f" , action="store_true", dest="full" ,default=False, help="Don't print annotation")

    (options, args) = parser.parse_args()

    estimate(options.f1, options.f2, options.f3, options.f4, options.f5, options.f6, options.nn, options.full,options.w)
