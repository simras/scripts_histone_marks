#!/usr/bin/python
#
#
# Simon H. Rasmussen, PhD
#

import sys

def make_feature_profile(fn1,fn2,flip):
    feature_annotation_file = open(fn1)
    feature_score_file = open(fn2)
    pot_first_exon = []
    first_exon = []
    pot_first = False
    old_chr = ""
    score_dic = {}
    ex_dic = {}
    for line in feature_annotation_file:
        line_list =  line.split()

        # load bed info
        my_chr = line_list[0]
        my_st = int(line_list[1])
        my_end = int(line_list[2])
        my_ID = line_list[3]
        expression = line_list[4]
        my_str = line_list[5].strip()
        genomic_coverage = [0] * (my_end - my_st)
        gen_cov_ID = my_ID
        
        # put in genomic score dictionary
        score_dic[gen_cov_ID] = genomic_coverage 
        ex_dic[gen_cov_ID] = expression 

    for f_line in feature_score_file:
        # 
        line_list =  f_line.split()
        
        # load bed info
        my_chr = line_list[0]
        my_st = int(line_list[1])
        my_end = int(line_list[2])
        my_score = int(float(line_list[4].strip()))
        my_ID = line_list[3]        
        my_str = line_list[5].strip()
        gen_cov_ID = my_ID
        
        # manipulate genomic coverage score
        genomic_coverage = score_dic[gen_cov_ID]
        ll=gen_cov_ID.split("_")
        if len(ll) == 6:
                (gene_ID, site_type, ann_chr, ann_st, ann_end, ann_str) = gen_cov_ID.split("_")
        else:
            print >> sys.stderr, ll
        for i in range(my_st - int(ann_st), my_end - int(ann_st)):
            genomic_coverage[i] =  genomic_coverage[i] + my_score
    for k in score_dic.keys():
        expression = ex_dic[k]                
        if flip == "1":
            ll = k.split("_")
            if len(ll) == 6:
                (gene_ID, site_type, ann_chr, ann_st, ann_end, ann_str) = ll
            else:
                print >> sys.stderr, ll
            if ann_str == "-":
                print "\t".join([k ,str(expression)]),
                print "\t".join((str(e) for e in (score_dic[k])[::-1]))
            else:
                print "\t".join([k ,(expression)]),
                print "\t".join((str(e) for e in (score_dic[k])))
        else:
            print "\t".join([k ,(expression)]),
            print "\t".join((str(e) for e in score_dic[k]))
        
make_feature_profile(sys.argv[1], sys.argv[2],sys.argv[3])


