#!/usr/bin/python
#
# Simon H. Rasmussen
# 


import sys
import math

def bin_profile(tss,bins):
    # 
    #
    #
    profile = open(tss)
    for l in profile:
        ex_dict = []
        ll = l.split()
        ID = ll[0]
        express = ll[1]
        pf = ll[2:(len(ll))]
        n = len(ll) - 2
        if n < bins:
            continue
        i = 0
        j = 0
        bin_sum = 0
        b = []
        num = int(n/bins)
        for bb in range(bins):
            
            if n % num > 0:
                rest = 1
            else:
                rest = 0
            b.append(num + rest)
            n = n - num - rest
        ii = 0
        bbb = b[ii]
        # Iterate over values of genomic coverage in profile
        for val in pf:
            i = i + 1
            j = j + 1
            bin_sum = bin_sum + float(val)
            # new bin
            if j == bbb:
                ex_dict.append(str(bin_sum/b[ii]))
                ii = ii + 1
                if ii < len(b):
                    bbb = bbb + b[ii]
                i = 0
                bin_sum = 0
        ex_dict = [ID , express] + ex_dict     
        print "\t".join(ex_dict)
        ex_dict = []
        
        
if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option("-1", action="store", type="string", dest="ffile",default="", help="TSS")
    parser.add_option("-b", action="store", type="int", dest="bins",default=200, help="Annotation region")

    (options, args) = parser.parse_args()

    bin_profile(options.ffile,options.bins)
