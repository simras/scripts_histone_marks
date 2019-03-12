#!/usr/bin/python
#
#
#

def trim(infile,trim_five,trim_three,ofile):
    import gzip
    o = open(ofile,"w")
    f = gzip.open(infile)
    i = 0
    
    for l in f:
        if i % 4 == 0:
            ID = l.split()[0]
        elif i % 4 == 1:
            ll = l.strip()
            o.write("\t".join([ID,ll[0:trim_five],"\n"]))
        if i % 2 == 1:
            ll = l.strip()
            new_l = len(ll) - trim_three
#            print new_l
            if trim_five + trim_three < len(ll):
                print ll[trim_five:new_l]
            else:
                print ""
        else:
                print l,
        i = i + 1  



if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option("-f", action="store", type="string", dest="ffile",default="", help="Input file")
    parser.add_option("-o", action="store", type="string", dest="ofile",default="", help="Input file")
    parser.add_option("-3", action="store", type="int", dest="trim_three",default=0, help="number of bases to trim if adapter has been removed")
    parser.add_option("-5", action="store", type="int", dest="trim_five",default=0, help="number of bases to trim if adapter has been removed")
 #   parser.add_option("-3", action="store", type="string", dest="ffile",default="", help="Termination")
#    parser.add_option("-r", action="store", type="string", dest="region",default="", help="Annotation region")

    (options, args) = parser.parse_args()

    trim(options.ffile, options.trim_five, options.trim_three,options.ofile)
