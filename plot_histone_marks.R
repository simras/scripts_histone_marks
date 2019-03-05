#!/usr/local/bin/Rscript
#                                        #
#                                        #
                                        #
require("oompaBase")
zscore = function(dat){
    dd=as.numeric(dat)
    return((dd - mean(dd))/sd(dd))
}
#white_list = c(1, 2, 19,20,21,22,23,24,26,27,28,30,33,34,35,36,37,39,40,41,42)
                                        #white_list = c("","","","","","","","","","","","","","","","","","","","","")
win1 = 50
win2 = 50
white_list = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
methylation_list = c(2,3,10)
acetylation_list = c(11,12,13,14,15,16,17,18)
black_list = c(-9,-10,-11,-12,-13,-14,-15,-16)
#d = read.csv("heatmap_matrix_reptransNET_bins_fc.txt",row.names=1,sep=" ",header=F)[white_list, ]
#d = read.csv("heatmap_matrix_reptransNET_bins.txt",row.names=1,sep=" ",header=F)[white_list, ]

#d = read.csv("heatmap_matrix_reptransRNAseq_bins.txt",row.names=1,sep=" ",header=F)[white_list, ]
#d = read.csv("heatmap_matrix_reptransRNAseq_bins_fc.txt",row.names=1,sep=" ",header=F)[white_list, ]

#d = read.csv("heatmap_matrix_reptransoldRNAseq_bins_fc.txt",row.names=1,sep=" ",header=F)[white_list, ]
d = read.csv("heatmap_matrix_reptrans50_NET_bins.txt",row.names=1,sep=" ",header=F)[white_list, ]
colnames(d) = 1:100
# Rainbow color scheme
cc = rainbow(n=length(d[,1]))
#cc = c("darkred","darkblue","darkolivegreen4","yellow","purple","black")
d3 = d
#d3 = rbind(ac_marks,d)
#d3 = rbind(mt_marks,d3)
#d = d[c(1,2,3,6,7,9,10),]
cc3 = jetColors(length(d[,1]))
#rownames(d2)[1] = "Methylation Marks"
#rownames(d2)[2] = "Acetylation Marks"


ac_marks = apply(d[acetylation_list,],2,mean)
mt_marks = apply(d[methylation_list,],2,mean)

d2 = rbind(ac_marks,d)
d2 = rbind(mt_marks,d2)
d2 = d2[c(1,2,3,6,7,9,10),]
cc2 = jetColors(length(d2[,1]))
rownames(d2)[1] = "Methylation Marks"
rownames(d2)[2] = "Acetylation Marks"

print(d)
d_no_rep = d[c(1,6,7,8,9,10,11,12,13,14,15,16,17,18),]
print(d)
d1 = rbind(ac_marks,d_no_rep)
#d1 = rbind(mt_marks,d1)
d1 = d1[1:7,]
cc1 = jetColors(length(d1[,1]))
#rownames(d1)[1] = "Methylation Marks"
rownames(d1)[1] = "Acetylation Marks"


# Plot with one grouped line
#pdf(file = "all_histone_marks_highlyexpressed_NET50.pdf",width=12,height=8)
#plot(zscore(d2[1,]),ylim=c(-5,5),col=cc[1],type="l",lwd=2.5)

#for(i in 2:length(d2[,1])){
#    lines(zscore(d2[i,]),col=cc[i],lwd=2.5)
#}
#legend(x=75,y=5,legend=rownames(d2),lwd=2.5,fill=cc)
#dev.off()
cc_no_rep = jetColors(length(d_no_rep[,1]))
print(zscore(d_no_rep[14,]))
    
pdf(file = "all_histone_marks_highlyexpressed_zscore_NET50_all.pdf",width=12,height=8)
plot(zscore(d3[1,]),ylim=c(-5,5),col=cc3[1],type="l",lwd=2.5,xaxt="n",yaxt="n", xlab="Gene Body", ylab="Relative enrichment")

for(i in 2:length(d3[,1])){
    lines(zscore(d3[i,]),col=cc3[i],lwd=2.5)
}
legend(x=75,y=5,legend=rownames(d3),lwd=2.5,fill=cc3)
dev.off()

pdf(file = "all_histone_marks_highlyexpressed_zscore_NET50_norep.pdf",width=12,height=8)
plot(zscore(d_no_rep[1,]),ylim=c(-5,5),col=cc_no_rep[1],type="l",lwd=2.5,xaxt="n",yaxt="n", xlab="Gene Body", ylab="Relative enrichment")

for(i in 2:length(d_no_rep[,1])){
    lines(zscore(d_no_rep[i,]),col=cc_no_rep[i],lwd=2.5)
}
legend(x=75,y=5,legend=rownames(d_no_rep),lwd=2.5,fill=cc_no_rep)
axis(side=1,at=c(0,win1+win2),labels=c("TSS","Term"),lwd=3)
axis(side=2,lwd=3)
dev.off()

pdf(file = "all_histone_marks_highlyexpressed_zscore_NET50_one_grouped.pdf",width=12,height=8)
plot(zscore(d1[1,]),ylim=c(-5,5),col=cc1[1],type="l",lwd=2.5,xaxt="n",yaxt="n", xlab="Gene Body", ylab="Relative enrichment")

for(i in 2:length(d1[,1])){
    lines(zscore(d1[i,]),col=cc1[i],lwd=2.5)
}
legend(x=75,y=5,legend=rownames(d1),lwd=2.5,fill=cc1)
axis(side=1,at=c(0,win1+win2),labels=c("TSS","Term"),lwd=3)
axis(side=2,lwd=3)
dev.off()

#pdf(file = "all_histone_marks_highlyexpressed_zscore_NET50_two_grouped.pdf",width=12,height=8)
#plot(zscore(d[1,]),ylim=c(-5,5),col=cc[1],type="l",lwd=2.5)

#for(i in 2:length(d[,1])){
#    lines(zscore(d[i,]),col=cc[i],lwd=2.5)
#}
#legend(x=75,y=5,legend=rownames(d),lwd=2.5,fill=cc)
#dev.off()

# plot with two grouped lines
#rownames(d2) = c("Methylation Marks", "Acetylation Marks", "H2AK121ub", "H3K36me2", "H3K36me3", "H3K4me2", "H3K4me3")
#d = d[a,]

#pdf(file = "all_histone_marks_highlyexpressed_NET50.pdf",width=12,height=8)
#plot(as.numeric(d[1,]),ylim=c(min(d),max(d)),col=cc[1],type="l")

#for(i in 2:length(d[,1])){
#    lines(as.numeric(d[i,]),col=cc[i])
#}
#dev.off()


pdf(file = "all_histone_marks_highlyexpressed_zscore_NET50_two_grouped.pdf",width=12,height=8)
plot(zscore(d2[1,]),ylim=c(-5,5),col=cc2[1],type="l",lwd=2.5,yaxt='n',xaxt='n',xlab="",ylab="")

for(i in 2:length(d2[,1])){
    lines(zscore(d2[i,]),col=cc2[i],lwd=2.5)
}
legend(x=75,y=5,legend=rownames(d2),lwd=2.5,fill=cc2)
axis(1,xlab="Gene body")
axis(2)
dev.off()

#exit()
#d = read.csv("heatmap_matrix_reptrans50NET_bins_fc.txt",row.names=1,sep=" ",header=F)[white_list, ]

# Rainbow color scheme
#cc = rainbow(n=length(d[,1]))

#pdf(file = "all_histone_marks_highlyexpressed_NET_fc50.pdf",width=12,height=8)
#par(mar=c(1,1,1,1))
                                        # plotting
#dz = zscore(d)
#plot(as.numeric(d[1,]),ylim=c(min(d),max(d)),col=cc[1],type="l")

#for(i in 2:length(d[,1])){
#    lines(as.numeric(d[i,]),col=cc[i])
#}
#legend(x=75,legend=row.names(d),lwd=2.5,fill=cc)
#dev.off()
#pdf(file = "all_histone_marks_highlyexpressed_zscore_NET_fc50.pdf",width=12,height=8)
#plot(zscore(d[1,]),ylim=c(-5,5),col=cc[1],type="l")

#for(i in 2:length(d[,1])){
#    lines(zscore(d[i,]),col=cc[i])
#}
#dev.off()
