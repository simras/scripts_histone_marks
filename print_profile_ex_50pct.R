#!/usr/bin/Rscript
#
#
# Simon H. Rasmussen, PhD
#

args = commandArgs(trailingOnly=TRUE)

win1 = as.numeric(args[1])
win2 = as.numeric(args[2])
dname = args[3]
outname = args[4]
title = args[5]
hmap_filename = args[6]

dat = read.table(dname, row.names = 1)

c11 = c()
c12 = c()
m1 = c()

c21 = c()
c22 = c()
m2 = c()

c31 = c()
c32 = c()
m3 = c()

c41 = c()
c42 = c()
m4 = c()

q= quantile(dat[! is.na(dat[,1]),1],probs=c(0.5),na.rm=F)[1]

for (i in 2:length(dat[1,])){
    se = sd(dat[dat[,1] <= q & ! is.na(dat[,1]),i])/sqrt(length((dat[dat[,1] <= q & ! is.na(dat[,1]),i])))
    my_mean = mean(dat[dat[,1] <= q & ! is.na(dat[,1]),i])
    m1 = c(m1,my_mean)
    c11 = c(c11, my_mean - 2* se)
    c12 = c(c12, my_mean + 2* se)
    
}
for (i in 2:length(dat[1,])){
    se = sd(dat[dat[,1] > q & ! is.na(dat[,1]),i])/sqrt(length((dat[dat[,1] > q & ! is.na(dat[,1]),i])))
    my_mean = mean(dat[dat[,1] > q & ! is.na(dat[,1]),i])
    m4 = c(m4,my_mean)
    c41 = c(c41, my_mean - 2* se)
    c42 = c(c42, my_mean + 2* se)

}

c1="darkolivegreen"
c2="deepskyblue3"
c3="yellow"
c4="pink"

pdf(file=paste(outname,"_",win1,"_",win2,".pdf",sep=""),width=16,height=8)
plot(c12 ,type = "l", col="lightgray", lwd=3.5, lty=3, xaxt="n", xlab="Gene Body", ylab="Mean Coverage",ylim=c(0,max(c12,m1,c11,c42,m4,c41)))
lines(m1, col=c1, lwd=2.5)
lines(c11, col="lightgray",lwd=3.5,lty=3)

lines(c42, main=title , col="lightgray", lwd=3.5, lty=3, xaxt="n")
lines(m4, col=c4, lwd=2.5)
lines(c41, col="lightgray",lwd=3.5,lty=3)


if(win1 == 50){
    write.table(rbind(c(title,m4)), append = T, file=paste(hmap_filename,"_bins.txt",sep=""), col.names=F, quote=F, row.names=F)
    write.table(rbind(c(title,log(m4/m1))), append = T, file=paste(hmap_filename,"_bins_fc.txt",sep=""), col.names=F, quote=F, row.names=F)
    axis(side=1,at=c(0,win1+win2),labels=c("TSS","Term"),lwd=3)
    axis(side=2,lwd=3)
    title(main=title)
}else{
    write.table(rbind(c(title,log2(c12))),append = T,file=paste(hmap_filename,"_",win1,"_",win2,".txt",sep=""),col.names=F,quote=F,row.names=F)
    abline(v=win1+1,lwd=2.5)
    abline(v=win1,lwd=2.5)
    axis(side=1,at=seq(0,win1+win2,50),labels=seq(-win1,win2+1,50),lwd=3)
    axis(side=2,lwd=3)
    title(main=title) 
}
legend(y=3*max(c12,m1,c11,c42,m4,c41)/4,x=(3/4)*(win1+win2),fill=c(c4,c1),legend=c("Highest 50%","Lowest 50%"))
dev.off()


