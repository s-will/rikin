#!/usr/bin/env Rscript
#
# draw dot plots for all states in a rikin (i)pp file
#

args <- commandArgs(trailingOnly=T)

if (length(args)!=2) {
    cat("USAGE: plot_pps.R pp-file output.pdf\n");
    quit()
}

infile=args[1];
outfile=args[2];

lenA=81
lenB=96

pdf(outfile);

## iterate over lines of infile
##  + plot dot plot for each line

plot(c(),c(),main="",xlab="A",ylab="B",xlim=c(1,lenA),ylim=c(1,lenB))

lines<-readLines(infile)
lines<-lines[c(1:3,5,6)]

indices=c()

color=1
colors=c()

for (line in lines) {
    entries <- unlist( strsplit(line,"[ \t]") );
    index<-as.integer(entries[1]);
    pf<-as.numeric(entries[2]);
    entries<-entries[3:length(entries)]
    indices<-c(indices, index)
    
    for ( k in seq(1,length(entries),3)) {
        i = as.numeric(entries[k+0]);
        j = as.numeric(entries[k+1]);
        p = sqrt(as.numeric(entries[k+2]));
        points(i,j,cex=p,col=color,pch=16);
    }
    colors=c(colors,color)
    color=color+1
}

legend("topright",legend=indices+1,fill=colors,border=F)

dev.off();
