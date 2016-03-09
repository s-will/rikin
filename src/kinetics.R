#!/usr/bin/Rscript

### parse command line using python-like argparse

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()

parser$add_argument("-o", "--output", default="kinetics.pdf",
                    help="output pdf file [default:\"%(default)s\"]")
parser$add_argument("--info", default="",
                    help="text to be shown in the plot [optional]")
parser$add_argument("--pps", default="",
                    help="file with pair probabilities; if given, add dot plot  [optional]")
parser$add_argument("--seqA", default="",
                    help="sequence A [required for --pps]")
parser$add_argument("--seqB", default="",
                    help="sequence B [required for --pps]")
parser$add_argument("--nameA", default="",
                    help="sequence A name")
parser$add_argument("--nameB", default="",
                    help="sequence B name")
parser$add_argument("--restrict-dp", action="store_true", default=FALSE, help="restrict dot plot to region of probabilities above threshold [default: off]")

parser$add_argument("--evalplot", action="store_true", default=FALSE, help="draw evaluation plot [default: off]")

parser$add_argument("input",nargs=1,help="input file [required]")

args <- parser$parse_args()

input    <- args$input
output   <- args$output
ppfile   <- args$pps
info     <- args$info
pps      <- args$pps
seqA     <- args$seqA
seqB     <- args$seqB
nameA     <- args$nameA
nameB     <- args$nameB
evalplot <- args$evalplot
restrictdp <- args$restrict_dp

if (pps!="") {
    if (seqA=="" || seqB=="") {
        cat("error: missing sequence(s) [required for --pps]\n");
        quit();
    }
}


if ( !file.exists(input) ) {
    inputgz <- paste(sep="",input,".gz");
    if (file.exists(inputgz)) {
        input<-inputgz
    } else {
        cat("error: input file does not exist\n");
        quit();
    }
}

### end parse command line

############################################################
### functions

kinetics <- function(tab) {
    time<-tab[[1]]
    n <- length(tab)-1 # number of states

    plot(c(),c(),
         xlim=c(min(time),max(time)),
         ylim=c(0,1),
         xlab="Time",
         ylab="p",
         log="x",
         main=sub(".kin","",basename(input))
         )

    niceidxs <- c()
    nicemaxs <- c()
    
    threshold <- 0.025
    
    idx=0
    for( i in 1:n ) {
        column <- tab[[i+1]]
        mx <- max(column)
        
        if (mx > threshold) {
            idx = idx+1
            niceidxs[idx]=i
            nicemaxs[idx]=round(mx,3)
            lines(time,column,col=idx,lwd=2)
        }
    }
    
    legend("left",legend=paste(niceidxs-1," (",nicemaxs,")"),fill=palette(),inset=0.01,border=F)

    ## return
    niceidxs
}

evaluation_plot <- function(tab,info) {
    time<-tab[[1]]
    n <- length(tab)-1 # number of states
    
    plot(c(),c(),
         xlim=c(min(time),max(time)),
         ylim=c(0,1),
         xlab="time",
         ylab="p",
         log="x",
         main="evaluation")
    
    T0<-time[1]
    T8<-time[length(time)]
    
    timescale=time[2]/time[1] ## scale of one time step
    
    p0<-tab[[3]]
    pi0 = p0[length(p0)]
    p1<-tab[[2]]
    pi1 = p1[length(p1)]
    
    y0 <- 1-(1-pi0)^(-1)*p0
    y1 <- pi1^(-1)*p1
    
    polygon(c(time,rev(time)),c(y0,rev(y1)),col="lightblue")
    lines(time, y1, col=1, lwd=4)
    lines(time, y0, col=2, lwd=3)
    
    lines(time, y0-y1, col="blue", lty=2)
    
    intervalwidth=log(10)/log(timescale)
    score <- sum(y0-y1) / intervalwidth  ## area between the curves (in log_10 time)
    convergencetime <- (log(T8)-log(T0))/log(10)
    nscore <- score/(log(T8)-log(T0))*log(10)
    
    text(T8,0.3,paste(sep="; ",round(score,2),round(convergencetime,2),round(nscore,2)),adj=1);
    
    text(T8,0.4,info,adj=1)
    
    cat(paste(sep="\t",basename(input),score,convergencetime,nscore,info),"\n");
    
    legend("bottomright",
           legend=c("1-p0 (norm'd)","p1 (norm'd)","difference"),
           fill=c(1,2,"blue"),
           inset=0.01)
}

ppdotplot <- function(ppfile,idxs,seqA,seqB,nameA,nameB,restrictdp) {
    lenA=nchar(seqA)
    lenB=nchar(seqB)
    
    ## iterate over lines of ppfile
    ##  + plot dot plot for each line

    lines<-readLines(ppfile)
    lines<-lines[idxs]

    if (restrictdp) {
        minxlim=lenA;
        maxxlim=1;
        minylim=lenB;
        maxylim=1;
        
        for (line in lines) {
            print(line);
            entries <- unlist( strsplit(line,"[ \t]") );
            if (length(entries)<3) {next;}
            entries<-entries[3:length(entries)]
            for ( k in seq(1,length(entries),3)) {
                i = as.numeric(entries[k+0]);
                j = as.numeric(entries[k+1]);
                print(c(i,j));
                minxlim=min(minxlim,i)
                minylim=min(minylim,j)
                maxxlim=max(maxxlim,i)
                maxylim=max(maxylim,j)
            }
        }
    } else {
        minxlim=1; maxxlim=lenA;
        minylim=1; maxylim=lenB;
    }

    print(c(minxlim,maxxlim,minylim,maxylim));

    addxlim=if (maxxlim-minxlim<40) 1 else 0
    addylim=if (maxylim-minylim<40) 1 else 0
    
    plot(c(),c(),main="ipps",xlab=nameA,ylab=nameB,xlim=c(minxlim-addxlim,maxxlim+addxlim),ylim=c(minylim-addylim,maxylim+addylim),axes=F)

    box()
    axis(1,at=c(1,seq(10,lenA,10)),cex.axis=0.6)
    axis(2,at=lenB+1-c(1,seq(10,lenB,10)),labels=c(1,seq(10,lenB,10)),cex.axis=0.6)
    

    
    indices=c()
    color=1
    colors=c()

       
    for (i in seq(5,lenA,5)) {
        lwd <- if (i%%10==0) 0.8 else 0.4
        abline(v=i,lty=2,lwd=lwd,col="grey")
    }
    for (i in seq(5,lenB,5)) {
        lwd <- if (i%%10==0) 0.8 else 0.4
        abline(h=lenB+1-i,lty=2,lwd=lwd,col="grey")
    }

    for (line in lines) {
        entries <- unlist( strsplit(line,"[ \t]") );
        index<-as.integer(entries[1]);
        # pf<-as.numeric(entries[2]);
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
    
    for (i in minxlim:maxxlim) {
        c <- substr(seqA,i,i)
        text(i,minylim-1,c,adj=0.5,cex=0.25)
        text(i,maxylim+1,c,adj=0.5,cex=0.25)
    }

    for (i in minylim:maxylim) {
        c <- substr(seqB,i,i)
        text(minxlim-1,i,c,adj=0.5,cex=0.25)
        text(maxxlim+1,i,c,adj=0.5,cex=0.25)
    }

    #legend("left",legend=indices,fill=colors,border=F,inset=0.01)
    
}

### end functions
############################################################

## setup graphics output
numofplots=1;
if (evalplot) {numofplots=numofplots+1;}
if (pps!="") {numofplots=numofplots+1;}

pdf(output,width=6*numofplots,height=7)

par(mfrow=c(1,numofplots),mar=c(4.25,4.25,2,1))

palette(c(rgb(0,0,0.4,1),
          adjustcolor(rainbow(8),alpha.f=0.6),
          adjustcolor(rainbow(8), alpha.f=0.3)));

## get input table
input_tab<-read.table(input)


########################################
## draw plots

niceidxs <- kinetics(input_tab)

if (evalplot) {
    evaluation_plot(input_tab,info)
}

### optionally draw dot plot (for niceidxs)
if (ppfile != "") {
    ppdotplot(ppfile,niceidxs,seqA,seqB,nameA,nameB,restrictdp)
}

## close graphics output
dev.off();
