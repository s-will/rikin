#!/usr/bin/Rscript

library("shape")

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
parser$add_argument("--restrict-dp", action="store_true", default=FALSE,
                    help="restrict dot plot to region of probabilities above threshold [default: off]")

parser$add_argument("--relative-dp", action="store_true", default=FALSE,
                    help="additionally weight dots by maximum probability of macrostate [default: off]")

parser$add_argument("--dproot", default="2",
                    help="transformation of probability to dot-size [default: 2 (circle area); e.g. use 3 for sphere volume]")

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
relativedp <- args$relative_dp
dproot <- as.numeric(args$dproot)

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

########################################
## define color palette

library("RColorBrewer")

mypalette <- c(rgb(0,0,0.4,1),
               rgb(0.8,0,0,1),
               adjustcolor(brewer.pal(8,"Dark2"),alpha.f=0.9),
               adjustcolor(brewer.pal(8,"Set2"),alpha.f=0.8),
               adjustcolor(brewer.pal(8,"Pastel2"),alpha.f=0.7))

## has to be the same as mypalette but without transparency
myopaquepalette <- c(rgb(0,0,0.4,1),
                     rgb(0.8,0,0,1),
                     brewer.pal(8,"Dark2"),
                     brewer.pal(8,"Set2"),
                     brewer.pal(8,"Pastel2"))

plotBg=grey(0.975)

############################################################
### functions


clearPlotBg <- function(color,log) {
    lim=par("usr")
    xmin=lim[1]
    xmax=lim[2]
    ymin=lim[3]
    ymax=lim[4]
    if (log=="x" || log=="xy") {
        xmin=10^xmin
        xmax=10^xmax
    }
    if (log=="y" || log=="xy") {
        ymin=10^ymin
        ymax=10^ymax
    }

    ## print (c(xmin,xmax,ymin,ymax))
    rect(xmin,ymin,xmax,ymax,col = color)
}

kinetics_plot <- function(tab,title) {
    time<-tab[[1]]
    n <- length(tab)-1 # number of states

    plot(c(),c(),
         xlim=c(min(time),max(time)),
         ylim=c(0,1),
         xlab="Time",
         ylab="State Probability",
         log="x",
         main=title
         )

    clearPlotBg(plotBg,"x")
    
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
            lines(time,column,col=mypalette[idx],lwd=2)
        }
    }
    
    legend("left",
           legend=paste(niceidxs-1," (",nicemaxs,")"),
           fill=mypalette,inset=0.02,
           border=grey(0.7),
           ## cex=0.85,
           bg=grey(0.9,0.8)
           )

    ## return
    niceidxs
}

evaluation_plot <- function(tab,info,instanceName) {
    time<-tab[[1]]
    n <- length(tab)-1 # number of states
    
    plot(c(),c(),
         xlim=c(min(time),max(time)),
         ylim=c(0,1),
         xlab="Time",
         ylab="State Probability",
         log="x",
         main="Evaluation")
    clearPlotBg(plotBg,"x")
    
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
    lines(time, y1, col=mypalette[1], lwd=4)
    lines(time, y0, col=mypalette[2], lwd=3)
    
    lines(time, y0-y1, col="blue", lty=2)
    
    intervalwidth=log(10)/log(timescale)
    score <- sum(y0-y1) / intervalwidth  ## area between the curves (in log_10 time)
    convergencetime <- (log(T8)-log(T0))/log(10)
    nscore <- score/(log(T8)-log(T0))*log(10)
    
    text(T8,0.3,paste(sep="; ",round(score,2),round(convergencetime,2),round(nscore,2)),adj=1);
    
    text(T8,0.4,info,adj=1)
    
    cat(paste(sep="\t","#evaluation",instanceName,score,convergencetime,nscore,info),"\n");
    
    legend("bottomright",
           legend=c("1-p0 (norm'd)","p1 (norm'd)","difference"),
           fill=c(1,2,"blue"),
           inset=0.01)
}

ppdotplot <- function(ppfile,idxs,weights,seqA,seqB,nameA,nameB) {
    lenA=nchar(seqA)
    lenB=nchar(seqB)
    
    ## iterate over lines of ppfile
    ##  + plot dot plot for each line

    lines<-readLines(ppfile)
    ## lines<-lines[idxs]

    if (restrictdp) {
        minxlim=lenA;
        maxxlim=1;
        minylim=lenB;
        maxylim=1;
        
        for (idx in idxs) {
            line <- lines[idx]
            ## print(line);
            entries <- unlist( strsplit(line,"[ \t]") );
            if (length(entries)<3) {next;}
            entries<-entries[3:length(entries)]
            for ( k in seq(1,length(entries),3)) {
                i = as.numeric(entries[k+0]);
                j = as.numeric(entries[k+1]);
                ## print(c(i,j));
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

    ## print(c(minxlim,maxxlim,minylim,maxylim));

    addxlim=if (maxxlim-minxlim<40) 1 else 0
    addylim=if (maxylim-minylim<40) 1 else 0
    
    plot(c(),c(),
         main="Interaction pair probabilities",
         xlab=nameA,
         ylab=nameB,
         xlim=c(minxlim-addxlim,maxxlim+addxlim),
         ylim=c(minylim-addylim,maxylim+addylim),
         axes=F)
    clearPlotBg(plotBg,"")

    box()
    axis(1,at=c(1,seq(10,lenA,10)),cex.axis=0.6)
    axis(2,at=lenB+1-c(1,seq(10,lenB,10)),labels=c(1,seq(10,lenB,10)),cex.axis=0.6)
    
    ## indices=c()
    ##idxOfIdx=1
    ## colors=c()

       
    for (i in seq(5,lenA,5)) {
        lwd <- if (i%%10==0) 0.8 else 0.4
        abline(v=i,lty=2,lwd=lwd,col="grey")
    }
    for (i in seq(5,lenB,5)) {
        lwd <- if (i%%10==0) 0.8 else 0.4
        abline(h=lenB+1-i,lty=2,lwd=lwd,col="grey")
    }

    dimx <- maxxlim-minxlim+1
    dimy <- maxylim-minylim+1
    matrixP <- as.list(rep(0,dimx*dimy))
    dim(matrixP)=c(dimx,dimy)
    matrixCol <- as.list(rep(0,dimx*dimy))
    dim(matrixCol)=c(dimx,dimy)
    
    for (idxOfIdx in 1:length(idxs)) {
        idx=idxs[idxOfIdx]
        
        line <- lines[idx]
        
        entries <- unlist( strsplit(line,"[ \t]") );
        ## index<-as.integer(entries[1]);

        ## pf<-as.numeric(entries[2]);
        if (length(entries)>2) {

            entries<-entries[3:length(entries)]
            ## indices<-c(indices, index)
            
            for ( k in seq(1,length(entries),3)) {
                i = as.numeric(entries[k+0]);
                j = as.numeric(entries[k+1]);
                p = as.numeric(entries[k+2]) * weights[idxOfIdx];

                ##registerDot(i,j,p,idxOfIdx);
                ii=i-minxlim+1
                jj=j-minylim+1
                
                matrixP[[ii,jj]] <- c(matrixP[[ii,jj]],p)
                matrixCol[[ii,jj]] <- c(matrixCol[[ii,jj]],idxOfIdx) 
                
                ## p=sqrt(p)
                ## points(i,j,cex=p,col=idxOfIdx,pch=16);
                ## print(c(i,j,p))
            }
        }
        ## colors=c(colors,idxOfIdx)
    }

    for (i in minxlim:maxxlim) {
        for (j in minylim:maxylim) {
            ii=i-minxlim+1
            jj=j-minylim+1
            ps <- matrixP[[ii,jj]]
            cols <- matrixCol[[ii,jj]]

            ps <- tail(ps,-1)
            cols <- tail(cols,-1)
            
            if (length(ps)>=1) {
                acc <- c(0)
                for (k in 1:length(ps)) {
                    acc=c(acc,acc[k]+ps[k])
                }
                total=acc[length(ps)+1]
                
                for (k in 1:length(ps)) {
                    rx = total^(1/dproot)*sqrt(2);

                    e <-
                        getellipse(
                            mid=c(i,j),
                            rx=rx,
                            ry=rx*dimy/dimx,
                            dr=0.1,
                            from=acc[k]/total*2*pi,
                            to=acc[k+1]/total*2*pi);

                    e <- rbind(c(i,j),e,c(i,j)) ## connect to center! ==> sector

                    polygon(e,col=mypalette[cols[k]],lty=0)
                }
            }
        }
    }
    
    ## ------------------------------------------------------------
    ## write sequence A (x-axis)
    for (i in minxlim:maxxlim) {
        c <- substr(seqA,i,i)
        text(i,minylim-1,c,adj=0.5,cex=0.25)
        text(i,maxylim+1,c,adj=0.5,cex=0.25)
    }

    ## ------------------------------------------------------------
    ## write sequence B (y-axis)
    for (i in minylim:maxylim) {
        c <- substr(seqB,i,i)
        text(minxlim-1,i,c,adj=0.5,cex=0.25)
        text(maxxlim+1,i,c,adj=0.5,cex=0.25)
    }

    ## legend("left",legend=indices,fill=colors,border=F,inset=0.01)
    
}

## plot interaction probability of each position of RNA A vs. each state
interaction_probability_plot <- function(ppfile,idxs,weights,seqA,nameA) {
    lenA=nchar(seqA)
    lenIdxs=length(niceidxs)
    
    lines<-readLines(ppfile)
    ## lines<-lines[idxs]

    minxlim=1; maxxlim=lenA;
    minylim=1; maxylim=lenIdxs;
    
    ## print(c(minxlim,maxxlim,minylim,maxylim));

    addxlim=if (maxxlim-minxlim<40) 1 else 0
    addylim=if (maxylim-minylim<40) 1 else 0
    
    plot(c(),c(),
         main="Conditional interaction probabilities",
         xlab=nameA,
         ylab="State",
         xlim=c(minxlim-addxlim,maxxlim+addxlim),
         ylim=c(minylim-addylim,maxylim+addylim),
         axes=F)
    clearPlotBg(plotBg,"")

    box()
    axis(1,at=c(1,seq(10,lenA,10)),cex.axis=0.6)
    axis(2,at=seq(1,lenIdxs,1),labels=rev(idxs-1),cex.axis=0.6,las=2)
    
       
    for (i in seq(0,lenA,5)) {
        lwd <- if (i%%10==0) 0.8 else 0.4
        abline(v=i,lty=2,lwd=lwd,col="grey")
    }

    for (i in seq(0,lenIdxs+1,1)) {
        lwd <- 0.4
        abline(h=lenIdxs+1-i-0.5,lty=2,lwd=lwd,col="grey")
    }

    dimx <- maxxlim-minxlim+1

    for (idxOfIdx in 1:length(idxs)) {
        idx <- idxs[idxOfIdx]
        line <- lines[idx]
        accp <- rep(0,dimx)
             
        entries <- unlist( strsplit(line,"[ \t]") )
        ## index<-as.integer(entries[1]);

        ## pf<-as.numeric(entries[2]);
        if (length(entries)>2) {

            entries <- entries[3:length(entries)]
            ## indices<-c(indices, index)
            
            for ( k in seq(1,length(entries),3)) {
                i = as.numeric(entries[k+0])
                j = as.numeric(entries[k+1])
                p = as.numeric(entries[k+2]) * weights[idxOfIdx]

                ## compute 1-based index
                ii=i-minxlim+1
                
                accp[ii] <- accp[ii] + p
            }

            for (i in minxlim:maxxlim) {
                ii=i-minxlim+1
                p=accp[ii]
                r=p^(1/dproot)
                e <- getellipse(
                    mid=c(i,length(idxs)-idxOfIdx+1),
                    rx=r*1.1, ## somewhat more overlap for high probs
                    ry=r/2,
                    dr=0.1
                    )
                polygon(e,
                        lty=0,
                        col=myopaquepalette[idxOfIdx])
            }
        }
    }

    
    ## ------------------------------------------------------------
    ## write sequence A (x-axis)
    for (i in minxlim:maxxlim) {
        c <- substr(seqA,i,i)
        text(i,minylim-1,c,adj=0.5,cex=0.25)
        text(i,maxylim+1,c,adj=0.5,cex=0.25)
    }

}


## plot interaction probability of each position of RNA A vs. time
time_interaction_probability_plot <- function(kintab,ppfile,idxs,seqA,nameA) {
    lenA=nchar(seqA)
    lenIdxs=length(niceidxs)
    
    time<-kintab[[1]] ## vector of time points in kinetics table
    statesnum <- length(kintab)-1 # number of states in kinetics table

    lines<-readLines(ppfile)

    minxlim=time[1]; maxxlim=time[length(time)];
    minylim=1; maxylim=lenA; ## put sequence A on y-axis!
    
    addxlim=0
    addylim=if (maxylim-minylim<40) 1 else 0
    
    plot(c(),c(),
         main="Interaction probabilities vs. time",
         xlab="Time",
         ylab=nameA,
         xlim=c(minxlim-addxlim,maxxlim+addxlim),
         ## ylim=c(minylim-addylim,maxylim+addylim),
         ylim=c(maxylim+addylim,minylim-addylim),
         axes=F,
         log="x")
    clearPlotBg(plotBg,"x")

    box()
    axis(1)
    axis(2,at=c(1,seq(10,lenA,10)),cex.axis=0.6,las=2)  

    ########################################
    ## collect conditional interaction probabilities for each state in matrix ips

    ips = matrix(0,nrow=length(idxs),ncol=lenA)
    
    for (idxOfIdx in 1:length(idxs)) {
        
        idx=idxs[idxOfIdx]
        line <- lines[idx]

        accp <- rep(0,lenA)
             
        entries <- unlist( strsplit(line,"[ \t]") )

        if (length(entries)>2) {
            
            entries <- entries[3:length(entries)]
            ## indices<-c(indices, index)
            
            for ( k in seq(1,length(entries),3)) {
                i = as.numeric(entries[k+0])
                j = as.numeric(entries[k+1])
                p = as.numeric(entries[k+2])

                accp[i] <- accp[i] + p
            }
            
            for (i in 1:lenA) {
                p=accp[i]
                ips[idxOfIdx,i]=p
            }
        }
    }

    opal = c(rgb(1,1,1),brewer.pal(9,"YlOrRd"))
    pal = opal # adjustcolor(opal,alpha.f=0.2)
    
    ## ----------------------------------------
    ## for each time point compute mixture distribution
    ## and plot
    for (timeIdx in 1:length(time)) {
        theTime=time[timeIdx]
        
        for (i in 1:lenA) {
            accp <- 0
            for(idxOfIdx in 1:length(idxs)) {
                idx=idxs[idxOfIdx]
                accp <- accp + ips[idxOfIdx,i] * kintab[timeIdx,idx+1]
            }
            p=accp^(1/dproot)
            mycol=as.integer(p*9+1)
            if (mycol>1) {
                #print(c(theTime,i,p))
                points(theTime,i,cex=1,pch=15,col=pal[mycol])
            }
        }
    }

        
    ## for (i in seq(0,lenA,5)) {
    ##     lwd <- if (i%%10==0) 0.8 else 0.4
    ##     abline(v=i,lty=2,lwd=lwd,col="grey")
    ## }

    for (i in seq(0,lenA,5)) {
        lwd <- if (i%%10==0) 0.8 else 0.4
        abline(h=i,lty=2,lwd=lwd,col=grey(0.5,0.5))
    }

    
    ## ------------------------------------------------------------
    ## write sequence A (x-axis)
    for (i in minylim:maxylim) {
        c <- substr(seqA,i,i)
        text(minxlim/1.4,i,c,adj=0.5,cex=0.25)
        text(maxxlim*1.4,i,c,adj=0.5,cex=0.25)
    }

    legend("left",
           legend=c("1.0",rep("",8),"0.0"),
           fill=rev(opal),
           inset=0.05,
           border=grey(0.7),
           ## cex=0.85,
           bg=grey(0.9,0.8)
           )
    
}



### end functions
############################################################

## setup graphics output
numofplots=1;
if (evalplot) {numofplots=numofplots+1;}
if (pps!="") {numofplots=numofplots+1;}

pdf(output,width=7*numofplots,height=7)

par(mfrow=c(1,numofplots),mar=c(4.25,4.25,2,1))


## get input table
input_tab<-read.table(input)

title <- "Kinetics"

instanceName<-""

if (nameA!="") {
    instanceName=paste(instanceName,nameA)
}
if (nameB!="") {
    if (nameA!="") {
        instanceName=paste(instanceName,"-",sep="")
    }
    instanceName=paste(instanceName,nameB,sep="")
}

title=paste(title,instanceName)


########################################
## draw plots

niceidxs <- kinetics_plot(input_tab,title)

if (relativedp) {
    niceweights <- c()
    for (i in 1:length(niceidxs)) {
        niceweights[i]=max(input_tab[[niceidxs[i]+1]])
    }
} else {
    niceweights <- rep(1.0,length(niceidxs))
}

if (evalplot) {
    evaluation_plot(input_tab,info,instanceName)
}

if (ppfile != "") {
    time_interaction_probability_plot(input_tab,ppfile,niceidxs,seqA,nameA)
}

### optionally draw dot plots (for niceidxs)
if (ppfile != "") {
    par(mfrow=c(1,2))

    ppdotplot(ppfile,niceidxs,niceweights,seqA,seqB,nameA,nameB)
    interaction_probability_plot(ppfile,niceidxs,niceweights,seqA,nameA)
}

## close graphics output
dev.off();

