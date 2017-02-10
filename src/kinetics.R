#!/usr/bin/env Rscript

library("shape")

## find the directory where this script resides (emulates perl's bindir)
bindir <- function() {
    cA <- commandArgs(trailingOnly=F)
    script.file <- (substring(grep("--file",cA,value=T),8))
    normalizePath(dirname(script.file))
}


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

parser$add_argument("--dotplots", action="store_true", default=FALSE, help="draw dotplots [default: off]")

parser$add_argument("input",nargs=1,help="input file [required]")

parser$add_argument("--int-startA", default="0",
                    help="start of interaction in RNA A")

parser$add_argument("--int-endA", default="0",
                    help="end of interaction in RNA A")

parser$add_argument("--int-startB", default="0",
                    help="start of interaction in RNA B")

parser$add_argument("--int-endB", default="0",
                    help="end of interaction in RNA B")

parser$add_argument("--interaction", default="",
                    help="interaction pattern (dot-bracket, order B&A !)")


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

intStartA <- as.integer(args$int_startA);
intEndA <- as.integer(args$int_endA);
intStartB <- as.integer(args$int_startB);
intEndB <- as.integer(args$int_endB);

interaction <- args$interaction;

dotplots <- args$dotplots

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

########################################
## source kinetics plotting and evaluation functions
#
source((file.path(bindir(),"kinetics_lib.R")))
########################################

## setup graphics output

numofplots=1;

if (evalplot) {numofplots=numofplots+1;}
if (pps!="") {numofplots=numofplots+1;}

pdf(output,width=7*numofplots,height=7)

par(mfrow=c(1,numofplots),mar=c(4.25,4.25,2,2))


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
    sink(paste(sep="",output,".txt"))
    time_interaction_probability_plot(input_tab,ppfile,niceidxs,seqA,nameA,intStartA,intEndA)
    sink()
}

### optionally draw dot plots (for niceidxs)
if (dotplots && ppfile != "") {
    par(mfrow=c(1,2))

    ppdotplot(ppfile,niceidxs,niceweights,seqA,seqB,nameA,nameB,intStartA,intEndA,intStartB,intEndB)
    interaction_probability_plot(ppfile,niceidxs,niceweights,seqA,nameA)
}
