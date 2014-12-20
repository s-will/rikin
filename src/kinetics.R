#!/usr/bin/Rscript

### parse command line
args <- commandArgs(trailingOnly = TRUE)
if (!(length(args)==2 || length(args)==3)) {
    print("USAGE: kinetics.R input output [info]");
    quit();
}
input <- args[1]
output <- args[2]
info <- if (length(args>=3)) args[3] else "";
### end parse command line

if ( !file.exists(input) ) {
    inputgz <- paste(sep="",input,".gz");
    if (file.exists(inputgz)) {
        input<-inputgz
    }
}

tab<-read.table(input)

time<-tab[[1]]
n <- length(tab)-1 # number of states

pdf(output,width=12,height=6)

par(mfrow=c(1,2))

plot(c(),c(),
     xlim=c(min(time),max(time)),
     ylim=c(0,1),
     xlab="Time",
     ylab="p",
     log="x",
     main=basename(input))

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

legend("right",legend=paste(niceidxs," (",nicemaxs,")"),fill=palette(),inset=0.01)


plot(c(),c(),
     xlim=c(min(time),max(time)),
     ylim=c(0,1),
     xlab="Time",
     ylab="p",
     log="x",
     main=basename(input))

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

legend("bottomright",legend=c("1-p0 (norm'd)","p1 (norm'd)","difference"),fill=c(1,2,"blue"),inset=0.01)
