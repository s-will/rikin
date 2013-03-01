t<-read.table("kinetics")

time<-t[[1]]
n <- length(t)-1 # number of states

colors = rainbow(5)

pdf("kinetics.pdf")

plot(c(),c(),xlim=c(1,max(time)),ylim=c(0,1),xlab="Time",ylab="p",log="x")

niceidxs <- c()
nicemaxs <- c()

threshold <- 0.025

idx=0
for( i in 1:n ) {
  column <- t[[i+1]]
  mx <- max(column)
  
  if (mx > threshold) {
    idx = idx+1
    niceidxs[idx]=i
    nicemaxs[idx]=round(mx,3)
    lines(time,t[[i+1]],col=idx,lwd=2)
  }
}

legend("right",legend=paste(niceidxs," (",nicemaxs,")"),fill=palette(),inset=0.01)
