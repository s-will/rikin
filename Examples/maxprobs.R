t<-read.table("kinetics")

time<-t[[1]]

n <- length(t)-1 # number of states

pdf("maxprobs.pdf")

maxima <- c()

threshold <- 0.02
nicemaxima <- c()
niceidxs <- c()

idx=0
for( i in 1:n ) {
  column <- t[[i+1]]
  mx <- max(column)
  maxima[i]=c(mx)

  if (mx > threshold) {
    idx = idx+1
    nicemaxima[idx]=mx
    niceidxs[idx]=i
  }
}

print(niceidxs)
print(nicemaxima)

plot(nicemaxima,axes=F)

axis(1,labels=niceidxs,at=1:idx)
axis(2)
box()
