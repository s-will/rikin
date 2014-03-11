t<-read.table("kinetics")

openstate<-2

time<-t[[1]]

helix<-t[[2]]
open<-t[[openstate+1]]

pdf("kinetics.pdf")

plot(time,open,t="l",col="red",ylim=c(0,1),xlab="Time",ylab="p",log="x")
lines(time,helix,col="green")

lines(time,t[[4+1]],col="blue")
lines(time,t[[14+1]],col="magenta")
lines(time,t[[27+1]],col="darkgreen")
lines(time,t[[23+1]],col="orange")
lines(time,t[[25+1]],col="cyan")
