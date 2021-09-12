args = commandArgs(trailingOnly = TRUE)
hybrid_pchis=args[1]
x <- read.table(hybrid_pchis,header=FALSE)
x=as.matrix(x)
y=matrix(1:nrow(x),ncol=1)
x=cbind(x,y)
x=cbind(x,y)
for(i in 1:nrow(x)){
    a=as.numeric(x[i,2]) +as.numeric(x[i,3])+as.numeric(x[i,4])+as.numeric(x[i,5])
    x[i,6]=a/4
    e1=as.numeric(x[i,6])
    x2 <- (as.numeric(x[i,2]) -e1)^2/e1+(as.numeric(x[i,3])-e1)^2/e1+(as.numeric(x[i,4])-e1)^2/e1+(as.numeric(x[i,5])-e1)^2/e1
    x[i,7] <- 1 - pchisq(x2,3)
}
write.table(x,args[2],quote = FALSE,row.names=F,col.names=F,sep="\t")

