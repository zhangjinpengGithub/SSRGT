args = commandArgs(trailingOnly = TRUE)
pure_pchis=args[1]
x <- read.table(pure_pchis,header=FALSE)
x=as.matrix(x)
y=matrix(1:nrow(x),ncol=1)
x=cbind(x,y)
x=cbind(x,y)
for(i in 1:nrow(x)){
    a=as.numeric(x[i,2]) +as.numeric(x[i,3])
    x[i,4]=a/2
    e1=as.numeric(x[i,4])
    x2 <- (as.numeric(x[i,2])-e1)^2/e1 + (as.numeric(x[i,3]) -e1)^2/e1
    x[i,5] <- 1 - pchisq(x2,1)
}
write.table(x,args[2],quote = FALSE,row.names=F,col.names=F,sep="\t")
