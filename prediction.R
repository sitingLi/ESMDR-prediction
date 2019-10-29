"lik1"<-
function(score,time,surv.st)
{
    n<-length(score)  
    r<-rank(time)
    ita<-score
    epita<-exp(score)
    d<-rep(0,n)
    dono<-rep(0,n)
    for(i in 1:n){	
	   d[i]<-sum(surv.st[r==r[i]])
	   dono[i]<-sum(epita[r>=r[i]])
   }

   lik<-sum((ita-log(dono))*surv.st)	
	return(lik)
}



files1 <- dir(pattern="4.txt")
files2 <- dir(pattern="csv")

i=8
data1 <- read.delim(files1[i],header=F)
data2 <- read.csv(files2[i])
y <- data2[,2649]
id1 <- as.vector(unlist(data1))
id1 <- unlist(strsplit(id1, " "))
id1 <- id1[2*(1:(length(id1)/2))]
id1 <- as.numeric(id1)
response <- read.table("response.txt", sep=" ", header=F)
yy <- response[,5]
asd <- match(id1, response[,1])
sum(abs(response[asd,5] -y))

library(glmnet)
library(survival)

cov <- response[asd, ]
fit <- coxph(Surv(cov[,2], cov[,3])~factor(data2[,1])+ cov[,4])
fit <- survfit(Surv(cov[,2], cov[,3])~factor(data2[,1]))
fit <- survdiff(Surv(cov[,2], cov[,3])~factor(data2[,1]))

library(Hmisc)

for(i in 1:dim(data2)[2]){

asd <- data2[,i]
if (sum(is.na(asd))){
randpool <- as.numeric(names(table(asd)))
p <- table(asd)/sum(table(asd))
if (length(p)==1) asd[is.na(asd)] <- rep(randpool, sum(is.na(asd)))
else{ 
replac <- as.numeric(rMultinom(rbind(p,p),sum(is.na(asd)))[1,])
asd[is.na(asd)] <- replac
}
data2[,i] <- asd
}


}



smoke <- rep(0, dim(data2)[1])
smoke[cov[,4]==2] <- 1
x <- as.matrix(cbind(smoke, data2))
snpname <- c("smoke", names(data2)[-dim(data2)[2]])
x <- x[,-dim(x)[2]]
#tempx <- apply(x, 1, is.double)

#cv.fit <- cv.glmnet(x, Surv(cov[,2], cov[,3]), family="cox", lambda= (90:150)/10000, penalty.factor=c(0,rep(1,dim(data2)[2]-1)), maxit = 2000)




i=7
data1 <- read.delim(files1[i],header=F)
data2 <- read.csv(files2[i])
y <- data2[,2649]
id1 <- as.vector(unlist(data1))
id1 <- unlist(strsplit(id1, " "))
id1 <- id1[2*(1:(length(id1)/2))]
id1 <- as.numeric(id1)
response <- read.table("response.txt", sep=" ", header=F)
yy <- response[,5]
asd <- match(id1, response[,1])
sum(abs(response[asd,5] -y))
cov1 <- response[asd, ]
smoke <- rep(0, dim(data2)[1])
smoke[cov1[,4]==2] <- 1
x1 <- as.matrix(cbind(smoke, data2))
x1 <- x1[,-dim(x1)[2]]


fit <- coxph(Surv(cov1[,2], cov1[,3])~factor(data2[,1])+ cov1[,4])
fit <- survfit(Surv(cov1[,2], cov1[,3])~factor(data2[,1]))
fit <- survdiff(Surv(cov1[,2], cov1[,3])~factor(data2[,1]))

pdf('124.pdf')
par(mfrow=c(2,2))
for(i in 1:4){

fit <- survfit(Surv(cov1[,2], cov1[,3])~factor(data2[,i]))
fit1 <- survdiff(Surv(cov1[,2], cov1[,3])~factor(data2[,i]))
plot(fit,main=names(data2)[i])
text(30, 0.2, paste("Test statsitcs is", round(fit1$chisq,1)))

}

dev.off()

i=8
data1tr <- read.delim(files1[i],header=F)
data2tr <- read.csv(files2[i])
id1 <- as.vector(unlist(data1tr))
id1 <- unlist(strsplit(id1, " "))
id1 <- id1[2*(1:(length(id1)/2))]
id1 <- as.numeric(id1)
response <- read.table("response.txt", sep=" ", header=F)
yy <- response[,5]
asd <- match(id1, response[,1])
cov <- response[asd, ]

tt2 <- ssimupred(data2tr[,1:2], data2[1,2], cov[,2:3] , 102)
fit1 <- coxph(Surv(cov1[,2], cov1[,3])~factor(tt2))
summary(fit1)$coef[2]


toplist <- read.delim("snpname-all.txt")
asd1 <- match(toplist[,2], names(data2))

result <- matrix(-10, length(asd1), 4)
for(j in 1:length(asd1)){
i <- asd1[j]
if(!is.na(i)){
if(length(unique(data2[!is.na(data2[,i]),i])) > 1) {fit1 <- survdiff(Surv(cov1[,2], cov1[,3])~factor(data2[,i]))
result[j,1] <- fit1$chisq
result[j,2] <- -pchisq(fit1$chisq,(length(unique(data2[!is.na(data2[,i]),i]))-1), lower.tail=F, log.p=T)
}
if(length(unique(x[!is.na(x[,(i+1)]),(i+1)])) > 1) {fit1 <- survdiff(Surv(cov[,2], cov[,3])~factor(x[,(i+1)]))
result[j,3] <- fit1$chisq
result[j,4] <- -pchisq(fit1$chisq,(length(unique(x[!is.na(x[,(i+1)]),(i+1)]))-1), lower.tail=F, log.p=T)
}

}
}

write(t(result), "temp.txt", ncol=4, sep='\t')



for(i in 1:dim(data2)[2]){

asd <- data2[,i]
if (sum(is.na(asd))){
randpool <- as.numeric(names(table(asd)))
p <- table(asd)/sum(table(asd))
if (length(p)==1) asd[is.na(asd)] <- rep(randpool, sum(is.na(asd)))
else{ 
replac <- as.numeric(rMultinom(rbind(p,p),sum(is.na(asd)))[1,])
asd[is.na(asd)] <- replac
}
data2[,i] <- asd
}


}

res <- matrix(0, dim(x)[2], 1000)

for(i in 1:1000){

fit <- glmnet(x, Surv(cov[,2], cov[,3]), family = "cox", lambda= ((90+i))/10000, penalty.factor=c(0,rep(1,dim(data2)[2])), maxit = 2000)
res[, i] <- as.vector(fit$beta)
cat("i=", i, "\n")

}

bic <-0

for(j in 1:80){

i <- 10*j 
score <- x%*%res[, i]
bic[j] <- -2*(lik1(score, cov[,2], cov[,3])) + log(dim(x)[1])*sum(res[,i]!=0)


}


plot for lars-lym

ins<-c(100,218,250,300,400,500)
lexin.score<-matrix(0,dim(x)[1],6)
for(i in 1:6){
	lexin.score[,i]<-x%*%res[,ins[i]]
}

library(survivalROC)



ins<-c(100,218,250,300,400,500)
lexin.score<-matrix(0,dim(x1)[1],6)
for(i in 1:6){
	lexin.score[,i]<-x1%*%res[,ins[i]]
}



lym2roc2<-matrix(0,9,6)
for (j in 1:6){
for (i in 1:9) {
roc12 <- survivalROC(cov1[,2], cov1[,3],lexin.score[,j], predict.time = (10*i+14), method="KM")
lym2roc2[i,j] <- roc12$AUC
}
}

temproc <- 0
for (i in 1:9) {
roc12 <- survivalROC(cov1[,2], cov1[,3],x1[,1], predict.time = (10*i+14), method="KM")
temproc[i] <- roc12$AUC
}

lym2roc2 <- cbind(lym2roc2, temproc)

timep <- 10*(1:9) + 14

matplot(timep,lym2roc2,lty=1:7,type="l",col=1:7,xlab="year", ylab="Area under the curve")
#par(new=T)
#plot(1:20, lym2roc2[,4], lwd=4, ylim=range(lym2roc2[,c(1,3,4,,6,8,11)]),type="l",xlab="time", ylab="Area under the curve")
legend("bottomright", c("183 gene","29 gene","19 gene","13 gene","4 gene","2 gene", "smoking only"),lty=1:7, col=1:7) 

write(t(lym2roc2), "aucall.txt", ncol=7, sep="\t")
savePlot("aucall.pdf", type="pdf")