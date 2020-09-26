library(edgeR)
library(limma)
library(qvalue)
load("TIL.Rdata")
ls()
head(count.data)
head(info)
Counts=count.data
level = info$TIL 
hpv = info$HPV

names = rownames(Counts)
y= DGEList(counts = Counts, genes = names)
boxplot(split(log10(y$samples$lib.size),hpv),main = "log10 library size")

isexpr = rowSums(cpm(y)>1) >=20
hasannot = rowSums(is.na(y$genes))==0
y = y[isexpr& hasannot, , keep.lib.sizes = FALSE]

gene.mean = apply(y$count,1,mean)
cv.sq = apply(y$count,1,var)/gene.mean^2
plot(log(gene.mean), log(cv.sq), cex=0.3,col="grey", ylim = c(-5,3))
abline(c(0,-1))

y <- calcNormFactors(y)
sf <- colSums(y$count); sf <- sf/mean(sf)
eff.sf <- sf*y$samples$norm.factors
norm.y <- y$count %*% diag(1/eff.sf)
gene.mean <- apply(norm.y,1,mean)
cv.sq <- apply(norm.y,1,var)/gene.mean^2
plot(log(gene.mean), log(cv.sq),cex=0.3,col="grey",ylim=c(-5,3))
abline(c(0,-1))

design = model.matrix(~hpv)
y = estimateDisp(y, design, robust=TRUE)
fit =glmFit(y, design, robust=TRUE)
lrt = glmLRT(fit)
topTags(lrt,n=15)
plotBCV(y,main="BCV plot")

pval = lrt$table$PValue
hist(pval, main = "histogram of HPV p-values")
sprintf("there are %i genes with p-values <=0.05.",sum(pval <=0.05))

fdr.pass <- qvalue(lrt$table$PValue,fdr.level=0.05)
hist(fdr.pass$qvalue,main="FDR, HPV Data")
plot(fdr.pass$lambda,fdr.pass$pi0.lambda,xlim=c(0,1)) 
lines(fdr.pass$lambda,fdr.pass$pi0.smooth) 
abline(h=fdr.pass$pi0,col=2)
title("Storey pi0 estimates for HPV data")
plotBCV(y,main="BCV plot",col.tagwise=as.numeric(fdr.pass$qvalue<=0.05)+1,cex=0.5)

sprintf("there are %i genes with q-values <=0.05.",sum(fdr.pass$qvalue <=0.05))

genenames <- row.names(y)[fdr.pass$significant]
head(genenames)
##############################################################################

mydata = load("TIL.Rdata")
levelmod <- info[info$TIL=="Moderate", ]
head(levelmod)
Counts=count.data%>% select(levelmod$Patient)
View(Counts)
level = levelmod$HPV


y= DGEList(counts = Counts)
boxplot(split(log10(y$samples$lib.size),level),main = "log10 library size")


gene.mean <- apply(y$counts,1,mean)

CV.sq <- apply(y$counts,1,var)/gene.mean^2

plot(log(gene.mean), log(CV.sq),cex=0.3,col="grey",ylim=c(-5,3))
abline(c(0,-1))

y <- calcNormFactors(y)
sf <- colSums(y$count); sf <- sf/mean(sf)
eff.sf <- sf*y$samples$norm.factors
norm.y <- y$count %*% diag(1/eff.sf)
gene.mean <- apply(norm.y,1,mean)
cv.sq <- apply(norm.y,1,var)/gene.mean^2
plot(log(gene.mean), log(cv.sq),cex=0.3,col="grey",ylim=c(-5,3))
abline(c(0,-1))

design = model.matrix(~level)
y = estimateDisp(y, design, robust=TRUE)
fit =glmFit(y, design, robust=TRUE)
lrt = glmLRT(fit)
topTags(lrt,n=15)
plotBCV(y,main="BCV plot")

pval = lrt$table$PValue
hist(pval, main = "histogram of hpv p-values")
sprintf("there are %i genes with p-values <=0.05.",sum(pval <=0.05))
genenames <- row.names(y)[fdr.pass$significant]
fdr.pass <- qvalue(lrt$table$PValue,fdr.level=0.05)
hist(fdr.pass$qvalue,main="FDR, hpv Data")
plot(fdr.pass$lambda,fdr.pass$pi0.lambda,xlim=c(0,1)) 
lines(fdr.pass$lambda,fdr.pass$pi0.smooth) 
abline(h=fdr.pass$pi0,col=2)
title("Storey pi0 estimates for hpv data")
plotBCV(y,main="BCV plot",col.tagwise=as.numeric(fdr.pass$qvalue<=0.05)+1,cex=0.5)

sprintf("there are %i genes with q-values <=0.05.",sum(fdr.pass$qvalue <=0.05))

genenames <- row.names(y)[fdr.pass$significant]
head(genenames)
