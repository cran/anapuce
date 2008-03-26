`LBI` <-
function(infile,name.M="M.norm",ind.array=1:2,graph=TRUE,graphout="FigM1M2") {

ind <- which(names(infile) %in% paste(name.M,ind.array,sep=""))
cond <- as.data.frame(apply(apply(infile[,ind],1,is.na),2,sum))
nbgene1 <- table(cond)[unique(cond)==1]
condb <- as.data.frame(apply(apply(infile[,ind],1,is.na),2,sum))
nbgene2 <- table(condb)[unique(condb)==2]
cat("Number of genes with only one observation :",nbgene1,"\n")
cat("Number of genes wihtout any observation :",max(nbgene2,0),"\n")

VarParGene2 <- apply(infile[,ind],1,var,na.rm=TRUE)
VarParGene1 <- apply(cbind(infile[,ind[1]],-infile[,ind[2]]),1,var,na.rm=TRUE)

Var2 <- mean(VarParGene2,na.rm=TRUE)
Var1 <- mean(VarParGene1,na.rm=TRUE)
cat("Variance per gene (Mean Log-ratio) : ", Var1,"\n")
cat("Variance per gene (Dye bias) : ", Var2,"\n")
cat("LBI : ", Var1/Var2,"\n")

if (graph) {
postscript(file=paste(graphout,".ps",sep=""))  
plot(infile[,ind[1]],infile[,ind[2]],xlim=c(-3,3),ylim=c(-3,3),xlab="log2(R1/G1)",ylab="log2(R2/G2)")
dev.off()
}
# (c) 2007 Institut National de la Recherche Agronomique
 
}

