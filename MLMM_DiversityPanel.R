#MLMM for analysing Data
#This script was made using the MLMM package by Segura et al. 2012 who owns the copyright
#I only wrote this code to assist people who need to run MLMM on diversity panel without having to got hrough the trouble of adjusting the original code by themselves
#All gratitude to Dr. Segura et al for their immense contribution to GWAS analysis
 
source("https://raw.githubusercontent.com/Gregor-Mendel-Institute/MultLocMixMod/master/R/mlmm.r")
source("https://raw.githubusercontent.com/Gregor-Mendel-Institute/MultLocMixMod/master/R/mlmm_cof.r")
source("https://raw.githubusercontent.com/Gregor-Mendel-Institute/MultLocMixMod/master/R/plot_mlmm.r")

install.packages("https://github.com/Gregor-Mendel-Institute/mlmm/files/1356516/emma_1.1.2.tar.gz", repos = NULL)
library(emma) 

library(devtools) 

library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) 
library("scatterplot3d")

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

#Load your hapmap data and use GAPT to convert it to numeric format
data_hmp <- read.delim("data.hmp.txt", head=F)
myG <- data_hmp
myGAPIT <- GAPIT(G=myG, output.numerical=TRUE)
dat=myGAPIT$GD

nan=dat[,1]
nan <- as.vector(as.matrix(nan))
dat <- dat[,2:ncol(dat)]
dat=as.matrix(dat)
dat=dat-1
rownames(dat)=nan
##
nam <- data.frame(nan)
names(nam) <- "Taxa"



## Creating Genetic Map
da<-data_hmp[c(1,3:4)]
SNP_INFO<-da[-1,]
colnames(SNP_INFO)<-c("SNP","Chr","Pos")
C<-as.character(SNP_INFO$Chr) 
Ch<-as.integer(C)
Pi<-as.character(SNP_INFO$Pos)
Po<-as.integer(Pi)
SNP_info<-cbind(SNP_INFO,Ch,Po)
snp_info<-SNP_info[c(1,4:5)] 
colnames(snp_info)<-c("SNP","Chr","Pos")

## Load Your Phenotypic Data
Pheno <- read.table("Phenotype.txt", sep="\t", header=T, stringsAsFactor=F) 

ph <- Pheno

## Extract your Kinship data from GAPIT output
kins <- myGAPIT$kinship
#kins <- as.matrix(kins)
k <- kins[,-1]
rownames(k) <- kins[,1]
namkin <- as.vector(as.matrix(kins[,1]))
names(k) <- namkin

# Load your PCA data, you can estimate PCA of your genetic data in TASSEL 5.0 GUI and import it here
myCV <- read.table("PCA.txt", sep="", header=T, stringsAsFactor=F)
colnames(myCV)[1] <- "Taxa"
# Remove first column bearing Genotype IDs
cof_fam2 <- myCV[,-1]
cof_fam2 <- as.matrix(cof_fam2)


rownames(cof_fam2) <- as.vector(as.matrix(myCV$Taxa))

setwd("/homes/omo/2014_NAM_Genomic_Data/BREAD/MLMM0.03")

for(i in 2:ncol(ph)){

  pheno2 <- ph[,c(1,i)]
trait.name <- names(pheno2[2])

print(paste("-------------- Performing GWAS on trait ", trait.name, "!!!!!!!---------------", sep = "")) 
  Y<- as.vector(pheno2[,2]) 
  names(Y) <- pheno2$Taxa
  Y=Y[!is.na(Y)]
  genot.trait=dat[match(names(Y),rownames(dat)),]
  
  k2 <-k[match(names(Y), rownames(k)),]
  k22 <- as.matrix(k2)
  k3 <-k22[,match(names(Y), colnames(k22))]
  #K <- as.matrix(k3)
  
  cov_fam2=cof_fam2[match(names(Y),rownames(cof_fam2)),]
  class(cov_fam2)


mygwas_trait<-mlmm_cof(Y=Y,X=genot.trait, cofs=cov_fam2, K=k3,2,10)

res.Trait=mygwas_trait$opt_mbonf$out
write.table(res.Trait,paste("GWAS_mlmm_results", colnames(Pheno)[i],".csv",sep=""), sep=",", quote=F, row.names=F, col.names=T)
step_table=mygwas_trait$step_table
write.table(step_table,paste("GWAS_mlmm_step_table", colnames(Pheno)[i],".csv",sep=""), sep=",", quote=F, row.names=F, col.names=T)
pval_step=mygwas_trait$pval_step
write.table(step_table,paste("GWAS_mlmm_pval_step", colnames(Pheno)[i],".csv",sep=""), sep=",", quote=F, row.names=F, col.names=T)
  
pdf(paste("step_table_extBIC_",colnames(Pheno)[i],".pdf",sep=""),width=15,height=4,paper='special')
plot_step_table(mygwas_trait,'extBIC')
dev.off()
pdf(paste("step_table_maxpval_",colnames(Pheno)[i],".pdf",sep=""),width=15,height=4,paper='special')
plot_step_table(mygwas_trait,'maxpval')
dev.off()
pdf(paste("step_RSS_",colnames(Pheno)[i],".pdf",sep=""),width=15,height=4,paper='special')
plot_step_RSS(mygwas_trait)
dev.off()
pdf(paste("fwd_GWAS1_",colnames(Pheno)[i],".pdf",sep=""),width=15,height=4,paper='special')
plot_fwd_GWAS(mygwas_trait,1,snp_info,0.1)
dev.off()

pdf(paste("fwd_GWAS2_",colnames(Pheno)[i],".pdf",sep=""),width=15,height=4,paper='special')
plot_fwd_GWAS(mygwas_trait,2,snp_info,0.1)
dev.off()
pdf(paste("fwd_GWAS3_",colnames(Pheno)[i],".pdf",sep=""),width=15,height=4,paper='special')
plot_fwd_GWAS(mygwas_trait,3,snp_info,0.1)
dev.off()
pdf(paste("opt_GWAS_extBIC_",colnames(Pheno)[i],".pdf",sep=""),width=15,height=4,paper='special')
plot_opt_GWAS(mygwas_trait,'extBIC',snp_info,0.1)
dev.off()

pdf(paste("opt_GWAS_mbonf_",colnames(Pheno)[i],".pdf",sep=""),width=15,height=4,paper='special')
plot_opt_GWAS(mygwas_trait,'mbonf',snp_info,0.1)
dev.off()
pdf(paste("qqplot_fwd_GWAS_",colnames(Pheno)[i],".pdf",sep=""))
qqplot_fwd_GWAS(mygwas_trait,5)
dev.off()
pdf(paste("qqplot_opt_GWAS_extBIC_",colnames(Pheno)[i],".pdf",sep=""))
qqplot_opt_GWAS(mygwas_trait,'extBIC')
dev.off()
pdf(paste("qqplot_opt_GWAS_mbonf_",colnames(Pheno)[i],".pdf",sep=""))
qqplot_opt_GWAS(mygwas_trait,'mbonf')
dev.off()

}

