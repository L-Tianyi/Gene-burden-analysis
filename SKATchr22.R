library(data.table)
chr22 <- as.data.frame(fread("~/NAS/joint_genotyping_anno_chr22.csv"))


###investigate the data set##
cadd <- ggplot(chr22, aes(x = cadd_phred, y = AF)) + geom_point(color = "blue") + labs(x="CADD Score", y = "Allele frequency")
cadd
dann <- ggplot(chr22, aes(x = dann, y = AF)) + geom_point(color = "red") + labs(x="DANN Score", y = "Allele frequency")
dann
categories <- unique(chr22$most_severe_consequence)
length(categories)
table(Genotype$most_severe_consequence)
genotype_rm <- Genotype[chr22$most_severe_consequence != "?",c(21,33)]
box <- ggplot(genotype_rm, aes(x=most_severe_consequence, y = cadd_phred, color = most_severe_consequence)) + geom_boxplot()+ labs(x = "consequences", y = "CADD Score")


###subset data###
chr22 <- chr22[,c(7,29,34:ncol(chr22))] ##allele frequencies, gene symbols, and genotypes are kept##
chrM <- as.matrix(chr22) ##convert into matrix##


###data filtering###
##calculate missing rate for each sample##
M <- array(0, ncol(chrM)-2)
for (i in 1:ncol(chrM)-2) {
M[i]  <- length(which(is.na(chrM[,i+2])))
}
M <- M/nrow(chrM)
chrM <- chrM[,c(1,2,which(M < 0.2)+2)] ##samples with missing > 20% is removed, 1001 samples left
##remove variants that has a high missing rate#
Mv <- array(0, nrow(chrM))
for (i in 1:nrow(chrM)) {
Mv[i]  <- length(which(is.na(chrM[i,3:1003])))
}
Mv <- Mv/(ncol(chrM)-2)
chrM <- chrM[which(Mv<0.2),] ##72554 variant left##
dim(chrM)  ##72554x1003##


###PCA analysis###
#PCA for Population Structure#
missing <- rep(0,nrow(chrM))
for (i in 1:length(missing)) {
  missing[i] <- length(which(is.na(chrM[i,])))
}


