###################
####import data####
###################
geno <-as.data.frame(fread("/media/pontikos_nas2/AnitaSzabo/projects/exome_pipeline/data/chr_of_UCLex/unique_databs/new_ctrl_list_FECD_unique_w_REVEL_Kaviar.csv"))

case <- read.table("/media/pontikos_nas2/PCA/Fuchs/varying_filters/postQCcases.tab")
ctrl <- read.table("/media/pontikos_nas2/PCA/Fuchs/varying_filters/postQCcontrols.tab")

geno <- geno[,c(4,9,10,20,29,41:(ncol(geno)-4))] #4=CHROM, 9=AF, 10=AC, 20=gene_symbol 29=cadd_phred


###############
###PCA of PS###
###############

#calculate missing rate for each variant#
M <- array(0, nrow(geno)) #number of samples

subgeno <- geno[,6:ncol(geno)] #2347180 X 1246
subgeno <- as.matrix(subgeno)

for (i in 1:nrow(subgeno)) {
M[i]  <- length(which(is.na(subgeno[i,])))
}

PCgeno <- geno[which(M == 0),] #left the varaiant without missingness for PCA
dim(PCgeno) #1143943 x 1251

#PCA#
matri <- PCgeno[,6:ncol(PCgeno)] #remove annotation columns
matri <- as.matrix(matri) #1143943 x 1246
matri <- t(matri) #1246 x 1143943
class(matri[,1])
matri <- as.numeric(matri)
matri <- matrix(matri,nrow = 1246, ncol = 1143943, byrow = FALSE)
matri <- matri[ , which(apply(matri, 2, var) != 0)] #1246 X 359770 
pcpop <- prcomp(matri, scale. = TRUE)

#make plot#
library(ggplot2)
library(factoextra)
library(plotly)
library(stats)

#screeplot#
pdf(file = "~/fast-storage/FECD_screepop.pdf") 
fviz_eig(pcpop)
dev.off()

#PCA plot#
loadingpop = pcpop$x
loadingpop <- as.data.frame(loadingpop)

fig <- plot_ly(loadingpop, x = loadingpop$PC1, y = loadingpop$PC2, type = 'scatter', mode = 'markers')%>%
  layout(
    plot_bgcolor='#e5ecf6',
    xaxis = list(
      title = "PC1",
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    yaxis = list(
      title = "PC2",
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'))
fig


#########################
###PCA for missingness###
#########################

cutgeno <- geno[,6:ncol(geno)]
cutgeno <- as.matrix(cutgeno)
cutgeno <- t(cutgeno) 
dim(cutgeno) #1246 X 2347180
cutgeno <- as.numeric(cutgeno)
cutgeno <- matrix(cutgeno, nrow = 1246, ncol = 2347180)

#transform the genotype matrix into a missingness matrix
for (i in 1:ncol(cutgeno)) {
  cutgeno[which(cutgeno[,i] == 0),i] <- 1
}
for (i in 1:ncol(cutgeno)) {
  cutgeno[which(cutgeno[,i] == 2),i] <- 1
}
for (i in 1:ncol(cutgeno)) {
  cutgeno[which(is.na(cutgeno[,i])),i] <- 0
}
table(cutgeno, useNA = "always")

#PCA
cutgeno <- cutgeno[ , which(apply(cutgeno, 2, var) != 0)]
pctech <- prcomp(cutgeno, scale. = TRUE)

#plot
loadingtech = pctech$x
loadingtech <- as.data.frame(loadingtech)

fig <- plot_ly(loadingtech, x = loadingtech$PC1, y = loadingtech$PC2, type = 'scatter', mode = 'markers')%>%
  layout(
    plot_bgcolor='#e5ecf6',
    xaxis = list(
      title = "PC1",
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'),
    yaxis = list(
      title = "PC2",
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor='#ffff'))
fig


##################
###creating GRM###
##################

library(AGHmatrix)
grm <- Gmatrix(matri, method = "VanRaden")


##########################
###preparation for SKAT###
##########################

##create vector of phenotype##
y <- c(rep(1,nrow(case)), rep(0,nrow(ctrl)))

#subset the rare variants##
rare <- geno[which(geno$AF < 0.01 & geno$cadd_phred > 20 & geno$AC < 40),]#422624 x1251
raregeno <- rare[,6:ncol(rare)] 
raregeno <- as.matrix(raregeno)
raregeno <- t(raregeno)
dim(raregeno) #1246 422624

raregeno <- as.numeric(raregeno)
raregeno <- matrix(raregeno, nrow = 1246,ncol = 422624)

#########################
###SKAT with PCA of PS###
#########################
PS <- pcapop$x
PS <- PS[,1:10]

obj <- SKAT_Null_Model(y ~ PS , out_type = "D")

gene <- rare[,4]
genenames <- unique(gene) #18606

ppop <- rep(0,length(genenames))

for (i in 1:length(genenames)) {
  ppop[i] <- SKAT(matrix(raregeno[,c(which(gene == genenames[i]))],nrow = 1246),obj,missing_cutoff=0.2)
}
ppop <- as.numeric(ppop)


#############################
###SKAT without covariates###
#############################

obj <- SKAT_Null_Model(y ~ 1 , out_type = "D")
p <- rep(0,length(genenames))
for (i in 1:length(genenames)) {
  p[i] <- SKAT(matrix(raregeno[,c(which(gene == genenames[i]))],nrow = 1246),obj,missing_cutoff=0.2)
}
p <- as.numeric(p)

pvalues <- data.frame("p"=p, "ppop"=ppop)

#########################
####SKAT with kinship####
#########################

obj<-SKAT_NULL_emmaX(y ~ 1, K = grm)

pK <- rep(0,length(genenames))
for (i in 1:length(genenames)) {
  pK[i] <- SKAT(matrix(raregeno[,c(which(gene == genenames[i]))],nrow = 1246),obj,missing_cutoff=0.2)
}
length(which(pK<0.05)) #1569

pK <- as.numeric(pK)

pvalues <- cbind(pvalues, pK)


##################################
###SKAT with all the covariates###
##################################
PCpop <- loadingpop[,1:10]
obj <- SKAT_NULL_emmaX(y ~ PCpop, K = grm)

pall <- rep(0,length(genenames))
for (i in 1:length(genenames)) {
  pall[i] <- SKAT(matrix(raregeno[,c(which(gene == genenames[i]))],nrow = 1246),obj,missing_cutoff=0.2)
}


pall <- as.numeric(pall)

pvalues <- cbind(pvalues, pall)

write.csv(pvalues, "~/NAS/FECD_pvalues.csv")


############
###qqplot###
############
list1<-list("no_cov" = p,"Kin"= pK)

pdf(file = "~/fast-storage/FECD_qq_K.pdf")
qqunif.plot(list1, auto.key=list(corner=c(.95,.05)))
dev.off()


##########################
###get the top 10 genes###
##########################
result2 <- cbind(p,genenames)
result2 <- result2[order(result1[,1]),] 
length(which(p < 0.05)) #1203
write.csv(result2, file = "~/NAS/FECD_results2.csv")

##################
###Venn diagram###
##################
p <- cbind(pvalues$p,genenames)
pK <- cbind(pvalues$pK,genenames)
p <- p[which(pvalues$p < 0.05),]
pK <- pK[which(pvalues$pK < 0.05),]
noCov <- p[,2]
Kinship <- pK[,2]
x = list(no_Cov = noCov, Kinship=Kinship)

library(ggvenn)
pdf("~/fast-storage/FECD_Venn.pdf")
ggvenn(
  x, 
  fill_color = c("red","blue"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()


######################
###plot the results###
######################
plots <- data.frame("p"= p,"pK"= pK,"genes"=genenames)
plots$significant <- rep(0,18606)
plots$significant[which(plots$p<0.05 & plots$pK<0.05)] <- 1   #1065
plots <- plots[order(plots$p),]
plots$significant <- as.factor(plots$significant)

plots$p <- -log10(plots$p)
plots$pK <- -log10(plots$pK)

library(ggrepel)
pdf(file = "~/fast-storage/FECD_labelplot.pdf")
ggplot(plots,aes(x=p, y=pK,colour=significant, label = genenames))+
  geom_point()+
  geom_text_repel(aes(label=ifelse(p>=5.084629,as.character(genenames),'')))+
  geom_vline(xintercept = 1.3,linetype="dashed", color = "grey")+
  geom_hline(yintercept = 1.3,linetype="dashed",color = "grey")+
  xlim(0,20)

dev.off()

