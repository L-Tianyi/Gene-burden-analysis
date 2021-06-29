#Language:R#
#load the data file#
Genotype <- read.csv("C:\\Users\\ltian\\Desktop\\GENE0014\\joint_genotyping_anno_chr22.csv")


#Get the information of the data file#
dim(Genotype)
colnames(Genotype)
Genotype[1:20,1:33] #the first 33 columns are the properties of the variants
#cadd = combined annotation dependent depletion
#dann = deleterious annotation of genetic variants using neural networks
#tommo = Tohoku Medical Megabank Organization
#KRGDB = Korean Reference Genome Database
#hgvd = Human genetic variation database


#Make some plots#
library(ggplot2)

cadd <- ggplot(Genotype, aes(x = cadd_phred, y = AF)) + geom_point(color = "blue") + labs(x="CADD Score", y = "Allele frequency")
cadd
#allele frequencies against CADD score#

dann <- ggplot(Genotype, aes(x = dann, y = AF)) + geom_point(color = "red") + labs(x="DANN Score", y = "Allele frequency")
dann
#allele frequencies against DANN score#

table(Genotype$most_severe_consequence)
genotype_rm <- Genotype[Genotype$most_severe_consequence != "?",c(21,33)]
box <- ggplot(genotype_rm, aes(x=most_severe_consequence, y = cadd_phred, color = most_severe_consequence)) + geom_boxplot()+ labs(x = "consequences", y = "CADD Score")
box
#the distribution of CADD score of each category of the consequence of genetic variation#
