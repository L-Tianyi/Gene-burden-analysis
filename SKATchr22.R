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
pcapop <- chrM[which(missing == 0),] #kep the variants that have a full coverage
pcapop <- pcapop[,3:ncol(pcapop)]
pcapop_t <- t(pcapop)
pcapop_t <- as.numeric(pcapop_t) #change into numeric because prcomp only accept numeric data
pcapop_t <- matrix(pcapop_t,nrow = 1001, ncol = 36644,byrow = FALSE)
pcapop_t <- pcapop_t[ , which(apply(pcapop_t, 2, var) != 0)] #remove the columns have 0 variance
pcpop <- prcomp(pcapop_t, scale. = TRUE)

library(ggplot2)
library(factoextra)

fviz_eig(pcpop) #get scree plot
loading = pcpop$x
loading <- as.data.frame(loading)
ggplot(loading, aes(x = PC1, y = PC2)) +
  geom_point(alpha=0.5) + 
  labs(x = "PC1", y = "PC2")+
  theme_classic()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

#PCA for missing#
chrM_t <- t(chrM)
miss <- chrM_t[3:nrow(chrM_t),]
miss <- as.numeric(miss)
miss <- matrix(miss, nrow = 1001, ncol = 72554,byrow = FALSE)
for (i in 1:ncol(miss)) {
  miss[which(miss[,i] == 0),i] <- 1
}
for (i in 1:ncol(miss)) {
  miss[which(miss[,i] == 1),i] <- 1
}
for (i in 1:ncol(miss)) {
  miss[which(miss[,i] == 2),i] <- 1
}
for (i in 1:ncol(miss)) {
  miss[which(is.na(miss[,i])),i] <- 0
}

miss <- miss[ , which(apply(miss, 2, var) != 0)]
pctech <- prcomp(miss, scale. = TRUE)

loading = pctech$x
loading <- as.data.frame(loading)
experiments <- c("4216",rep("B",16),rep("BGI",127),rep("Blizard",7),rep("Broad",29),rep("Brogan",5),rep("CEGC",28),rep("Czech",35),
                 rep("Davina",48),rep("Hardcastle",15),rep("IRDC",128),rep("IoN",7),rep("IoO",208),rep("Kaoru",5),rep("kelsell",6),
                 "kinki","kinki",rep("LIMM",9),"Levine",rep("Nejentsev",30),"Neringa",rep("OneKG",4),rep("QMUL",32),rep("A",7),
                 rep("C",4),rep("D",8),rep("Shamima",16),rep("Sisodiya",5),rep("UCL",7),rep("VanitaBerry",36),rep("Villiamy",172),
                 "gosgene","gosgene") #creat a vector for experiments#
loading = cbind(loading,experiments)
loading$experiments <- as.factor(loading$experiments)
ggplot(loading, aes(x = PC1, y = PC2, color = experiments)) +
  geom_point(alpha=0.5) + 
  labs(x = "PC1", y = "PC2")


##simulate phenotype and covariates##
#choose a subset of causal variants#
maf <- chrM[,1]
maf <- as.numeric(maf)
rare <- which(maf <= 0.05)
causal <- sample(rare,0.05*nrow(chr22)) #choose 5% of rare variants (maf<0.05) from the sample as causal variants
causal_variants <- chrM[causal, 3:ncol(chrM)]
causal_variants <- as.numeric(causal_variants)
causal_variants <- matrix(causal_variants, nrow = 3627, ncol = 1001,byrow = FALSE)
#impute missing values in the subset of variants#
causal_maf <- chrM[causal, 1]
causal_maf <- as.numeric(causal_maf)
causal_NA <- rep(0,length(causal))
for (i in 1:length(causal_NA)) {
  causal_NA[i] <- length(which(is.na(causal_variants[i,])))
}
for (i in 1:length(causal)) {
  causal_variants[i,which(is.na(causal_variants[i,]))] <- rbinom(causal_NA[i],2,causal_maf[i])
}
##get effect size (beta)#
beta <- causal_maf
beta <- abs(log(beta))
#simulate covariates
age <- sample(10:60,1001,replace = TRUE)
sex <- rbinom(1001,1,0.5)

#simulate phenotypes using linear combination#
y <- rep(0,1001)
for (i in 1:1001) {
  y[i] = 1 + sum(beta*causal_variants[,i])+ 0.2*age[i] + 0.2*sex[i]
}
summary(y)
##turn y into binary
mediany <- mediant(y)
y[which(y < mediany)] <- 0
y[which(y >= mediany)] <- 1



###SKAT analysis###
#create vactors and matrix needed for the analysis#
gene <- chrM[,2]
genenames <- unique(chrM[,2])
p <- rep(0,length(genenames))

chrSKAT <- t(chrM[,3:ncol(chrM)])
chrSKAT <- as.numeric(chrSKAT)
chrSKAT <- matrix(chrSKAT, nrow = 1001, ncol = 72554, byrow = FALSE)

#SKAT1#
obj <- SKAT_Null_Model(y ~ sex + age, out_type = "D")
for (i in 1:length(genenames)) {
  p[i] <- SKAT(matrix(chrSKAT[,c(which(gene == genenames[i]))],nrow = 1001),obj,missing_cutoff=0.2)
}
p <- as.numeric(p)

#SKAT2 (with PCtech)#
Tech <- pctech$x[,1:6] #choose first few PCs
obj <- SKAT_Null_Model(y ~ sex + age + Tech, out_type = "D")

pT <- rep(0,length(genenames))
for (i in 1:length(genenames)) {
  pT[i] <- SKAT(matrix(chrSKAT[,c(which(gene == genenames[i]))],nrow = 1001),obj,missing_cutoff=0.2)
}
pT <- as.numeric(pT)

p_values <- cbind(p,pT)

#creat Kinship matrix#
missing <- rep(0,nrow(chrM)) 
for (i in 1:length(missing)) {
  missing[i] <- length(which(is.na(chrM[i,])))
}
chrKinship <- chrM[which(missing == 0),] #36644x1003 #use the well-covered variants to calculate kinship
chrKinship <- chrKinship[,3:ncol(chrKinship)]
chrKinship <- t(chrKinship) #1001x36644
chrKinship <- as.numeric(chrKinship)
chrKinship <- matrix(chrKinship, nrow = 1001, ncol = 36644,byrow = FALSE)

library(AGHmatrix)
grm <- Gmatrix(chrKinship, method = "VanRaden")

##SKAT3 with Kinship##
obj<-SKAT_NULL_emmaX(y ~ age + sex, K = grm)
pK <- rep(0,length(genenames))
for (i in 1:length(genenames)) {
  pK[i] <- SKAT(matrix(chrSKAT[,c(which(gene == genenames[i]))],nrow = 1001),obj,missing_cutoff=0.2)
}
pK <- as.numeric(pK)

p_values <- cbind(p_values,pK)


###make qqplot###
#set a function, codes are from https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R##
library(lattice)
qqunif.plot<-function(pvalues, 
	should.thin=T, thin.obs.places=2, thin.exp.places=2, 
	xlab=expression(paste("Expected (",-log[10], " p-value)")),
	ylab=expression(paste("Observed (",-log[10], " p-value)")), 
	draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
	already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
	par.settings=list(superpose.symbol=list(pch=pch)), ...) {
	
	#error checking
	if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
	if(!(class(pvalues)=="numeric" || 
		(class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
		stop("pvalue vector is not numeric, can't draw plot")
	if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
	if (already.transformed==FALSE) {
		if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
	} else {
		if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
	}
	
	grp<-NULL
	n<-1
	exp.x<-c()
	if(is.list(pvalues)) {
		nn<-sapply(pvalues, length)
		rs<-cumsum(nn)
		re<-rs-nn+1
		n<-min(nn)
		if (!is.null(names(pvalues))) {
			grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
			names(pvalues)<-NULL
		} else {
			grp=factor(rep(1:length(pvalues), nn))
		}
		pvo<-pvalues
		pvalues<-numeric(sum(nn))
		exp.x<-numeric(sum(nn))
		for(i in 1:length(pvo)) {
			if (!already.transformed) {
				pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
				exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
			} else {
				pvalues[rs[i]:re[i]] <- pvo[[i]]
				exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
			}
		}
	} else {
		n <- length(pvalues)+1
		if (!already.transformed) {
			exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
			pvalues <- -log10(pvalues)
		} else {
			exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
		}
	}

	#this is a helper function to draw the confidence interval
	panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
		require(grid)
		conf.points = min(conf.points, n-1);
		mpts<-matrix(nrow=conf.points*2, ncol=2)
        	for(i in seq(from=1, to=conf.points)) {
            		mpts[i,1]<- -log10((i-.5)/n)
            		mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
            		mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
            		mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
        	}
        	grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
    	}

	#reduce number of points to plot
	if (should.thin==T) {
		if (!is.null(grp)) {
			thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
				exp.x = round(exp.x, thin.exp.places),
				grp=grp))
			grp = thin$grp
		} else {
			thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
				exp.x = round(exp.x, thin.exp.places)))
		}
		pvalues <- thin$pvalues
		exp.x <- thin$exp.x
	}
	gc()
	
	prepanel.qqunif= function(x,y,...) {
		A = list()
		A$xlim = range(x, y)*1.02
		A$xlim[1]=0
		A$ylim = A$xlim
		return(A)
	}

	#draw the plot
	xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
		prepanel=prepanel, scales=list(axs="i"), pch=pch,
		panel = function(x, y, ...) {
			if (draw.conf) {
				panel.qqconf(n, conf.points=conf.points, 
					conf.col=conf.col, conf.alpha=conf.alpha)
			};
			panel.xyplot(x,y, ...);
			panel.abline(0,1);
		}, par.settings=par.settings, ...
	)
}

##make my qqplots##
my.pvalue.list1<-list("sex+age" = p_values[,1],"sex+age+Tech"= p_values[,2])
my.pvalue.list2<-list("sex+age"= p_values[,1],"sex+age+Kin" = p_values[,3])

qqunif.plot(my.pvalue.list1, auto.key=list(corner=c(.95,.05)))
qqunif.plot(my.pvalue.list2, auto.key=list(corner=c(.95,.05)))
