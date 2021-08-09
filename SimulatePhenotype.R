##simulate a binary phenotype##
#choose a subset of causal variants#
maf <- chrM[,2]
maf <- as.numeric(maf)
rare <- which(maf <= 0.05)
causal <- sample(rare,0.05*nrow(chr22)) #choose 5% of rare variants (maf<0.05) from the sample as causal variants#
causal_variants <- chrM[causal, 3:ncol(chrM)]
causal_variants <- as.numeric(causal_variants)
causal_variants <- matrix(causal_variants, nrow = 3627, ncol = 1001,byrow = FALSE)

#impute missing values in the subset of variants#
causal_maf <- chrM[causal, 2]
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
##turn Z into binary
mediany <- mediant(y)
y[which(y < mediany)] <- 0
y[which(y >= mediany)] <- 1
