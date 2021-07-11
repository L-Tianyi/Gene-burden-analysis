install.packages("SKAT")

#1.Analysis with example data
```{r}
library(SKAT)
data(SKAT.example)
names(SKAT.example)
attach(SKAT.example)
# Z=Genotype(2000,67) y.c=ContinuousPhenotype y.b=BinaryPhenotype X=CovariatesMatrix
```

#1.1 SKAT
```{r}
#continuous trait
obj<-SKAT_Null_Model(y.c ~ X, out_type="C") #the null model for no associations
SKAT(Z,obj)$p.value #0.002877041
```
```{r}
#dichotomouse trait
obj<-SKAT_Null_Model(y.b ~ X, out_type="D")
SKAT(Z, obj)$p.value #0.1401991
```
```{r}
#Assign weights for each SNP
SKAT(Z, obj, kernel = "linear.weighted", weights.beta=c(0.5,0.5))$p.value #0.4931639
```

#1.2 SKAT-O
```{r}
SKAT(Z, obj, method="SKATO")$p.value #0.1008976
```

#1.3 CommonRare
```{r}
SKAT_CommonRare(Z, obj)$p.value #0.2238025
```

#1.4 Adjust for kinship
```{r}
data(SKAT.fam.example)
attach(SKAT.fam.example)
```
```{r}
obj<-SKAT_NULL_emmaX(y ~ X, K=K)
SKAT(Z, obj)$p.value #0.2123192
SKAT(Z, obj, method="SKATO")$p.value #0.352943
```
```{r}
detach(SKAT.fam.example)
```

#2.analyse with my onw data
```{r}
load("C:\\Users\\ltian\\Desktop\\BIOL0050\\assigment1\\snps.RData")
load("C:\\Users\\ltian\\Desktop\\BIOL0050\\assigment1\\mat.gtex.RData")
```
```{r}
genotypes <- snps
genotypes_round <- round(genotypes)
expression <- mat.gtex
PIK3CA <- expression[,1]
CDKN2A <- expression[,2]
TP53 <- expression[,3]
SMAD4 <- expression[,4]
covariates <- X
covariates <- covariates[1:53,2]
```
```{r}
obj1<-SKAT_Null_Model(PIK3CA ~ covariates, out_type="C")
SKAT(genotypes_round,obj1)$p.value #1.729861e-07
SKAT(genotypes_round,obj1,method = "SKATO")$p.value #2.523092e-07
SKAT_CommonRare(genotypes_round, obj1)$p.value #0.4095536
```
```{r}
obj2<-SKAT_Null_Model(CDKN2A ~ covariates, out_type="C")
SKAT(genotypes_round,obj2)$p.value #0.7777197
SKAT(genotypes_round,obj2,method = "SKATO")$p.value #0.3338262
SKAT_CommonRare(genotypes_round, obj2)$p.value #0.4815765
```
```{r}
obj3<-SKAT_Null_Model(TP53 ~ covariates, out_type="C")
SKAT(genotypes_round,obj3)$p.value #0.0001729987
SKAT(genotypes_round,obj3,method = "SKATO")$p.value #0.0003085336
SKAT_CommonRare(genotypes_round, obj3)$p.value #0.3797676
```
```{r}
obj4<-SKAT_Null_Model(SMAD4 ~ covariates, out_type="C")
SKAT(genotypes_round,obj4)$p.value #2.107267e-06
SKAT(genotypes_round,obj4,method = "SKATO")$p.value #3.451335e-06
SKAT_CommonRare(genotypes_round, obj4)$p.value #0.4945894
```

#2.1 Adjust for kinship
```{r}
install.packages("AGHmatrix")
```
```{r}
library(AGHmatrix)
```
```{r}
grm <- Gmatrix(genotypes_round)
grm_rd <- round(grm, digits = 2)
grm_rd <- pmax(grm_rd,0)
```
```{r}
obj1<-SKAT_NULL_emmaX(PIK3CA ~ covariates, K=grm_rd)
SKAT(genotypes_round, obj1)$p.value #2.141998e-07
SKAT(genotypes_round, obj1, method="SKATO")$p.value #3.927592e-07
```
```{r}
obj2<-SKAT_NULL_emmaX(CDKN2A ~ covariates, K=grm_rd)
SKAT(genotypes_round, obj2)$p.value #0.752798
SKAT(genotypes_round, obj2, method="SKATO")$p.value #0.3708531
```
```{r}
obj3<-SKAT_NULL_emmaX(TP53 ~ covariates, K=grm_rd)
SKAT(genotypes_round, obj3)$p.value #0.752798
SKAT(genotypes_round, obj3, method="SKATO")$p.value #0.3708531
```
```{r}
obj4<-SKAT_NULL_emmaX(SMAD4 ~ covariates, K=grm_rd)
SKAT(genotypes_round, obj4)$p.value #2.10725e-06
SKAT(genotypes_round, obj4, method="SKATO")$p.value #3.451309e-06
```
