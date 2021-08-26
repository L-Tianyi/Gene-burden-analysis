#####################################
####simulate genotype matrices(G)####
#####################################

###G1####
G <- matrix(0, nrow = 600, ncol = 2000)
for (i in 1:600) {
  G[i,] = rbinom(2000,2,0.005)
}

###G2###
G <- matrix(0, nrow = 600, ncol = 2000)
for (i in 1:600) {
  G[i,] = rbinom(2000,2,0.025)
}


###G3###
G <- matrix(0, nrow = 600, ncol = 2000)
for (i in 1:200) {
    G[i,] = rbinom(2000,2,0.05)
}
for (i in 201:400) {
  G[i,] = rbinom(2000,2,0.025)
}
for (i in 401:600) {
  G[i,] = rbinom(2000,2,0.005)
}

###G4###
G <- matrix(0, nrow = 600, ncol = 20000)
for (i in 1:200) {
    G[i,] = rbinom(20000,2,0.005)
}
for (i in 201:400) {
  G[i,] = rbinom(20000,2,0.0025)
}
for (i in 401:600) {
  G[i,] = rbinom(20000,2,0.0005)
}

###G5###
G <- matrix(0, nrow = 600, ncol = 2000)
for (i in 1:150) {
    G[i,] = rbinom(2000,2,0.005)
}
for (i in 151:300) {
  G[i,] = rbinom(2000,2,0.0025)
}
for (i in 301:450) {
  G[i,] = rbinom(2000,2,0.0005)
}
for (i in 451:600) {
  G[i,] = rbinom(2000,2,0.025)
}

###G6###
G <- matrix(0, nrow = 600, ncol = 2000)
for (i in 1:600) {
  G[i,] = rbinom(2000,2,0.005)
}

###############################
###create missingness matrix###
###############################

###simulate missing genptypes in G6###
G_M <- G
for (i in 1:200) {
  for (j in 1:2000) {
  x = rbinom(1,1,0.01)
  if (x == 1) G_M[i,j] <- -1
  }
}
for (i in 201:400) {
  for (j in 1:2000) {
  x = rbinom(1,1,0.02)
  if (x == 1) G_M[i,j] <- -1
  }
}
for (i in 401:600) {
  for (j in 1:2000) {
  x = rbinom(1,1,0.05)
  if (x == 1) G_M[i,j] <- -1
  }
}

##create a minssingness matrix##
miss <- G_M
miss[miss == 0] <- 1
miss[miss == 1] <- 1
miss[miss == 2] <- 1
miss[miss == -1] <- 0

#########
###PCA###
#########
G_rm <- G[ , which(apply(G, 2, var) != 0)] #remove the variants that does not have variation

G_pca <- prcomp(G_rm,scale. = TRUE) #PCA with genotype matrices
miss_pca <- prcomp(miss,scale. = TRUE) #PCA with the missingness matrix

##############
###PCA plot###
##############
library(ggplot2)
library(factoextra)

fviz_eig(G_pca) #generate the scree plot using the "factoextra" package
fviz_eig(miss_pca)

##plot for G1 and G2##
G_loading = G_pca$x
G_loading <- as.data.frame(G_loading)
ggplot(G_loading, aes(x = PC1, y = PC2 )) +
  geom_point(alpha=0.5) +
  labs(x = "PC1", y = "PC2")+
  theme_classic()+
  scale_color_brewer(palette="Set1")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

##plot for G3##
G_loading = pca$x
#create a column of minor allele frequency (MAF)
G_loading = cbind(G_loading,c(rep(0.1,200),rep(0.05,200),rep(0.01,200))) 
colnames(G_loading) <- c(colnames(pca$x),"MAF")
G_loading <- as.data.frame(G_loading, col.names = colnames(G_loading))
G_loading$MAF <- as.factor(G_loading$MAF)

ggplot(G_loading, aes(x = PC1, y = PC2, color = MAF)) +
  geom_point(alpha=0.5) + 
  labs(x = "PC1", y = "PC2")+
  theme_classic()+
  scale_color_brewer(palette="Set1")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

##plot for G4##
G_loading = pca$x
G_loading = cbind(G_loading,c(rep(0.01,200),rep(0.005,200),rep(0.001,200)))
colnames(G_loading) <- c(colnames(pca$x),"MAF")
G_loading <- as.data.frame(G_loading, col.names = colnames(G_loading))
G_loading$MAF <- as.factor(G_loading$MAF)

ggplot(G_loading, aes(x = PC1, y = PC2, color = MAF)) +
  geom_point(alpha=0.5) + 
  labs(x = "PC1", y = "PC2")+
  theme_classic()+
  scale_color_brewer(palette="Set1")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
        
##plot for G5##
G_loading = pca$x
G_loading = cbind(G_loading,c(rep(0.01,150),rep(0.005,150),rep(0.001,150),rep(0.05,150)))
colnames(G_loading) <- c(colnames(pca$x),"MAF")
G_loading <- as.data.frame(G_loading, col.names = colnames(G_loading))
G_loading$MAF <- as.factor(G_loading$MAF)

ggplot(G_loading, aes(x = PC1, y = PC2, color = MAF)) +
  geom_point(alpha=0.5) + 
  labs(x = "PC1", y = "PC2")+
  theme_classic()+
  scale_color_brewer(palette="Set1")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

##plot for missingness##
G_loading = pca$x
G_loading = cbind(G_loading,c(rep(0.05,200),rep(0.1,200),rep(0.15,200)))
colnames(G_loading) <- c(colnames(pca$x),"Missing")
G_loading <- as.data.frame(G_loading, col.names = colnames(G_loading))
G_loading$Missing <- as.factor(G_loading$Missing)

ggplot(G_loading, aes(x = PC1, y = PC2, color = Missing)) +
  geom_point(alpha=0.5) + 
  labs(x = "PC1", y = "PC2")+
  theme_classic()+
  scale_color_brewer(palette="Set1")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))




