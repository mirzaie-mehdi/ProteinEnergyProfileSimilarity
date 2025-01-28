# This script generates the CPE based on Sequence from Sequence_Example.csv file.
# Th output is stored in the CPE data frame.
# ------------------------------------------------------------
# ---------- load requirement functions -------------
# All required functions are available at Functions.R script
source('Functions.R')
# ------------------------------------------------------------
library(bio3d)
library(geometry)
exData <- read.csv('Example/Sequence_Example.csv') # This file should include a column that contains the Sequences.
Nprotein <- nrow(exData)

CPE <- data.frame(exData,
                  matrix(NA, nrow = Nprotein, ncol = 210))

colnames(CPE)[-c(1:2)] <- pair_inter


for (i in 1:Nprotein){
  print(i)
  seq <- CPE$seq[i]
  aa210_p <- Energy_CPE210(seq)
  CPE[i,-c(1:2)] <- aa210_p
}
# -------------------------------------------------
dis <- as.matrix(proxy::dist(CPE[,-c(1:2)], method = manhat))

ump <- umap(d = dis,
            n_neighbors = 13,
            min_dist = 0.1,
            input="dist")

ump <- data.frame(exData, 
                  UMAP1=ump$layout[,1],
                  UMAP2=ump$layout[,2])

library(ggplot2)
ggplot(ump, aes(UMAP1, UMAP2, color= globin)) +
  geom_point(alpha=.7,size=3) +
  theme_bw() +
  ggtitle('CPE')
