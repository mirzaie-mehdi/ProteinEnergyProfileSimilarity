# This script generates the SPE based on PDBIDs from PDBIDs_Example.csv file.
# Th output is stored in the SPE data frame.
# ------------------------------------------------------------
# ---------- load requirement functions -------------
# All required functions are available at Functions.R script
source('Functions/Functions.R')
# ------------------------------------------------------------
library(bio3d)
library(geometry)
exData <- read.csv('Example/PDBIDs_Example.csv') # This file should include a column that contains the PDBIDs.
Nprotein <- nrow(exData)

SPE <- data.frame(exData,
                     matrix(NA, nrow = Nprotein, ncol = 210))

colnames(SPE)[-c(1:2)] <- pair_inter


for (i in 1:Nprotein){
  print(i)
  pdbname <- substr(SPE$pdbID[i],1,4)
  ch <- substr(SPE$pdbID[i],6,6)
  # ----------------------------------
  # ----- read pdb file ------------

  pdb0 <- read.pdb(pdbname,verbose = FALSE)
  # ----------------------------------
  # ----- filter ------------
  sele1 <- atom.select(pdb0,"protein", type = "ATOM", chain = ch)
  sele2 <- atom.select(pdb0, "protein", elety=c("H","OXT"),inverse=TRUE)
  sele4 <- atom.select(pdb0, "noh")
  sele <- combine.select(sele1, sele2, sele4, operator="AND")
  new_pdb <- trim.pdb(pdb0, sele,verbose = FALSE)
  # ----------------------------------
  # ----- Extract the atomic neighbors ------------
  Net_Frame <- NetworkFrame(PDB = new_pdb)
  # ----------------------------------
  # ----- SPE -----------
  energy210 <- Energy_SPE210(Net_Frame, energy_dell_dunbrack)
  SPE[i,-c(1:2)] <- energy210
}
# -------------------------------------------------
dis <- as.matrix(proxy::dist(SPE[,-c(1:2)], method = manhat))

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
  ggtitle('SPE')
