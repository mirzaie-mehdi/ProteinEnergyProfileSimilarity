# This script generates the SPE and CPE based on PDBIDs from PDBIDs_Example.csv file.
# Th outputs are stored in the SPE and CPE data frames.
# ------------------------------------------------------------
# ---------- load requirement functions -------------
# All required functions are available at Functions.R script
source('Functions.R')
# ------------------------------------------------------------
library(bio3d)
library(geometry)
exData <- read.csv('Data/PDBIDs_Example.csv') # This file should include a column that contains the PDBIDs.
Nprotein <- nrow(exData)

SPE <- data.frame(exData,
                     matrix(NA, nrow = Nprotein, ncol = 210))

colnames(SPE)[-c(1:2)] <- pair_inter

CPE <- SPE


for (i in 1:Nprotein){
  print(i)
  pdbname <- substr(SPE$pdbID[i],1,4)
  ch <- substr(SPE$pdbID[i],5,5)
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
  # ---------------------------------------
  # ----------------------------------
  # ----- extract protein sequence from PDB file ------------
  seq <- new_pdb$seqres
  seq <- seq[names(seq)==ch]
  seq <- paste0(aa321(seq), collapse = '')
  # ----------------------------------
  # ----- CPE ------------
  aa210_p <- Energy_CPE210(seq)
  CPE[i,-c(1:2)] <- aa210_p
}

# -------------------------------------------------
# ------- visual total energy ---------
library(ggplot2)
library(ggpubr)

df <- data.frame(exData,
                 SPE=rowSums(SPE[,-c(1:3)]),
                 CPE=rowSums(CPE[,-c(1:3)])/2)

ggplot(df, aes(SPE, CPE, col=species))+
  geom_point(size=3, show.legend = F)+
  geom_text(aes(label = species), size=4, 
            hjust=ifelse(df$species=='Rat', .5, -.2),
            vjust=ifelse(df$species=='Rat', 1.5, .2), show.legend = F)+
  geom_smooth(method = 'glm', col='skyblue3', linewidth = .7, se = F)+
  stat_cor(method = "pearson", size=5,
           r.accuracy = 1e-2,
           color = "red", geom = "label")+
  theme_bw(base_line_size = .3)+
  ggtitle('Total Energy based on SPE and CPE')+
  xlab('Structural Profile of Energy (SPE)')+
  ylab('Compositional Profile of Energies (CPE)')
