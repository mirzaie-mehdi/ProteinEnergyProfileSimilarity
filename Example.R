source('Functions.R')
pdb_chain <- c('2dn1A','2dhbA','3gouA','3gdjA','3dhtA','2raoA','1fsxA')
sp <- c('Human', 'Horse','Dog','Camel','Rat','Rabbit','Cow')

Nprotein <- length(sp)

df_SPE <- data.frame(pdbID=pdb_chain, species= sp,
                     matrix(NA, nrow = Nprotein, ncol = 210))

colnames(df_SPE)[-c(1,2)] <- pair_inter

df_CPE <- df_SPE


for (i in 1:Nprotein){
  print(i)
  pdbname <- substr(df_SPE$pdbID[i],1,4)
  ch <- substr(df_SPE$pdbID[i],5,5)
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
  # ----- contact map ------------
  Net_Frame <- NetworkFrame(PDB = new_pdb)
  # ----------------------------------
  # ----- SPE ------------
  energy210 <- Energy_SPE210(Net_Frame, energy_dell_dunbrack)
  df_SPE[i,-c(1:2)] <- energy210
  # ---------------------------------------
  # ----------------------------------
  # ----- extract sequence ------------
  seq <- new_pdb$seqres
  seq <- seq[names(seq)==ch]
  seq <- paste0(aa321(seq), collapse = '')
  # ----------------------------------
  # ----- CPE ------------
  aa210_p <- Energy_CPE210(seq)
  df_CPE[i,-c(1:2)] <- aa210_p
}
