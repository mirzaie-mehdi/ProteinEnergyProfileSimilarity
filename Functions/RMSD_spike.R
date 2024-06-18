# RMSD_spike
#This script computes the RMSD for spike glycoprotein proteins. The PDB IDs of these proteins are located in the 'Data/csv/spike.csv' file.
# ----------------------------------------------------------
spike <- read.csv("Data/csv/spike.csv")
Np <- nrow(spike)
#----------------------------------------------
rms <- matrix(NA,Np,Np)

tictoc::tic()
for (i in 1:Np) {
  print(i)
  tryCatch({
    chi <- substr(spike$pdbID[i],5,5)
    pdbnamei <- substr(spike$pdbID[i],1,4)
    pdb0 <- read.pdb(tolower(pdbnamei))
    sele <- atom.select(pdb0,chain=chi)
    new_pdbi <- trim.pdb(pdb0, sele,verbose = FALSE)
    for (j in 1:Np) {
      if(i>j){
        chj <- substr(spike$pdbID[j],5,5)
        pdbnamej <- substr(spike$pdbID[j],1,4)
        pdb0 <- read.pdb(tolower(pdbnamej))
        sele <- atom.select(pdb0,'protein',chain=chj)
        new_pdbj <- trim.pdb(pdb0, sele,verbose = FALSE)
        pdbs <- list(new_pdbi, new_pdbj)
        aln_pdbs <- pdbaln(pdbs,fit = T)
        xyz <- pdbfit(aln_pdbs)
        rms[i,j] <- rms[j,i] <- rmsd(xyz)[1,2]
      }
    }
  }, error = function(e) {
    cat("Error downloading file:", e$message, "\n")
  })
}
file.remove('aln.fa')
tictoc::toc()
rownames(rms) <- colnames(rms) <- spike$pdbID
saveRDS(rms,'Data/covidPDB/RMSD_spike.rds')
