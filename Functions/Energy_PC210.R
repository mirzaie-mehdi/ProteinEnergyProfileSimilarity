Energy_PC210 <- function(NetworkFrame,energy_dell_dunbrack){ # Pairwise Contact potential
  ## ***************read file trained energy**************
  #energy_dell_dunbrack <- read.csv("C:/Users/Choopanian/Desktop/structure_energy/energy.csv",header = FALSE,sep ="" )
  #energy_dell_dunbrack <- energy_dell_dunbrack[,5:33]
  # *****************Energy file Indexing**************
  cn <- 1; index167 <- matrix(0,ncol = 167,nrow = 167)
  for (i in 1:167){
    for (j in i: 167){
      index167[i,j] <- cn
      cn <- cn+1
    }
  }
  cn <- 1; index20 <- matrix(0,ncol = 20,nrow = 20)
  for (i in 1:20){
    for (j in i: 20){
      index20[i,j] <- cn
      cn <- cn+1
    }
  }
  # ***************atomtype2num calculation***********************
  v <- function(Atom,AminoAcid) {  
    atomtype2num167(as.character(Atom),as.character(AminoAcid))
  }
  # ***********
  atomi <- NetworkFrame$atif
  amacidi <- NetworkFrame$amif
  numi <- mapply(v,atomi,amacidi)
  # ***********
  atomj <- NetworkFrame$atjf
  amacidj <- NetworkFrame$amjf
  numj <- mapply(v,atomj,amacidj)
  # **************************************************************
  nrf <- nrow(NetworkFrame)
  energy <- rep(0,nrf)
  for (i in 1:nrf){
    vi <- numi[[i]];vj <- numj[[i]]
    if( vi > 0 & vj > 0){
      bin <- NetworkFrame$distance_binf[i];
      s1 <- min(vi,vj); t1 <- max(vi,vj)
      ei <- index167[s1,t1]
      energy[i] <- energy_dell_dunbrack[ei,bin]
    }
  }
  NetworkFrame <- data.frame(NetworkFrame,energy)
  es <- matrix(0, nrow = 20, ncol = 20)
  ami <- NetworkFrame$amif;amj <- NetworkFrame$amjf
  am2numi <- sapply(as.character(ami), amacid2num)
  am2numj <- sapply(as.character(amj), amacid2num)
  NetworkFrame <- data.frame(NetworkFrame,am2numi,am2numj,numi,numj)
  energy_pair_am <- aggregate(NetworkFrame$energy,by=list(NetworkFrame$am2numi,NetworkFrame$am2numj),FUN=sum, na.rm=TRUE)
  colnames(energy_pair_am) <- c("aminoacid1","aminoacid2",'energy')
  energy210 <- rep(0,210); ninter <- dim(energy_pair_am)[1]
  for (i in 1:ninter){
    r <-energy_pair_am$aminoacid1[i]; c <-energy_pair_am$aminoacid2[i]
    rw=min(r,c); cl=max(r,c)
    ind210 <- index20[rw,cl]
    energy210[ind210] <- energy_pair_am$energy[i]
  }
  return(energy210)
}
