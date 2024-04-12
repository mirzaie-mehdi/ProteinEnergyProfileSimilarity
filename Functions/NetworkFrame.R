NetworkFrame <- function(PDB){  #, AtomSelect
  Tdist <- 5; start <- 0.75; binw <- 0.5; Tseqsep <- 1
  pdb <- PDB$atom#[AtomSelect$atom,]
  resnum <- length(unique(pdb$resno))
  xyz3 <- PDB$xyz#[,AtomSelect$xyz]
  coordinate <- matrix(xyz3,ncol=3,byrow = TRUE)
  atnumber <- pdb$eleno
  atomname <- pdb$elety
  amacidname <- pdb$resid
  amnumber <- pdb$resno
  #***********Contact extraction based on Delaunay Triangulation********
  DT <- delaunayn(coordinate)
  x12 <- DT[,1:2]; x13 <- DT[,c(1,3)]; x14 <- DT[,c(1,4)]
  x23 <- DT[,2:3]; x24 <- DT[,c(2,4)]; x34 <- DT[,3:4]
  edges <-rbind(x12,x13,x14,x23,x24,x34)
  edges <- edges[!duplicated(edges),]# remove duplicated edges
  cs <- abs(amnumber[edges[,1]]-amnumber[edges[,2]])
  edges <- edges[cs>=Tseqsep,] #remove contacts within sequence seperation less than 1
  cs <- cs[cs>=Tseqsep]
  dis <- sqrt(rowSums((coordinate[edges[,1],]-coordinate[edges[,2],])^2))
  #**********************************************
  # remove contacts with distance > Tdist
  edges <- edges[dis<= Tdist,];cs <- cs[dis<= Tdist];dis <- dis[dis<= Tdist]
  Ff <- trunc((dis-start)/binw)+1;Ff[Ff<1]<- 1
  # remove interaction of peptide bond
  at1 <- atomname[edges[,1]]
  at2 <- atomname[edges[,2]]
  ind <- !((at1 =='N' & at2 =='C' & cs ==1)|(at2 =='N' & at1 =='C' & cs ==1))
  seqsep <- cs[ind]
  distance <- dis[ind]
  distance_bin <- Ff[ind]
  edge_list <- edges[ind,]
  ati <- at1[ind];atj <- at2[ind]
  ami <- amacidname[edge_list[,1]];amj <- amacidname[edge_list[,2]]
  resnoi <- amnumber[edge_list[,1]];resnoj <- amnumber[edge_list[,2]]
  # This section creates a data.frame including contact information
  #Net_Frame <- data.frame(edge_list,distance,distance_bin,ati,atj,ami,amj,seqsep,resnoi,resnoj)
  edge_list2 <- cbind(edge_list[,2],edge_list[,1])
  edge_listf <- rbind(edge_list,edge_list2)
  distancef <- c(distance,distance)
  distance_binf <- c(distance_bin,distance_bin)
  atif <- c(ati, atj)
  atjf <- c (atj, ati)
  amif <- c(ami, amj)
  amjf <- c (amj, ami)
  seqsepf <- c(seqsep,seqsep)
  resnoif <- c(resnoi,resnoj)
  resnojf <- c(resnoj,resnoi)
  Net_Frame2 <- data.frame(edge_listf,distancef,distance_binf,atif,atjf,amif,amjf,seqsepf,resnoif,resnojf)
  Net_Frame <- Net_Frame2[!duplicated(edge_listf),]
  return(Net_Frame)
}
