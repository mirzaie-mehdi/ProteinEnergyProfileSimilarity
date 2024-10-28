# -----------------------------------------------------------------------------
rm(list = ls())
# -----------------------------------
# -------------------------------------------------------------
# ------------------------------
# ---- Load Knowledge_based potential from structure and sequence
# --------------------------------------------
# --------------------------------------------------------
# -------------------------------------------------------------------------------------
energy_dell_dunbrack <- read.csv("Data/csv/energy_rm_Olaps.csv",header = FALSE,sep ="," )
energy_dell_dunbrack <- energy_dell_dunbrack[,5:33]
# ------Pairwise Amino Acid names-------------
AA <- c("PHE","LEU","ILE","VAL","TRP","TYR","MET","CYS","HIS","THR",
        "ARG","ALA","ASN","GLN","PRO","SER","ASP","GLY","LYS","GLU")
letters_list <- bio3d::aa321(AA)
pair_inter <- character(0)
for (i in 1:20) {
  for (j in i:20) {
    pair_inter <- c(pair_inter, paste0(letters_list[i], letters_list[j]))
  }
} 
# -----------------------------------
# -----  Sequence Energy Estimator (predictor)
aaenergy <- read.csv("Data/csv/Pij_rm_Olaps.csv", header = F,sep = ";")
aaenergy <- (aaenergy + t(aaenergy))/2
colnames(aaenergy) <- rownames(aaenergy) <- letters_list
# -------------------------------------------------------------
rm(i,j,AA)
# -----------
# -------------------------
# ------ This function code Amino Acid name into number.
# ------------------------------------------------------------
amacid2num <- function(amname){
  v <- switch (amname,
       "PHE" = 1,
       "LEU" = 2,
       "ILE" = 3,
       "VAL" = 4,
       "TRP" = 5,
       "TYR" = 6,
       "MET" = 7,
       "CYH" = 8,
       "CYS" = 8,
       "CYSS" =8,
       "HIS" = 9,
       "THR" = 10,
       "ARG" = 11,
       "ALA" = 12,
       "ASN" = 13,
       "GLN" = 14,
       "PRO" = 15,
       "SER" = 16,
       "ASP" = 17,
       "GLY" = 18,
       "LYS" = 19,
       "LYZ" = 19,
       "GLU" = 20,
       21
       )
  v
}
# -------------------------------------------------
# -----------------------
# ----------- This function code pair of Atom types into numbers.
# -----------------------------------------------------
atomtype2num167 <- function(atname,amname){
  v <- switch (amname,
               "GLY" = switch (atname,
                               "N" = 1,
                               "CA" = 2,
                               "CB" = 2,
                               "C" = 3,
                               "O" = 4,
                               0),
               "ALA" = switch (atname,
                               "N" = 5,
                               "CA" = 6,
                               "C" = 7,
                               "O" = 8,
                               "CB" = 9,
                               0),
               "VAL" = switch (atname,
                               "N" = 10,
                               "CA" = 11,
                               "C" = 12,
                               "O" = 13,
                               "CB" = 14,
                               "CG1" = 15,
                               "CG2"= 16,
                               0),
               "LEU" = switch (atname,
                               "N" = 17,
                               "CA" = 18,
                               "C" = 19,
                               "O" = 20,
                               "CB" = 21,
                               "CG" = 22,
                               "CD1" = 23,
                               "CD2" = 24,
                               0),
               "ILE" = switch (atname,
                               "N" = 25,
                               "CA" = 26,
                               "C" = 27,
                               "O" = 28,
                               "CB" = 29,
                               "CG1" = 30,
                               "CG2" = 31,
                               "CD1" = 32,
                               0),
               "MET" = switch (atname,
                               "N" = 33,
                               "CA" = 34,
                               "C" = 35,
                               "O" = 36,
                               "CB" = 37,
                               "CG" = 38,
                               "SD" = 39,
                               "CE" = 40,
                               0),
               "PRO" = switch (atname,
                               "N" = 41,
                               "CA" = 42,
                               "C" = 43,
                               "O" = 44,
                               "CB" = 45,
                               "CG" = 46,
                               "CD" = 47,
                               0),
               "HIS" = switch (atname,
                               "N" = 48,
                               "CA" = 49,
                               "C" = 50,
                               "O" = 51,
                               "CB" = 52,
                               "CG" = 53,
                               "ND1" = 54,
                               "CD2" = 55,
                               "CE1" = 56,
                               "NE2" = 57,
                               0),
               "PHE" = switch (atname,
                               "N" = 58,
                               "CA" = 59,
                               "C" = 60,
                               "O" = 61,
                               "CB" = 62,
                               "CG" = 63,
                               "CD1" = 64,
                               "CD2" = 65,
                               "CE1" = 66,
                               "CE2" = 67,
                               "CZ" = 68,
                               0),
               "TYR" = switch (atname,
                               "N" = 69,
                               "CA" = 70,
                               "C" = 71,
                               "O" = 72,
                               "CB" = 73,
                               "CG" = 74,
                               "CD1" = 75,
                               "CD2" = 76,
                               "CE1" = 77,
                               "CE2" = 78,
                               "CZ" = 79,
                               "OH" = 80,
                               0),
               "TRP" = switch (atname,
                               "N" = 81,
                               "CA" = 82,
                               "C" = 83,
                               "O" = 84,
                               "CB" = 85,
                               "CG" = 86,
                               "CD1" = 87,
                               "NE1" = 88,
                               "CD2" = 89,
                               "CE2" = 90,
                               "CE3" = 91,
                               "CZ3" = 92,
                               "CZ2" = 93,
                               "CH2" = 94,
                               0),
               "CYS" = switch (atname,
                               "N" = 95,
                               "CA" = 96,
                               "C" = 97,
                               "O" = 98,
                               "CB" = 99,
                               "SG" = 100,
                               0),
               "CYSS" = switch (atname,
                                "N" = 95,
                                "CA" = 96,
                                "C" = 97,
                                "O" = 98,
                                "CB" = 99,
                                "SG" = 100,
                                0),
               "CYH" = switch (atname,
                               "N" = 95,
                               "CA" = 96,
                               "C" = 97,
                               "O" = 98,
                               "CB" = 99,
                               "SG" = 100,
                               0),
               "SER" = switch (atname,
                               "N" = 101,
                               "CA" = 102,
                               "C" = 103,
                               "O" = 104,
                               "CB" = 105,
                               "OG" = 106,
                               0),
               "THR" = switch (atname,
                               "N" = 107,
                               "CA" = 108,
                               "C" = 109,
                               "O" = 110,
                               "CB" = 111,
                               "OG1" = 112,
                               "CG2" = 113,
                               0),
               "ASN" = switch (atname,
                               "N" = 114,
                               "CA" = 115,
                               "C" = 116,
                               "O" = 117,
                               "CB" = 118,
                               "CG" = 119,
                               "OD1" = 120,
                               "ND2" = 121,
                               0),
               "GLN" = switch (atname,
                               "N" = 122,
                               "CA" = 123,
                               "C" = 124,
                               "O" = 125,
                               "CB" = 126,
                               "CG" = 127,
                               "CD" = 128,
                               "OE1" = 129,
                               "NE2" = 130,
                               0),
               "ASP" = switch (atname,
                               "N" = 131,
                               "CA" = 132,
                               "C" = 133,
                               "O" = 134,
                               "CB" = 135,
                               "CG" = 136,
                               "OD1" = 137,
                               "OD2" = 138,
                               0),
               "GLU" = switch (atname,
                               "N" = 139,
                               "CA" = 140,
                               "C" = 141,
                               "O" = 142,
                               "CB" = 143,
                               "CG" = 144,
                               "CD" = 145,
                               "OE1" = 146,
                               "OE2" = 147,
                               0),
               "LYS" = switch (atname,
                               "N" = 148,
                               "CA" = 149,
                               "C" = 150,
                               "O" = 151,
                               "CB" = 152,
                               "CG" = 153,
                               "CD" = 154,
                               "CE" = 155,
                               "NZ" = 156,
                               0),
               "LYZ" = switch (atname,
                               "N" = 148,
                               "CA" = 149,
                               "C" = 150,
                               "O" = 151,
                               "CB" = 152,
                               "CG" = 153,
                               "CD" = 154,
                               "CE" = 155,
                               "NZ" = 156,
                               0),
               "ARG" = switch (atname,
                               "N" = 157,
                               "CA" = 158,
                               "C" = 159,
                               "O" = 160,
                               "CB" = 161,
                               "CG" = 162,
                               "CD" = 163,
                               "NE" = 164,
                               "CZ" = 165,
                               "NH1" = 166,
                               "NH2" = 167,
                               0),
               0)
  v
}

# ------------------------------------------------------
# ---------------------------
# ------------- This function extract Chain ID and Start and End positions for each domain.
# -------------------------------------------------
split_position <- function(position){
  # position means the chain and residue numbers
  input_string <- position
  split_strings <- unlist(strsplit(input_string, ",")) # split the string by comma
  chains <- starts <- ends <-character(0)
  # Process each part of the string
  for (part in split_strings) {
    parts <- strsplit(part, ":")[[1]]
    chain <- parts[1]
    if (length(parts) > 1) {
      range <- unlist(strsplit(parts[2], "-"))
      start <- ifelse(length(range) == 1, NA, as.numeric(range[1]))
      end <- ifelse(length(range) == 1, NA, as.numeric(range[2]))
    } else {
      start <- NA
      end <- NA
    }
    chains <- c(chains, chain)
    starts <- c(starts, start)
    ends <- c(ends, end)
  }
  df <- data.frame(chain = chains, start = starts, end = ends)
  df[is.na(df)] <- "NA"
  return(df)
}
# ----------------------------------------------------
# -------------------------------
# ------------ This function convert PDB structure into edge list of contact graph. 
# -----------------------------------------------
NetworkFrame <- function(PDB){
  Tdist <- 5; start <- 0.75; binw <- 0.5; Tseqsep <- 1
  pdb <- PDB$atom
  resnum <- length(unique(pdb$resno))
  xyz3 <- PDB$xyz
  coordinate <- matrix(xyz3,ncol=3,byrow = TRUE)
  atnumber <- pdb$eleno
  atomname <- pdb$elety
  amacidname <- pdb$resid
  amnumber <- pdb$resno
  # ---------------- Contact extraction based on Delaunay Triangulation -------------------
  DT <- delaunayn(coordinate)
  x12 <- DT[,1:2]; x13 <- DT[,c(1,3)]; x14 <- DT[,c(1,4)]
  x23 <- DT[,2:3]; x24 <- DT[,c(2,4)]; x34 <- DT[,3:4]
  edges <-rbind(x12,x13,x14,x23,x24,x34)
  edges <- edges[!duplicated(edges),]# remove duplicated edges
  cs <- abs(amnumber[edges[,1]]-amnumber[edges[,2]])
  edges <- edges[cs>=Tseqsep,] #remove contacts within sequence seperation less than 1
  cs <- cs[cs>=Tseqsep]
  dis <- sqrt(rowSums((coordinate[edges[,1],]-coordinate[edges[,2],])^2))
  # ----------------------------------------------------------
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
# ------------------------------------------------
# -----------------------
# ---------- This function extracts the positions of secondary structures Alpha Helix and Beta Sheet. 
# --------------------------------------------------
inedx_HelixSheet <- function(PDB){
  hs <- PDB$helix$start; he <- PDB$helix$end
  ss <- PDB$sheet$start; se <- PDB$sheet$end
  indhs<-c()
  if(!is.null(hs)|length(hs)!=0){
    for (i in 1:length(hs)) {
      h <- hs[i]:he[i]
      indhs <- c(indhs,h)
    }
  }
  if(!is.null(ss)|length(ss)!=0){
    for (i in 1:length(ss)) {
      s <- ss[i]:se[i]
      indhs <- c(indhs,s)
    }
  }
  return(indhs)
}
# --------------------------------------------------------
# --------------------------------
# ------------- This function convert the contact graph into the 210D profile energy.
# ----------------------------------------------------------
Energy_SPE210 <- function(NetworkFrame,energy_dell_dunbrack){ 
  # Pairwise Contact potential
  # -------------------- Energy file Indexing ----------------------------
  cn <- 1; index167 <- matrix(0,ncol = 167,nrow = 167)
  for (i in 1:167){
    for (j in i: 167){
      index167[i,j] <- cn
      cn <- cn+1
    }
  }
  cn <- 1; index20 <- matrix(0,ncol = 20,nrow = 20)
  for (i in 1:20){
    for (j in i:20){
      index20[i,j] <- cn
      cn <- cn+1
    }
  }
  # ------------------ atomtype2num calculation --------------------
  v <- function(Atom,AminoAcid) {  
    atomtype2num167(as.character(Atom),as.character(AminoAcid))
  }
  # ----------------------
  atomi <- NetworkFrame$atif
  amacidi <- NetworkFrame$amif
  numi <- mapply(v,atomi,amacidi)
  # ----------------------------
  atomj <- NetworkFrame$atjf
  amacidj <- NetworkFrame$amjf
  numj <- mapply(v,atomj,amacidj)
  # ---------------------------------------------------
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
  ami <- NetworkFrame$amif;amj <- NetworkFrame$amjf
  am2numi <- sapply(as.character(ami), amacid2num)
  am2numj <- sapply(as.character(amj), amacid2num)
  NetworkFrame <- data.frame(NetworkFrame,am2numi,am2numj,numi,numj)
  energy_pair_am <- aggregate(NetworkFrame$energy,by=list(NetworkFrame$am2numi,NetworkFrame$am2numj),FUN=sum, na.rm=TRUE)
  # ----------------------------------
  colnames(energy_pair_am) <- c("aminoacid1","aminoacid2",'energy')
  energy210 <- rep(0, max(index20)); ninter <- dim(energy_pair_am)[1]
  for (i in 1:ninter){
    r <-energy_pair_am$aminoacid1[i]; c <-energy_pair_am$aminoacid2[i]
    rw=min(r,c); cl=max(r,c)
    ind210 <- index20[rw,cl]
    energy210[ind210] <- energy_pair_am$energy[i]
  }
  return(energy210)
}
# --------------------------------------------------
# --------------------------------
# ------------- This function calculate the 210D profile energy from Amino Acid sequence.
# ---------------------------------------------------------
Energy_CPE210 <- function(sequenceAA){
  s <- toupper(unlist(strsplit(sequenceAA,"")))
  freq <- table(s)[letters_list]
  indxNA <- which(is.na(freq))
  names(freq)[indxNA] <- letters_list[indxNA]
  freq[indxNA] <- 0
  freq_matrix <- matrix(rep(freq, each = 20), nrow = 20)
  en_fr <- aaenergy * freq_matrix
  en_fr <- en_fr/length(s)
  pair_es <- diag(freq) %*% as.matrix(en_fr)
  aa210_p <- pair_es[lower.tri(pair_es,diag = T)]
  return(aa210_p)
}
# ---------------------------------------------------------------
# ------------------------------------
# ------------ This function calculates the distance between two profiles of energy.
# -----------------------------------------------------------
manhat <- function(X, Y){
  x<-X; y<-Y
  xy <- abs(x - y)
  xy <- xy[xy != 0]
  if(length(xy)==0){
    d <- 0
  }else{
    d <- sum(xy)/length(xy)
  }
  return(d)
}
# ---------------------------------------------------
# -----------------------------------
# ---------- Customized theme for figures.
# ---------------------------------------------------------
mytheme <- function(){
  # ------------  theme ----------------
  mythem <- theme_bw(base_line_size = .2)+
    theme(axis.text.x = element_text(size=13, face = "bold"),
          axis.text.y = element_text(size=13, face = "bold"),
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 15, colour = "black"),
          legend.title = element_blank(),
          title = element_text(size = 15)
    )
  return(mythem)
}
# --------------------------------
