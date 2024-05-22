library(bio3d)

rm(list = setdiff(ls(), lsf.str()))
#************Initialization***************************
inda <- function(Phi, Psi, Binw){
  indphi <- floor(Phi/Binw)+1
  indpsi <- floor(Psi/Binw)+1
  cbind(indphi,indpsi)
}
path <- '~/Desktop/proj/Dunbrack/'
binw <- 5;
filaname_chain <- read.delim(paste0(path,"list_with_chainID.txt"), header = F)
Nprotein <- length(filaname_chain[,1])

freq <- array(0,c(72,72,20))
for (NP in 1:Nprotein) {
  #NP <- 4
  print(NP)
  pdbname <- paste0(tolower(substr(filaname_chain[NP,1],1,4)),'.pdb')
  #********read pdb and select atom section, remove hydrogen atoms********
  pdb0 <- read.pdb(paste0(path,'PDBs/',pdbname))
  if(nchar(filaname_chain[NP,1])>4 & substr(filaname_chain[NP,1],5,5) != " " & !is.na(match(substr(filaname_chain[NP,1],5,5),unique(pdb0$atom$chain)))){
    ch <- substr(filaname_chain[NP,1],5,5)
    sele1 <- atom.select(pdb0,"protein",type = "ATOM", chain=ch)
  } else if(!is.na(pdb0$atom$chain[1])){
    for (ich in 1:length(unique(pdb0$atom$chain))) {
      if(!is.na(atom.select(pdb0,"protein",type = "ATOM", chain=pdb0$atom$chain[ich])$atom[ich])){
        ch <- pdb0$atom$chain[ich]
        sele1 <- atom.select(pdb0,"protein",type = "ATOM", chain=pdb0$atom$chain[ich])
        break
      }
    }
  } else if(is.na(pdb0$atom$chain[1])){
    ch <- ""
    sele1 <- atom.select(pdb0,"protein",type = "ATOM")
  } else if(is.na(atom.select(pdb0,"protein",type = "ATOM", chain=pdb0$atom$chain[ich])$atom[1])){
    ch <- pdb0$atom$chain[1]
    sele1 <- atom.select(pdb0,"protein",type = "ATOM")
  }
  sele2 <- atom.select(pdb0, "protein", elety=c("H","OXT"),inverse=TRUE)
  sele3 <- atom.select(pdb0, "noh")
  sele <- combine.select(sele1, sele2,sele3, operator="AND")
  pdb1 <- trim.pdb(pdb0,sele)
  pdb <- pdb1$atom
  pdb<-pdb[pdb$elety=='CA',]
  tor <- torsion.pdb(pdb1)
  tortbl <- tor$tbl
  annot <- strsplit(rownames(tortbl), split = "\\.")
  tortbl<-data.frame(tortbl,seqres=unlist(lapply(annot, function(x){return(x[3])})))
  df <- tortbl[,c(1,2,8)]
  ind <- which(!is.na(df$phi) & !is.na(df$psi))
  df <- df[ind,]
  df[df[,1]<0,1] <- df[df[,1]<0,1]+360
  df[df[,2]<0,2] <- df[df[,2]<0,2]+360
  am2num <- mapply(amacid2num, as.character(df$seqres))
  df <- data.frame(df, am2num)
  indarray <- inda(df$phi,df$psi,binw)
  df <- data.frame(df, indarray)
  for (i in 1:dim(df)[1]) {
    freq[df$indphi[i], df$indpsi[i], df$am2num[i]] <- freq[df$indphi[i], df$indpsi[i], df$am2num[i]]+1
  }
}
saveRDS(freq,'Data/Dunbrack/freq_torsion_Dunbrack.rds')
# **************  
freq_torsion_Dunbrack<-readRDS('Data/Dunbrack/freq_torsion_Dunbrack.rds')
freq <- freq_torsion_Dunbrack+1
RT <- 0.582
Sigma <- 0.01
E1 <- array(0, c(72,72,20))
E2 <- array(0, c(72,72,20))
fi <- array(0, c(72,72,20))
Mi <- rep(0,20)
for (i in 1:20) {
  Mi[i] <- sum(freq[,,i])
  fi[,,i] <- freq[,,i]/Mi[i]
}

fx <- rowSums(freq, dims = 2)/sum(Mi)
x<-which(fx==0,arr.ind = T)
for (i in 1:20) {
  E1[,,i] <- RT*log(1+Mi[i]*Sigma)-RT*log(1+Mi[i]*Sigma*(fi[,,i]/fx))
  E2[,,i] <- -RT*log(1+(freq[,,i]/mean(freq[,,i])))
}
for (i in 1:20) {
  print(i)
  for (j in 1:nrow(x)) {
    ro <- x[j,1]
    cl <- x[j,2]
    E1[ro, cl, i] <- 0
  }
}

saveRDS(E1,"Data/Dunbrack/energy_torsion_Dunbrack1.rds")
saveRDS(E2,"Data/Dunbrack/energy_torsion_Dunbrack2.rds")

