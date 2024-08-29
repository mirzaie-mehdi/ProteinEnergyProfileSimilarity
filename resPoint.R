df <- readRDS("Data/rds/E210.seq_spike.rds")
rownames(df) <- df$pdbID

# -------------------------------------------
clustalo_I20240806_105328_0545_97518050_p1m <- read_table("~/Downloads/clustalo-I20240806-105328-0545-97518050-p1m.pim",
                                                          col_names = FALSE, comment = "#")
dIdeseq <- data.frame(clustalo_I20240806_105328_0545_97518050_p1m[,-c(1,2)])
colnames(dIdeseq) <- rownames(dIdeseq) <- clustalo_I20240806_105328_0545_97518050_p1m$X2
dIdeseq <- dIdeseq/100
dIdeseq <- 1 - dIdeseq
dIdeseq <- dIdeseq[df$pdbID, df$pdbID]
# --------------------------
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
# -------------------

#x <- df[,-c(1:5)]/df$length
#x  <- as.matrix(x) 
#dis <- as.matrix(x %>% proxy::dist(method = manhat))
#rms <- data.frame(readRDS('Data/covidPDB/RMSD_spike.rds'))
#diag(rms)<-0
#dis <- as.matrix(as.dist(rms))

tm_vec <- read.csv('Data/covidPDB/disTM_vec_spike.csv')
colnames(tm_vec) <- rownames(tm_vec) <- df$pdbID
tm_vec <- 1 - tm_vec
dis <- as.matrix(as.dist(tm_vec))
dis <- as.matrix(as.dist(dIdeseq))

data <- df[, -c(1,2,4,5)]
colnames(data)[1]<-'group'
data$group<-factor(data$group)

# Function to calculate minimum pairwise distance between two groups
calculate_min_distance <- function(group1_indices, group2_indices, distance_matrix) {
  # Extract the submatrix for group1 and group2
  submatrix <- distance_matrix[group1_indices, group2_indices]
  return(min(submatrix))
}
# Extract groups
group_indices <- split(seq_len(nrow(data)), data$group)
# Calculate minimum distances between each pair of groups
min_distances <- matrix(NA, nrow = 3, ncol = 3, 
                        dimnames = list(levels(data$group), 
                                        levels(data$group)))
for (i in 1:3) {
  for (j in 1:3) {
    min_distances[i, j] <- calculate_min_distance(group_indices[[i]], 
                                                  group_indices[[j]],
                                                  dis)
  }
}
final_dis <- min_distances
final_dis <- (final_dis - min(final_dis))/(max(final_dis) - min(final_dis))
# Display the minimum distance matrix
print(round(final_dis,2))
# --------------------------------------------
library(cluster)

# Function to calculate medoid for a group
calculate_medoid <- function(dist_matrix, group_indices) {
  submatrix <- dist_matrix[group_indices, group_indices]
  medoid_index <- pam(submatrix, 1)$id.med
  return(grep(rownames(submatrix)[medoid_index], rownames(dist_matrix)))
}

distances_medoid <- matrix(NA, nrow = 3, ncol = 3, 
                           dimnames = list(levels(data$group), 
                                           levels(data$group)))
for (i in 1:3) {
  for (j in 1:3) {
    x <- calculate_medoid(dis, group_indices[[i]])
    y <- calculate_medoid(dis, group_indices[[j]])
    distances_medoid[i, j] <- dis[x,y]
  }
}

rownames(distances_medoid) <- colnames(distances_medoid) <- levels(data$group)

final_dis <- distances_medoid
final_dis <- (final_dis - min(final_dis))/(max(final_dis) - min(final_dis))
# Display the distance matrix
print(round(final_dis,2))
# ------------------------------------
my_sab <- matrix(0,3,3)
for (i in 1:3) {
  print(i)
  for (j in 1:3) {
    if(i>=j){
      disi <- dis[group_indices[[i]], group_indices[[i]]]
      dii <- mean(diag(as.matrix(disi[apply(disi,2,which.min),])))
      disj <- dis[group_indices[[j]], group_indices[[j]]]
      djj <- mean(diag(as.matrix(disj[apply(disj,2,which.min),])))
      disij <- dis[group_indices[[i]], group_indices[[j]]]
      disji <- dis[group_indices[[j]], group_indices[[i]]]
      dij <- mean(c(diag(as.matrix(disij[apply(disij,2,which.min),])),
                    diag(as.matrix(disji[apply(disji,2,which.min),]))))
      my_sab[i,j] <- dij - ((dii + djj)/2)
      my_sab[j,i] <- my_sab[i,j]
    }
  }
}
colnames(my_sab) <- rownames(my_sab) <- names(group_indices)

final_dis <- my_sab
final_dis <- (final_dis - min(final_dis))/(max(final_dis) - min(final_dis))
# Display the distance matrix
print(round(final_dis,2))

# -------------------------------
gt <- read.csv('~/Downloads/gt_training.CNN/train_df.csv')
gt$seq <- gt$rawseq
gt$seq <- gsub('-','', gt$seq)
gt$seq <- toupper(gt$seq)
gt$len <- nchar(gt$seq)
gt <- gt[gt$len>100,]

E210.seq_GT <- data.frame(gt[,c(2,3,10,11)],matrix(NA,nrow(gt),210))
colnames(E210.seq_GT)[-c(1:4)] <- pair_inter

hist(E210.seq_GT$len, 50)
tictoc::tic()
for (NP in 1:nrow(E210.seq_GT)) {
  #if(NP%%10==0)
  print(NP)
  seq <- unlist(strsplit(unlist(E210.seq_GT$seq[NP]),''))
  freq = c()
  for (j in 1:20){
    freq[j] <-  length(grep(tolower(letters_list[j]),tolower(seq)))
  }
  freq_matrix <- matrix(rep(freq, each = 20), nrow = 20)
  en_fr <- aaenergy * freq_matrix
  en_fr <- en_fr/length(seq)
  pair_es <- diag(freq) %*% as.matrix(en_fr)
  aa210_p<-pair_es[lower.tri(pair_es,diag = T)]
  E210.seq_GT[NP,-c(1:4)] <- aa210_p
}
df <- E210.seq_GT
saveRDS(df, '~/Downloads/gt_training.CNN/E210.seq_GT.rds')
tictoc::toc()
df <- readRDS('~/Downloads/gt_training.CNN/E210.seq_GT.rds')
df <- df[!df$fold%in%'u',]

df <- df %>%
  group_by(fold) %>%
  slice_sample(n = 4500)

# To remove the group_by structure, use ungroup
df <- data.frame(ungroup(df))

tictoc::tic()
dis <- as.matrix(proxy::dist(df[,-c(1:4)], method = manhat))
tictoc::toc()
ump <- umap(d = dis,
            n_neighbors = 20,
            min_dist = 0.01,
            input="dist")
ump_df <- data.frame(df[,1:4],
                     UMAP1 = ump$layout[, 1],
                     UMAP2 = ump$layout[, 2])

df2 <- cbind(df[,1:4], total=rowSums(df[,-c(1:4)]))
p1 <- ggplot(ump_df, aes(UMAP1, UMAP2, color= fold)) +
  geom_point(alpha=.7)

# ------------------------------
# ---- Nice data ------
path <- '~/Desktop/new_data/'
df <- readFASTA(paste0(path,'ARF-RowData/MSA_ARF.fasta'))
df <- data.frame(pdbID=names(df), align=t(data.frame(df)))
df$seq <- toupper(gsub('-','',df$align))
df$length <- nchar(df$seq)
arf <- df
arf$family <- 'ARF'

df <- readFASTA(paste0(path,'RAB-RowData/MSA_RAB.fasta'))
df <- data.frame(pdbID=names(df), align=t(data.frame(df)))
df$seq <- toupper(gsub('-','',df$align))
df$length <- nchar(df$seq)
rab <- df
rab$family <- 'RAB'

df <- readFASTA(paste0(path,'RAS-RowData/MSA_RAS.fasta'))
df <- data.frame(pdbID=names(df), align=t(data.frame(df)))
df$seq <- toupper(gsub('-','',df$align))
df$length <- nchar(df$seq)
ras <- df
ras$family <- 'RAS'

df <- readFASTA(paste0(path,'RHO-RowData/MSA_RHO.fasta'))
df <- data.frame(pdbID=names(df), align=t(data.frame(df)))
df$seq <- toupper(gsub('-','',df$align))
df$length <- nchar(df$seq)
rho <- df
rho$family <- 'RHO'

df <- rbind(arf, rab, ras, rho)

df <- data.frame(df, matrix(NA, nrow(df), 210))
colnames(df)[-c(1:5)] <- pair_inter

hist(df$length, 50)

for (NP in 1:nrow(df)) {
  print(NP)
  seq <- unlist(strsplit(unlist(df$seq[NP]),''))
  freq = c()
  for (j in 1:20){
    freq[j] <-  length(grep(tolower(letters_list[j]),tolower(seq)))
  }
  freq_matrix <- matrix(rep(freq, each = 20), nrow = 20)
  en_fr <- aaenergy * freq_matrix
  en_fr <- en_fr/length(seq)
  pair_es <- diag(freq) %*% as.matrix(en_fr)
  aa210_p<-pair_es[lower.tri(pair_es,diag = T)]
  df[NP,-c(1:5)] <- aa210_p
}

dftot <- data.frame(family=df$family, total=rowSums(df[,-c(1:5)])/df$length)
ggplot(dftot, aes(family, total))+
  geom_boxplot()

dis <- as.matrix(proxy::dist(df[,-c(1:5)], method = manhat))
ump <- umap(d = dis,
            n_neighbors = 13,
            min_dist = 0.001,
            input="dist")
ump_df <- data.frame(df[,1:5],
                     UMAP1 = ump$layout[, 1],
                     UMAP2 = ump$layout[, 2])

ggplot(ump_df, aes(UMAP1, UMAP2, color= family)) +
  geom_point(alpha=.7)


diag(dis) <- NA
fourF <- toupper(c('arf', 'rab', 'ras', 'rho'))
# ---------- 1NN --------------------
im <- apply(dis,1, FUN = which.min)
res <- data.frame(source=df$family, target=df$family[im])
## ----------
cm <- matrix(0,4,4)
colnames(cm) <- rownames(cm) <- fourF
for (sf in fourF) {
  x <- res[res$source%in%sf,]
  z <- table(x$target)
  cm[sf, names(z)]<-z
}
ac <- round(sum(diag(cm))/sum(cm),2)*100
pr <- round(diag(cm)/colSums(cm) , 2)*100
re <- round(diag(cm)/rowSums(cm) , 2)*100
f1 <- round((2*pr*re)/(pr+re))
rbind(ac,pr,re,f1)
# ------------ SPE ------------
# ------------ ARF ------------
arf2 <- data.frame(arf,
                   matrix(NA, ncol = 210, 
                          nrow = nrow(arf)))
colnames(arf2)[-c(1:5)] <- pair_inter


for (NP in 1:nrow(arf2)) {
  print(NP)
  pdbname <- arf2$pdbID[NP]
  pdbfile <- paste0(path,'ARF-RowData/',pdbname,'.pdb')
  tryCatch({
    pdb0 <- read.pdb(pdbfile)
    sele1 <- atom.select(pdb0,"protein",type = "ATOM")
    sele2 <- atom.select(pdb0, "protein", elety=c("H","OXT"),inverse=TRUE)
    sele4 <- atom.select(pdb0, "noh")
    sele <- combine.select(sele1, sele2,sele4, operator="AND")
    new_pdb <- trim.pdb(pdb0, sele,verbose = FALSE)
    Net_Frame <- NetworkFrame(PDB = new_pdb)
    energy210 <- Energy_PC210(Net_Frame,energy_dell_dunbrack)
    arf2[NP,-c(1:5)] <- energy210
  }, error = function(e) {
    cat("Error downloading file:", e$message, "\n")
  })
}
# ------------ RAB ------------
rab2 <- data.frame(rab,
                   matrix(NA, ncol = 210, 
                          nrow = nrow(rab)))
colnames(rab2)[-c(1:5)] <- pair_inter


for (NP in 1:nrow(rab2)) {
  print(NP)
  pdbname <- rab$pdbID[NP]
  pdbfile <- paste0(path,'RAB-RowData/',pdbname,'.pdb')
  tryCatch({
    pdb0 <- read.pdb(pdbfile)
    sele1 <- atom.select(pdb0,"protein",type = "ATOM")
    sele2 <- atom.select(pdb0, "protein", elety=c("H","OXT"),inverse=TRUE)
    sele4 <- atom.select(pdb0, "noh")
    sele <- combine.select(sele1, sele2,sele4, operator="AND")
    new_pdb <- trim.pdb(pdb0, sele,verbose = FALSE)
    Net_Frame <- NetworkFrame(PDB = new_pdb)
    energy210 <- Energy_PC210(Net_Frame,energy_dell_dunbrack)
    rab2[NP,-c(1:5)] <- energy210
  }, error = function(e) {
    cat("Error downloading file:", e$message, "\n")
  })
}

# ------------ RAS ------------
ras2 <- data.frame(ras,
                   matrix(NA, ncol = 210, 
                          nrow = nrow(ras)))
colnames(ras2)[-c(1:5)] <- pair_inter


for (NP in 1:nrow(ras2)) {
  print(NP)
  pdbname <- ras2$pdbID[NP]
  pdbfile <- paste0(path,'RAS-RowData/',pdbname,'.pdb')
  tryCatch({
    pdb0 <- read.pdb(pdbfile)
    sele1 <- atom.select(pdb0,"protein",type = "ATOM")
    sele2 <- atom.select(pdb0, "protein", elety=c("H","OXT"),inverse=TRUE)
    sele4 <- atom.select(pdb0, "noh")
    sele <- combine.select(sele1, sele2,sele4, operator="AND")
    new_pdb <- trim.pdb(pdb0, sele,verbose = FALSE)
    Net_Frame <- NetworkFrame(PDB = new_pdb)
    energy210 <- Energy_PC210(Net_Frame,energy_dell_dunbrack)
    ras2[NP,-c(1:5)] <- energy210
  }, error = function(e) {
    cat("Error downloading file:", e$message, "\n")
  })
}

# ------------ RHO ------------
rho2 <- data.frame(rho,
                   matrix(NA, ncol = 210, 
                          nrow = nrow(rho)))
colnames(rho2)[-c(1:5)] <- pair_inter


for (NP in 1:nrow(rho2)) {
  print(NP)
  pdbname <- rho2$pdbID[NP]
  pdbfile <- paste0(path,'RHO-RowData/',pdbname,'.pdb')
  tryCatch({
    pdb0 <- read.pdb(pdbfile)
    sele1 <- atom.select(pdb0,"protein",type = "ATOM")
    sele2 <- atom.select(pdb0, "protein", elety=c("H","OXT"),inverse=TRUE)
    sele4 <- atom.select(pdb0, "noh")
    sele <- combine.select(sele1, sele2,sele4, operator="AND")
    new_pdb <- trim.pdb(pdb0, sele,verbose = FALSE)
    Net_Frame <- NetworkFrame(PDB = new_pdb)
    energy210 <- Energy_PC210(Net_Frame,energy_dell_dunbrack)
    rho2[NP,-c(1:5)] <- energy210
  }, error = function(e) {
    cat("Error downloading file:", e$message, "\n")
  })
}


df2 <- rbind(arf2, rab2, ras2, rho2)
dftot <- data.frame(family=df2$family, total=rowSums(df2[,-c(1:5)])/df2$length)
ggplot(dftot, aes(family, total))+
  geom_boxplot()

dis <- as.matrix(proxy::dist(df2[,-c(1:5)], method = manhat))
ump <- umap(d = dis,
            n_neighbors = 13,
            min_dist = 0.001,
            input="dist")
ump_df <- data.frame(df2[,1:5],
                     UMAP1 = ump$layout[, 1],
                     UMAP2 = ump$layout[, 2])

ggplot(ump_df, aes(UMAP1, UMAP2, color= family)) +
  geom_point(alpha=.7)


diag(dis) <- NA
fourF <- toupper(c('arf', 'rab', 'ras', 'rho'))
# ---------- 1NN --------------------
im <- apply(dis,1, FUN = which.min)
res <- data.frame(source=df2$family, target=df2$family[im])
## ----------
cm <- matrix(0,4,4)
colnames(cm) <- rownames(cm) <- fourF
for (sf in fourF) {
  x <- res[res$source%in%sf,]
  z <- table(x$target)
  cm[sf, names(z)]<-z
}
ac <- round(sum(diag(cm))/sum(cm),2)*100
pr <- round(diag(cm)/colSums(cm) , 2)*100
re <- round(diag(cm)/rowSums(cm) , 2)*100
f1 <- round((2*pr*re)/(pr+re))
rbind(ac,pr,re,f1)
# --------------------------------
E210.seq <- function(seq){
  seq <- unlist(strsplit(unlist(seq),''))
  freq = c()
  for (j in 1:20){
    freq[j] <-  length(grep(tolower(letters_list[j]),tolower(seq)))
  }
  freq_matrix <- matrix(rep(freq, each = 20), nrow = 20)
  en_fr <- aaenergy * freq_matrix
  en_fr <- en_fr/length(seq)
  pair_es <- diag(freq) %*% as.matrix(en_fr)
  aa210_p<-pair_es[lower.tri(pair_es,diag = T)]
  return(aa210_p)
}


path <- '~/Desktop/new_data/'
faFiles <- list.files(paste0(path,"Bacteria"),full.names = TRUE)
#fas <- purrr::map(faFiles, readFASTA)

nSeq <- unlist(lapply(fas, length))
hist(nSeq,50)
a <- unlist(lapply(fas, FUN = function(x){unlist(lapply(x, nchar))}))
nID <- length(faFiles)

df <- data.frame(annotID=faFiles, nProt=NA,
                 matrix(NA, nrow = nID, ncol = 210))
colnames(df)[-c(1:2)] <- pair_inter
tictoc::tic()
for (NP in 1:nID) {
  #x <- fas[[NP]]
  x <- readFASTA(faFiles[NP])
  print(paste0(NP, ': ', length(x)))
  df$nProt[NP] <- length(x)
  #tictoc::tic()
  x <- lapply(x,E210.seq)
  #tictoc::toc()
  x <- do.call(rbind, x)
  aa210_p <- colMeans(x)
  rm(x)
  df[NP,-c(1:2)] <- aa210_p
}
tictoc::toc()

saveRDS(df, paste0(path,'Bacteria.rds'))
# ----------------------------
source('Functions.R')
library(protr)
astral <- readFASTA("~/Desktop/cath_scope/cath-domain-seqs-S35.fa")
length(astral)

df <- matrix(NA, nrow = length(astral), ncol = 210)
for (NP in 1:nrow(df)) {
  print(NP)
  seq <- unlist(strsplit(astral[[NP]],''))
  freq <- c()
  for (j in 1:20){
    freq[j] <-  length(grep(tolower(letters_list[j]),tolower(seq)))
  }
  freq_matrix <- matrix(rep(freq, each = 20), nrow = 20)
  en_fr <- aaenergy * freq_matrix
  en_fr <- en_fr/length(seq)
  pair_es <- diag(freq) %*% as.matrix(en_fr)
  aa210_p <- pair_es[lower.tri(pair_es,diag = T)]
  df[NP,] <- aa210_p
}

annot <- read_table("~/Desktop/cath_scope/cath-domain-list-S35.txt",col_names = F)
annot <- annot[,1:5]
colnames(annot) <- c('cathID', 'Class', 'Architecture', 'Topology', 'Superfamily')



df <- astral[!is.na(astral$FF),]
idx_NA <- which(is.na(astral$FF))
library(data.table)
tm <- data.frame(fread('../disTM_vec_astral95.csv',header = T))
tm <- tm[-idx_NA, -idx_NA]
tm <- 1 - tm
diag(tm) <- NA

sf <- table(df$family)
sf <- sf[sf>100]
#sf <- sf[names(sf)!='b.1.2.1']
idx_sf <- which(df$family%in%names(sf))
df <- df[df$family%in%names(sf),]

dis <- as.matrix(proxy::dist(df[,-c(1:8)], method = manhat))
diag(dis) <- NA


tm <- tm[idx_sf, idx_sf]

im <- apply(dis, 1, FUN = which.min)
CPE <- data.frame(df[,5:8],df[im,5:8])
CPE$TP_class       <- ifelse(CPE$class==CPE$class.1,1,0)
CPE$TP_fold        <- ifelse(CPE$fold==CPE$fold.1,1,0)
CPE$TP_superfamily <- ifelse(CPE$superfamily==CPE$superfamily.1,1,0)
CPE$TP_family      <- ifelse(CPE$family==CPE$family.1,1,0)
# Accuracy
round(colMeans(CPE[,9:12]),2)


im <- apply(tm, 1, FUN = which.min)
TM <- data.frame(df[,5:8],df[im,5:8])
TM$TP_class       <- ifelse(TM$class==TM$class.1,1,0)
TM$TP_fold        <- ifelse(TM$fold==TM$fold.1,1,0)
TM$TP_superfamily <- ifelse(TM$superfamily==TM$superfamily.1,1,0)
TM$TP_family      <- ifelse(TM$family==TM$family.1,1,0)
# Accuracy
round(colMeans(TM[,9:12]),2)

res_CPE <- aggregate(TP_family~family,CPE, mean, rm.na=T)
colnames(res_CPE)[2] <- 'CPE'
res_TM <- aggregate(TP_family~family,TM, mean, rm.na=T)
colnames(res_TM)[2] <- 'TH'
res <- left_join(res_CPE, res_TM, 'family')
sf2 <- res$family[res$CPE>.7]

# --------------------------------------
# superfamily greater than 300
scope <- readFASTA('~/Desktop/cath_scope/astral-scopedom-seqres-all-2.08-stable.fa')
annot <- read_table('~/Desktop/cath_scope/astral-scopedom-seqres-all-2.08-stable.fa', col_names = F)
annot <- annot[grep("^>", annot$X1),1:3]
colnames(annot) <- c("scopeID","ID_scope","Position")
annot$scopeID <- gsub('>','',annot$scopeID)
annot$pdbID <- substr(annot$scopeID, start = 2, stop = 5)
#----------Split scopeID into class, fold, super and family---------
split_values <- strsplit(annot$ID_scope, "\\.")
matrix_values <- matrix("", nrow = nrow(annot), ncol = 4)
for (i in 1:4) {
  matrix_values[, i] <- sapply(split_values, function(x) paste(x[1:i], collapse = "."))
}
df_split <- as.data.frame(matrix_values)
colnames(df_split) <- c("class", "fold", "superfamily", "family")
annot$Position <- gsub("[()]", "", annot$Position)
annot$position <- substr(annot$Position, start = regexpr(":", annot$Position) + 1, stop = nchar(annot$Position))
annot$chainID <- substr(annot$Position,start = 1, stop = 1)
annot$start <- NA
annot$end <- NA
split_values <- strsplit(annot$position, "-")
for (i in seq_along(split_values)) {
  print(i)
  if (length(split_values[[i]]) == 2) {
    annot$start[i] <- as.numeric(split_values[[i]][1])
    annot$end[i] <- as.numeric(split_values[[i]][2])
  }
}
annot <- cbind(annot,df_split)
colnames(annot)
annot <- annot[,-c(2,3,5)]
# --------------------
# sf100

source('Functions.R')
astral <- readRDS("Data/rds/E210.seq_astral_95.rds")
df <- astral[!is.na(astral$FF),]
idx_NA <- which(is.na(astral$FF))
library(data.table)
tm <- data.frame(fread('../disTM_vec_astral95.csv',header = T))
tm <- tm[-idx_NA, -idx_NA]
tm <- 1 - tm
diag(tm) <- NA

sf <- table(df$family)
sf <- sf[sf>100]
#sf <- sf[names(sf)!='b.1.2.1']
idx_sf <- which(df$family%in%names(sf))
df <- df[df$family%in%names(sf),]

dis <- as.matrix(proxy::dist(df[,-c(1:8)], method = manhat))
diag(dis) <- NA


tm <- tm[idx_sf, idx_sf]

im <- apply(dis, 1, FUN = which.min)
CPE <- data.frame(df[,5:8],df[im,5:8])
CPE$TP_class       <- ifelse(CPE$class==CPE$class.1,1,0)
CPE$TP_fold        <- ifelse(CPE$fold==CPE$fold.1,1,0)
CPE$TP_superfamily <- ifelse(CPE$superfamily==CPE$superfamily.1,1,0)
CPE$TP_family      <- ifelse(CPE$family==CPE$family.1,1,0)
# Accuracy
round(colMeans(CPE[,9:12]),2)


im <- apply(tm, 1, FUN = which.min)
TM <- data.frame(df[,5:8],df[im,5:8])
TM$TP_class       <- ifelse(TM$class==TM$class.1,1,0)
TM$TP_fold        <- ifelse(TM$fold==TM$fold.1,1,0)
TM$TP_superfamily <- ifelse(TM$superfamily==TM$superfamily.1,1,0)
TM$TP_family      <- ifelse(TM$family==TM$family.1,1,0)
# Accuracy
round(colMeans(TM[,9:12]),2)

res_CPE <- aggregate(TP_family~family,CPE, mean, rm.na=T)
colnames(res_CPE)[2] <- 'CPE'
res_TM <- aggregate(TP_family~family,TM, mean, rm.na=T)
colnames(res_TM)[2] <- 'TH'
res <- left_join(res_CPE, res_TM, 'family')
sf2 <- res$family[res$CPE>.7]

# ----------------------------------------
ferritin <- read.csv('Data/csv/Ferritin_Like_seq.csv')
tm <- read.csv('Data/csv/disTM_vec_Ferritin_Like_seq.csv')
colnames(tm)<-rownames(tm)<-ferritin$pdbID
diag(tm) <- 1
tm <- as.dist(1 - tm)
tree2 <- neighborNet(tm)
library(ape)

# Compute pairwise distances between nodes
node_distances <- dist.nodes(nnet)
avg_distances <- apply(node_distances, 1, mean)
f1 <- c('BMA','BMB','FAA','RNR')
group1 <- df$pdbID[df$group%in%f1]
group1<-c("1lkoA", "1yuzA", "2fzfA", "3qhbA")
group1_indices <- match(group1, nnet$tip.label)

avg_dist_to_group1 <- apply(node_distances[, group1_indices], 1, mean, na.rm = TRUE)

# Find the node with the minimal average distance to both groups
pseudo_root <- which.min(avg_distances)

# Calculate distance from pseudo-root to each leaf
leaf_indices <- which(nnet$tip.label %in% nnet$tip.label)
leaf_distances <- node_distances[pseudo_root, leaf_indices]

# Combine with leaf names
leaf_distances_df <- data.frame(Leaf = nnet$tip.label, Distance_to_Root = leaf_distances)
colnames(leaf_distances_df)[1]<-'pdbID'
z<-left_join(df[,1:3], leaf_distances_df,'pdbID')
z<-aggregate(Distance_to_Root~group,z,mean)


library(vegan)
mantel_result <- mantel(dis, as.dist(tm))
print(mantel_result)
# ---------------------------------
# Atlas
E210.seq <- function(seq){
  seq <- unlist(strsplit(unlist(seq),''))
  freq = c()
  for (j in 1:20){
    freq[j] <-  length(grep(tolower(letters_list[j]),tolower(seq)))
  }
  freq_matrix <- matrix(rep(freq, each = 20), nrow = 20)
  en_fr <- aaenergy * freq_matrix
  en_fr <- en_fr/length(seq)
  pair_es <- diag(freq) %*% as.matrix(en_fr)
  aa210_p<-pair_es[lower.tri(pair_es,diag = T)]
  return(aa210_p)
}
# --------------
path <- '~/Desktop/new_data/'
faFiles <- list.files(paste0(path,"Eukaryota"),full.names = TRUE)
fas <- map(faFiles, readFASTA)

nSeq <- unlist(lapply(fas, length))
hist(nSeq,50)
nID <- length(fas)

df <- data.frame(annotID=faFiles, nProt=NA,
                 matrix(NA, nrow = nID, ncol = 210))
colnames(df)[-c(1:2)] <- pair_inter
tictoc::tic()
for (NP in 1:nrow(df)) {
  x <- fas[[NP]]
  print(paste0(NP, ': ', length(x)))
  df$nProt[NP] <- length(x)
  #tictoc::tic()
  x <- lapply(x,E210.seq)
  #tictoc::toc()
  x <- do.call(rbind, x)
  aa210_p <- colMeans(x)
  df[NP,-c(1:2)] <- aa210_p
}
tictoc::toc()

saveRDS(df, paste0(path,'Eukaryota.rds'))
#
bac <- readRDS('~/Desktop/new_data/Bacteria.rds')
euk <- readRDS('~/Desktop/new_data/Eukaryota.rds')
arc <- readRDS('~/Desktop/new_data/Archaea.rds')

df <- data.frame(group=rep(c('Archaea','Eukaryota','Bacteria'),
                           c(nrow(arc), nrow(euk), nrow(bac))),
                 rbind(arc, euk, bac))
rownames(df)<-basename(df$annotID)

df2 <- cbind(df[,1:3], total=rowSums(df[,-c(1:3)]))
ggplot(df2, aes(group, total)) +geom_boxplot()

dis <- as.matrix(df[,-c(1:3)] %>% proxy::dist(method = manhat))
diag(dis) <- NA
# ---------- 1NN --------------------
im <- apply(dis,1, FUN = which.min)
res <- data.frame(source=df$group, target=df$group[im])
## ----------
triGR <- unique(df$group)
cm <- matrix(0,3,3)
colnames(cm) <- rownames(cm) <- triGR
for (sf in triGR) {
  x <- res[res$source%in%sf,]
  z <- table(x$target)
  cm[sf,names(z)]<-z
}
round(sum(diag(cm))/sum(cm),3)*100
# ------------------
ump <- umap(d = dis,
            input="dist",
            min_dist = 0.01,
            n_neighbors = 100) # n_components=3

ump_df <- data.frame(df[,1:3],
                     UMAP1 = ump$layout[, 1],
                     UMAP2 = ump$layout[, 2])

ggplot(ump_df, aes(UMAP1, UMAP2, color= group)) +
  geom_point(alpha=.7)+
  scale_color_manual(values = c('skyblue3', 'red', 'green'))


plot_ly(x = ~ump_df$UMAP1,
        y = ~ump_df$UMAP2,
        z = ~ump_df$UMAP3,
        color = ~factor(ump_df$group),colors = c("red", "green", "blue"),
        type = 'scatter3d', mode = 'markers')


data <- df[, -c(2,3)]
library(cluster)
colnames(data)[1]<-'group'
data$group<-factor(data$group)
group_indices <- split(seq_len(nrow(data)), data$group)
# ------------------
# Function to calculate medoid for a group
calculate_medoid <- function(dist_matrix, group_indices) {
  submatrix <- dist_matrix[group_indices, group_indices]
  medoid_index <- pam(submatrix, 1)$id.med
  return(grep(rownames(submatrix)[medoid_index], rownames(dist_matrix)))
}
distances_medoid <- matrix(NA, nrow = 3, ncol = 3,
                           dimnames = list(levels(data$group),
                                           levels(data$group)))
for (i in 1:3) {
  for (j in 1:3) {
    print(paste0(i,': ',j))
    x <- calculate_medoid(dis, group_indices[[i]])
    y <- calculate_medoid(dis, group_indices[[j]])
    distances_medoid[i, j] <- dis[x,y]
  }
}

diag(distances_medoid) <- 0
rownames(distances_medoid) <- colnames(distances_medoid) <- levels(data$group)

final_dis <- distances_medoid
final_dis <- (final_dis - min(final_dis))/(max(final_dis) - min(final_dis))
round(final_dis,3)
# -------------------------------
gt <- read.csv('~/Downloads/gt_training.CNN/train_df.csv')
gt$seq <- gt$rawseq
gt$seq <- gsub('-','', gt$seq)
gt$seq <- toupper(gt$seq)
gt$len <- nchar(gt$seq)
gt <- gt[gt$len>100,]

E210.seq_GT <- data.frame(gt[,c(2,3,10,11)],matrix(NA,nrow(gt),210))
colnames(E210.seq_GT)[-c(1:4)] <- pair_inter

hist(E210.seq_GT$len, 50)
tictoc::tic()
for (NP in 1:nrow(E210.seq_GT)) {
  #if(NP%%10==0)
  print(NP)
  seq <- unlist(strsplit(unlist(E210.seq_GT$seq[NP]),''))
  freq = c()
  for (j in 1:20){
    freq[j] <-  length(grep(tolower(letters_list[j]),tolower(seq)))
  }
  freq_matrix <- matrix(rep(freq, each = 20), nrow = 20)
  en_fr <- aaenergy * freq_matrix
  en_fr <- en_fr/length(seq)
  pair_es <- diag(freq) %*% as.matrix(en_fr)
  aa210_p<-pair_es[lower.tri(pair_es,diag = T)]
  E210.seq_GT[NP,-c(1:4)] <- aa210_p
}
df <- E210.seq_GT
saveRDS(df, '~/Downloads/gt_training.CNN/E210.seq_GT.rds')
tictoc::toc()
df <- readRDS('~/Downloads/gt_training.CNN/E210.seq_GT.rds')
df <- df[!df$fold%in%'u',]

df <- df %>%
  group_by(fold) %>%
  slice_sample(n = 4500)

# To remove the group_by structure, use ungroup
df <- data.frame(ungroup(df))

tictoc::tic()
dis <- as.matrix(proxy::dist(df[,-c(1:4)], method = manhat))
tictoc::toc()
ump <- umap(d = dis,
            n_neighbors = 20,
            min_dist = 0.01,
            input="dist")
ump_df <- data.frame(df[,1:4],
                     UMAP1 = ump$layout[, 1],
                     UMAP2 = ump$layout[, 2])

df2 <- cbind(df[,1:4], total=rowSums(df[,-c(1:4)]))
p1 <- ggplot(ump_df, aes(UMAP1, UMAP2, color= fold)) +
  geom_point(alpha=.7)
# ---------------------------------------