ferritin <- read.csv('Data/csv/Ferritin_Like_seq.csv')
df <- read.csv('Data/csv/disTM_vec_Ferritin_Like_seq.csv')
rownames(df)<-colnames(df)<-paste0(ferritin$pdbID, ferritin$chainID)
df <- 1 - df
hc <- as.dist(df)

df <- E210.seq_ferritin
rownames(df)<-df$pdbID

x <- df[,-c(1:4)]
x  <- as.matrix(x) 
dis <- x %>% proxy::dist(method = manhat)
tri <- ape::bionjs(as.matrix(dis))
hc <- dis %>% hclust(method = "complete")
adjustedRandIndex(cutree(hc, 8), df$group)
for (i in 2:20) {
  print(adjustedRandIndex(cutree(hc, i), df$group))
}
classError(cutree(hc, 8), df$group)
dend <- hc %>% as.dendrogram

df <- df[with(df, order(df$group)), ]
row.names(df) <- df$pdbID
rn <- df$pdbID
tree_data <- as.phylo(dend)

grp <- list(Bac     = rn[1:5],
            BMA  = rn[6:7],
            BMB  = rn[8:10],
            Dps  = rn[11:25],
            FAA  = rn[26:29],
            Fer  = rn[30:40],
            RNA  = rn[41:48],
            Rubrerythrin  = rn[49:53])

myBoots <- boot.phylo(tree_data, as.matrix(dis),B = 1000,
                      function(e) bionj(proxy::dist(e,method=manhat)))

p_data <- ggtree(tri, layout = 'unrooted',#dendrogram
                 size = 2) 
p_data$data$boots<-c(rep(NA,53),myBoots)

t_CPE = groupOTU(p_data, grp, "groups") + aes(color=groups) +
  theme(legend.position="none",
        plot.title = element_text(size = 50, face="bold",hjust = 0.5,vjust = 2))+
  #geom_tiplab2(aes(label=label,angle=angle),
  #             hjust=-1,size=3.5,fontface="bold",align=F)+
  geom_tiplab(aes(label=label),
              hjust=1.3,size=3.5,fontface="bold",align=F)

# -------------------------------------------
clustalo_I20240806_105328_0545_97518050_p1m <- read_table("~/Downloads/clustalo-I20240806-105328-0545-97518050-p1m.pim",
                                                          col_names = FALSE, comment = "#")
dIdeseq<-data.frame(clustalo_I20240806_105328_0545_97518050_p1m[,-c(1,2)])
colnames(dIdeseq)<-rownames(dIdeseq)<-clustalo_I20240806_105328_0545_97518050_p1m$X2
dIdeseq<-dIdeseq/100
dIdeseq<- 1 - dIdeseq
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
df <- readRDS("Data/rds/E210.seq_spike.rds")

rownames(df) <- df$pdbID
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
df <- E210.seq_GT[!E210.seq_GT$fold%in%'u',]

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


