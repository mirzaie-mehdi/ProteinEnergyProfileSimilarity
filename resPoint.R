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
