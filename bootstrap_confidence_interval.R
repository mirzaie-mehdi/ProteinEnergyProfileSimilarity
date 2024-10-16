# Install necessary packages if not already installed
if (!require(ape)) install.packages("ape")
if (!require(binom)) install.packages("binom")

# Load the libraries
library(ape)
library(binom)

# Load example DNA sequence data from the ape package
data(woodmouse)
# woodmouse is a DNA multiple alignment of wood mouse mitochondrial DNA sequences

# Step 1: Construct a distance matrix from the DNA sequences
d <- dist.dna(woodmouse)

# Step 2: Build a neighbor-joining tree using the distance matrix
tree <- nj(d)

# Step 3: Define a function to reconstruct the tree from bootstrapped data
my_tree_func <- function(x) {
  # Compute distance matrix from bootstrapped data
  d_boot <- dist.dna(x)
  # Reconstruct the tree using neighbor-joining method
  tree_boot <- nj(d_boot)
  return(tree_boot)
}

# Step 4: Set the number of bootstrap replicates
B <- 1000  # Increase B for more accurate confidence intervals

# Step 5: Perform bootstrap analysis
set.seed(1234)  # For reproducibility
bs_values <- boot.phylo(tree, woodmouse, my_tree_func, B=B)

# Step 6: Calculate confidence intervals for bootstrap proportions
# bs_values contains the number of times each node was recovered
# Convert counts to proportions
bootstrap_proportions <- bs_values / B

# Compute confidence intervals using the Exact (Clopper-Pearson) method
ci <- binom.confint(bs_values[-1], B, methods="exact")

# Step 7: Prepare labels with bootstrap values and confidence intervals
node_labels <- paste0(
  round(bootstrap_proportions * 100, 1), "%\n(",
  round(ci$lower * 100, 1), "-", round(ci$upper * 100, 1), "%)"
)

# Step 8: Plot the phylogenetic tree with bootstrap support and confidence intervals
plot(tree, show.node.label=FALSE, main="Phylogenetic Tree with Bootstrap Support")
nodelabels(node_labels, adj = c(-0.1, 0.5), frame = "none", cex=0.7)
# ---------------------------------------
# -- PWdisGroup
# -----------------------------
group <- factor(df$virus)
group_indices <- split(seq_len(length(group)), group)
# --------------------------------------------------
# ---------minimum pairwise distance-------------
min_pairwise_dis <- function(dis, gr, gr_idx){
  # Function to calculate minimum pairwise distance between two groups
  calculate_min_distance <- function(group1_indices,
                                     group2_indices,
                                     distance_matrix) {
    submatrix <- distance_matrix[group1_indices, group2_indices]
    return(min(submatrix))
  }
  min_dis <- matrix(NA, nrow = 3, ncol = 3,
                    dimnames = list(levels(gr), 
                                    levels(gr)))
  for (i in 1:3) {
    for (j in 1:3) {
      min_dis[i, j] <- calculate_min_distance(gr_idx[[i]],
                                              gr_idx[[j]],
                                              dis)
    }
  }
  diag(min_dis)<-0
  final_dis <- min_dis
  final_dis <- (final_dis - min(final_dis))/
    (max(final_dis) - min(final_dis))
  final_dis <- round(final_dis,2)
  return(final_dis)
}
# ----------------------------------------------------
# -----------------------Medoid---------------------
medoid_dis <- function(dis, gr, gr_idx){
  # Function to calculate medoid for a group
  calculate_medoid <- function(dist_matrix, group_indices) {
    submatrix <- dist_matrix[group_indices, group_indices]
    medoid_index <- cluster::pam(submatrix, 1)$id.med
    return(grep(rownames(submatrix)[medoid_index],
                rownames(dist_matrix)))
  }
  dis_medoid <- matrix(NA, nrow = 3, ncol = 3,
                       dimnames = list(levels(gr),
                                       levels(gr)))
  for (i in 1:3) {
    for (j in 1:3) {
      x <- calculate_medoid(dis, gr_idx[[i]])
      y <- calculate_medoid(dis, gr_idx[[j]])
      dis_medoid[i, j] <- dis[x,y]
    }
  }
  final_dis <- dis_medoid
  final_dis <- (final_dis - min(final_dis))/
    (max(final_dis) - min(final_dis))
  final_dis <- round(final_dis,2)
  return(final_dis)
}
# -----------------------------------------------
# ----------SAB like drug similarity measure---------------
sAB_dis <- function(dis, gr, gr_idx){
  my_sab <- matrix(0,3,3)
  for (i in 1:3) {
    for (j in 1:3) {
      disi <- dis[gr_idx[[i]], gr_idx[[i]]]
      dii <- mean(diag(as.matrix(disi[apply(disi,2,which.min),])))
      disj <- dis[gr_idx[[j]], gr_idx[[j]]]
      djj <- mean(diag(as.matrix(disj[apply(disj,2,which.min),])))
      disij <- dis[gr_idx[[i]], gr_idx[[j]]]
      disji <- dis[gr_idx[[j]], gr_idx[[i]]]
      dij <- mean(c(diag(as.matrix(disij[apply(disij,2,which.min),])),
                    diag(as.matrix(disji[apply(disji,2,which.min),]))))
      my_sab[i,j] <- dij - ((dii + djj)/2)
    }
  }
  final_dis <- my_sab
  final_dis <- (final_dis - min(final_dis))/
    (max(final_dis) - min(final_dis))
  final_dis <- round(final_dis,2)
  return(final_dis)
}
# -------------------
dists <- list(CPE=dis_CPE,
              SPE=dis_SPE,
              seqid=dis_seqid,
              TM_vec=tm,
              RMSD=rms)

pair_dis_minPair <- lapply(dists, FUN = min_pairwise_dis,
                           gr=group,gr_idx=group_indices)

pair_dis_medoid <- lapply(dists, FUN = medoid_dis,
                          gr=group,gr_idx=group_indices)

pair_dis_sAB <- lapply(dists, FUN = sAB_dis,
                       gr=group,gr_idx=group_indices)

