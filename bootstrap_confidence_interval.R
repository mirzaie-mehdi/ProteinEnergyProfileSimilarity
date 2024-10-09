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
