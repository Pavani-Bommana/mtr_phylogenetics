install.packages(ggtree)
# Step 1: Install BiocManager if you don’t have it
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Step 2: Install ggtree from Bioconductor
BiocManager::install("ggtree")

# Step 3: Load the package
library(ggtree)
installed.packages()["ggtree","Version"]
library(ape)

# Load your tree
tree <- read.tree("C:/Users/pvani/Downloads/aligned_filtered_unique_overall_WGS_hits_withref.fasta.treefile")

# Plot without tip labels
plot(tree, show.tip.label = FALSE)
library(ggplot2)
geom_segment(size = 0.3)   # or linewidth = 0.3
library(ggplot2)

library(ggplot2)

library(ape)

# Load your tree
tree <- read.tree("C:\\Users\\pvani\\Downloads\\promoter_regions.fasta.treefile")

# Plot without tip labels
plot(tree, show.tip.label = FALSE)


ggtree(tree, layout="circular") + geom_tiplab(size=2)
# Circular tree without tip labels
plot(tree, type = "fan", show.tip.label = FALSE)

library(ape)

tree <- read.tree("C:/Users/pvani/Downloads/aligned_filtered_unique_overall_WGS_hits_withref.fasta.treefile")

# Circular tree with thinner lines and more spacing
plot(tree,
     type = "fan",          # circular layout
     show.tip.label = FALSE, # no labels
     edge.width = 0.5,       # thinner branches
     cex = 0.3)              # smaller tip symbols
################################################################################################################

tree <- read.tree("C:/Users/pvani/Downloads/aligned_filtered_unique_overall_WGS_hits_withref.fasta.treefile")
# All branch lengths
branch_lengths <- tree$edge.length

# Quick summary
summary(branch_lengths)

# Histogram
hist(branch_lengths, breaks = 50, col = "skyblue",
     main = "Distribution of Branch Lengths",
     xlab = "Branch length (substitutions/site)")

# Parse node labels into bootstrap and SH-aLRT
parse_support <- function(lbl) {
  if (is.na(lbl) || lbl == "") return(c(NA, NA))
  parts <- strsplit(lbl, "/", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    return(as.numeric(parts))
  } else {
    return(c(as.numeric(lbl), NA))
  }
}

supports <- t(sapply(tree$node.label, parse_support))
colnames(supports) <- c("bootstrap", "shalrt")

# Histograms
hist(supports[, "bootstrap"], breaks = 30, col = "lightgreen",
     main = "Bootstrap Support Distribution", xlab = "Bootstrap %")

hist(supports[, "shalrt"], breaks = 30, col = "orange",
     main = "SH-aLRT Support Distribution", xlab = "SH-aLRT %")


library(ape)
library(ggplot2)



# Parse support values
parse_support <- function(lbl) {
  if (is.na(lbl) || lbl == "") return(c(NA, NA))
  parts <- strsplit(lbl, "/", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    return(as.numeric(parts))
  } else {
    return(c(as.numeric(lbl), NA))
  }
}

supports <- t(sapply(tree$node.label, parse_support))
colnames(supports) <- c("bootstrap", "shalrt")

# Internal nodes are numbered Ntip+1 : Ntip+Nnode
Ntip <- length(tree$tip.label)
Nnode <- tree$Nnode
internal_nodes <- Ntip + seq_len(Nnode)

# Each edge connects parent -> child
edges <- tree$edge
edge_lengths <- tree$edge.length

# Find edges whose child is an internal node
internal_edges <- which(edges[,2] %in% internal_nodes)

# Build data frame: branch length + support for that child node
df <- data.frame(
  branch_length = edge_lengths[internal_edges],
  bootstrap = supports[, "bootstrap"],
  shalrt = supports[, "shalrt"]
)

# Scatter plots
ggplot(df, aes(x = branch_length, y = bootstrap)) +
  geom_point(color = "steelblue") +
  labs(title = "Branch length vs Bootstrap support",
       x = "Branch length", y = "Bootstrap (%)")

ggplot(df, aes(x = branch_length, y = shalrt)) +
  geom_point(color = "orange") +
  labs(title = "Branch length vs SH-aLRT support",
       x = "Branch length", y = "SH-aLRT (%)")


tree <- read.tree("C:/Users/pvani/Downloads/aligned_filtered_unique_overall_WGS_hits_withref.fasta.treefile")

library(ape)
library(ggplot2)


# Parse support values
parse_support <- function(lbl) {
  if (is.na(lbl) || lbl == "") return(c(NA, NA))
  parts <- strsplit(lbl, "/", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    return(as.numeric(parts))
  } else {
    return(c(suppressWarnings(as.numeric(lbl)), NA))
  }
}

supports <- t(sapply(tree$node.label, parse_support))
colnames(supports) <- c("bootstrap", "shalrt")

# Internal nodes are numbered after the tips
Ntip <- length(tree$tip.label)
internal_nodes <- Ntip + seq_along(tree$node.label)

# For each internal node, find the edge leading into it
edges <- tree$edge
edge_lengths <- tree$edge.length

node_to_edge <- match(internal_nodes, edges[,2])

df <- data.frame(
  node = internal_nodes,
  branch_length = edge_lengths[node_to_edge],
  bootstrap = supports[, "bootstrap"],
  shalrt = supports[, "shalrt"]
)

str(df)   # should show a proper data frame with equal-length columns
ggplot(df, aes(x = branch_length, y = bootstrap)) +
  geom_point(color = "steelblue") +
  labs(title = "Branch length vs Bootstrap support",
       x = "Branch length", y = "Bootstrap (%)")

ggplot(df, aes(x = branch_length, y = shalrt)) +
  geom_point(color = "orange") +
  labs(title = "Branch length vs SH-aLRT support",
       x = "Branch length", y = "SH-aLRT (%)")
###########################################################################

library(ape)

# Extract all branch lengths
branch_lengths <- tree$edge.length

# Calculate the 75th percentile
cutoff <- quantile(branch_lengths, 0.75)

cutoff


# --- 1. Branch length cutoff (75th percentile) ---
branch_lengths <- tree$edge.length
bl_cutoff <- quantile(branch_lengths, 0.75)   # ≈ 0.0001539713

# --- 2. Parse support values (bootstrap / SH-aLRT) ---
parse_support <- function(lbl) {
  if (is.na(lbl) || lbl == "") return(c(NA, NA))
  parts <- strsplit(lbl, "/", fixed = TRUE)[[1]]
  if (length(parts) == 2) {
    return(as.numeric(parts))
  } else {
    return(c(suppressWarnings(as.numeric(lbl)), NA))
  }
}

supports <- t(sapply(tree$node.label, parse_support))
colnames(supports) <- c("bootstrap", "shalrt")

# --- 3. Internal nodes and their incoming edge lengths ---
Ntip <- length(tree$tip.label)
internal_nodes <- Ntip + seq_along(tree$node.label)
edges <- tree$edge
edge_lengths <- tree$edge.length
node_to_edge <- match(internal_nodes, edges[,2])

df <- data.frame(
  node = internal_nodes,
  branch_length = edge_lengths[node_to_edge],
  bootstrap = supports[, "bootstrap"],
  shalrt = supports[, "shalrt"]
)

# --- 4. Apply thresholds ---
boot_thresh <- 70
shalrt_thresh <- 80

keepers <- subset(df,
                  branch_length >= bl_cutoff &
                    bootstrap >= boot_thresh &
                    shalrt >= shalrt_thresh)

# --- 5. Extract clades that meet criteria ---
selected_clades <- list()
for (nd in keepers$node) {
  cl <- try(extract.clade(tree, node = nd), silent = TRUE)
  if (inherits(cl, "phylo")) {
    selected_clades[[paste0("node_", nd)]] <- cl
  }
}

# --- 6. Plot selected clades ---
for (nm in names(selected_clades)) {
  plot(selected_clades[[nm]], type = "fan", show.tip.label = FALSE,
       edge.width = 0.5, main = nm)
}
par(mar = c(1, 1, 1, 1))   # shrink margins
for (nm in names(selected_clades)) {
  dev.new(width = 8, height = 8)   # open a new graphics window
  plot(selected_clades[[nm]], type = "fan", show.tip.label = FALSE,
       edge.width = 0.5, main = nm)
}
for (nm in names(selected_clades)) {
  png(paste0(nm, ".png"), width = 1600, height = 1600, res = 300)
  plot(selected_clades[[nm]], type = "fan", show.tip.label = FALSE,
       edge.width = 0.5, main = nm)
  dev.off()
}
par(mfrow = c(1,1))

################################################################################
library(ggtree)
library(treeio)
library(dplyr)

# ---- Load IQ-TREE tree file ----
tree <- read.tree("C:/Users/pvani/Downloads/aligned_filtered_unique_overall_WGS_hits_withref.fasta.treefile")
class(tree)
# ====== Clean IQ-TREE phylo tree by length + support and plot ======
# Works with plain "phylo" objects (tree$edge.length, tree$node.label)

library(ape)       # read.tree, di2multi
library(ggtree)    # plotting
library(ggplot2)   # for ggsave if you want
library(stringr)   # parsing support labels

# -------- User inputs (edit as needed) ----------
treefile <- "C:/Users/pvani/Downloads/aligned_filtered_unique_overall_WGS_hits_withref.fasta.treefile"
# path to your .treefile
boot_thresh   <- 70                 # UFBoot threshold (user-specified)
shalrt_thresh <- 80                 # SH-aLRT threshold (user-specified)
pctile        <- 0.75               # percentile for branch-length cutoff (75th)
out_png       <- "clean_tree.png"   # output image file
# ------------------------------------------------

# 1) Load tree (phylo)
tree <- read.tree(treefile)

# Basic checks
if (!inherits(tree, "phylo")) stop("Loaded object is not a phylo tree.")

Ntip <- length(tree$tip.label)
Nnode <- tree$Nnode   # number of internal nodes

# 2) Branch-length cutoff (75th percentile of all edge.lengths)
bl_all <- tree$edge.length
bl_cutoff <- as.numeric(quantile(bl_all, probs = pctile, na.rm = TRUE))
message("Branch-length 75th percentile cutoff = ", signif(bl_cutoff, 6))

# 3) Parse node support labels
# node labels in a phylo object correspond to internal nodes:
# tree$node.label[1] -> node number Ntip + 1, etc.
node_labels <- tree$node.label
if (is.null(node_labels)) {
  # no internal node labels present
  node_labels <- rep(NA, Nnode)
}

# Prepare vectors for parsed support values (length = Nnode)
SH_vals <- rep(NA_real_, Nnode)
UF_vals <- rep(NA_real_, Nnode)
library(stringr)
for (i in seq_len(Nnode)) {
  lbl <- as.character(node_labels[i])
  if (is.na(lbl) || lbl == "" || lbl == "NA") next
  
  # Common IQ-TREE formats:
  #  "SH/UFBoot"  e.g. "95/85"
  #  "95" (single value; could be bootstrap only)
  #  "95/NA" or "NA/85"
  if (str_detect(lbl, "/")) {
    parts <- str_split_fixed(lbl, "/", 2)
    # try to convert; if fails -> NA
    SH_vals[i] <- suppressWarnings(as.numeric(parts[1]))
    UF_vals[i] <- suppressWarnings(as.numeric(parts[2]))
  } else {
    # single value present: we'll conservatively treat it as UFBoot if typical,
    # but we can't be sure — so store it into UF_vals and leave SH NA.
    tmp <- suppressWarnings(as.numeric(lbl))
    UF_vals[i] <- tmp
  }
}


# 4) Determine which internal nodes PASS support thresholds
SH_pass <- ifelse(!is.na(SH_vals) & SH_vals >= shalrt_thresh, TRUE, FALSE)
UF_pass <- ifelse(!is.na(UF_vals) & UF_vals >= boot_thresh, TRUE, FALSE)

# Define "both_pass" meaning the node passed both tests.
# If a support measure is missing (NA), we treat that measure as failing.
both_pass <- (SH_pass & UF_pass)

# 5) For each edge, decide whether to collapse it:
#    - collapse edges with length < bl_cutoff
#    - collapse edges whose child node is an internal node that failed support
#
# In a phylo object, tree$edge is a matrix nrow = number of edges, columns:
#   parent (col 1), child (col 2)
# Internal node numbers: (Ntip + 1) .. (Ntip + Nnode)
edges <- tree$edge
edge_len <- tree$edge.length

# Create flag vector: TRUE => collapse (set length to 0)
collapse_flag <- rep(FALSE, length(edge_len))

for (ei in seq_along(edge_len)) {
  # 1) short branch
  if (!is.na(edge_len[ei]) && edge_len[ei] < bl_cutoff) {
    collapse_flag[ei] <- TRUE
    next
  }
  # 2) support fail at child internal node
  child <- edges[ei, 2]
  if (child > Ntip) {
    # index into node_labels is (child - Ntip)
    nidx <- child - Ntip
    # if this internal node did NOT pass both supports -> collapse
    if (!isTRUE(both_pass[nidx])) {
      collapse_flag[ei] <- TRUE
    }
  }
}

# Report summary
message("Edges total: ", length(edge_len))
message("Edges flagged for collapse (either short or low support): ", sum(collapse_flag))

# 6) Set flagged edge lengths to zero (so they can be collapsed)
tree2 <- tree
tree2$edge.length[which(collapse_flag)] <- 0

# 7) Collapse zero-length branches (convert to multifurcations)
# Use di2multi with small tol (zero lengths <= tol are collapsed). Use a tiny tol.
tree_collapsed <- di2multi(tree2, tol = 1e-12)

# 8) Update node labels: keep labels only for nodes that passed (and still exist)
# After collapsing, node numbering may change. A safer approach is to re-calc
# supports for nodes by mapping original internal nodes -> their MRCA in the new tree
# and copy label if the original node passed both thresholds.
#
# We'll try to map by tip sets: for each original internal node, get its descendant tips,
# find MRCA in the collapsed tree, and if that MRCA exists, set its label to the original label
# only if the original node passed.
orig_tree <- tree  # keep original
orig_node_desc_tips <- vector("list", Nnode)
install.packages("phangorn")   # only once
library(phangorn)

for (i in seq_len(Nnode)) {
  node_num <- Ntip + i
  desc_tips <- Descendants(orig_tree, node = node_num, type = "tips")[[1]]
  orig_node_desc_tips[[i]] <- sort(orig_tree$tip.label[desc_tips])
}

# Initialize collapsed labels as NA
new_node_labels <- rep(NA_character_, tree_collapsed$Nnode)

# For each original internal node that passed both supports, try to map it
for (i in seq_len(Nnode)) {
  if (!both_pass[i]) next
  tipset <- orig_node_desc_tips[[i]]
  if (length(tipset) == 0) next
  # find MRCA in collapsed tree
  mrca_node <- tryCatch({
    getMRCA(tree_collapsed, tipset)
  }, error = function(e) NA)
  if (is.na(mrca_node)) next
  # collapsed tree internal nodes are numbered Ntip_coll + 1 : ...
  idx_coll <- mrca_node - length(tree_collapsed$tip.label)
  if (idx_coll >= 1 && idx_coll <= tree_collapsed$Nnode) {
    # copy original label (if present) to collapsed node label
    lbl <- orig_tree$node.label[i]
    if (!is.null(lbl) && lbl != "" && !is.na(lbl)) {
      new_node_labels[idx_coll] <- lbl
    }
  }
}

# Assign reconstructed labels to collapsed tree
tree_collapsed$node.label <- new_node_labels

# 9) Plot with ggtree: show only nodes that have labels (i.e., those that passed)
p <- ggtree(tree_collapsed, layout = "rectangular", right = FALSE) +
  geom_tiplab(size = 3) +
  geom_text2(aes(subset = !is.na(label) & label != "", label = label),
             hjust = -0.3, size = 3) +
  theme_tree2()

# Save to file
ggsave(out_png, p, width = 10, height = 10, units = "in", dpi = 300)
message("Saved cleaned tree to: ", out_png)

# Also print the plot to R's current device
print(p)
# ====== End script ======
library(ape)
library(ggtree)
library(stringr)   # for str_detect, str_split
library(dplyr)

# -------------------- USER SETTINGS --------------------
treefile <- "C:/Users/pvani/Downloads/aligned_filtered_unique_overall_WGS_hits_withref.fasta.treefile"
boot_thresh   <- 70
shalrt_thresh <- 80
pctile        <- 0.75
# -------------------------------------------------------

# 1. Load tree
tree <- read.tree(treefile)

# 2. Extract branch lengths
bl <- tree$edge.length
bl_cutoff <- quantile(bl, pctile, na.rm = TRUE)

# 3. Extract internal node labels (IQ-TREE format: "SH/UFBoot")
labs <- tree$node.label

# Prepare storage
SH_vals  <- rep(NA, length(labs))
UF_vals  <- rep(NA, length(labs))

for (i in seq_along(labs)) {
  lbl <- labs[i]
  if (!is.na(lbl) && str_detect(lbl, "/")) {
    tmp <- str_split(lbl, "/")[[1]]
    SH_vals[i] <- as.numeric(tmp[1])
    UF_vals[i] <- as.numeric(tmp[2])
  }
}

# 4. Apply thresholds on support + branch length
keep_node <- rep(FALSE, length(labs))

for (i in seq_along(labs)) {
  node_num <- i + length(tree$tip.label)  # internal node number
  edge_index <- which(tree$edge[,2] == node_num)
  
  if (!is.na(SH_vals[i]) &&
      !is.na(UF_vals[i]) &&
      SH_vals[i] >= shalrt_thresh &&
      UF_vals[i] >= boot_thresh &&
      bl[edge_index] >= bl_cutoff) {
    keep_node[i] <- TRUE
  }
}

# New node labels: keep only passing nodes
clean_labels <- ifelse(keep_node, tree$node.label, NA)
tree_clean <- tree
tree_clean$node.label <- clean_labels



p_nolabels <- ggtree(tree_clean)
print(p_nolabels)
ggsave("clean_tree_no_labels.png", p_nolabels, width=10, height=10, dpi=300)


library(ape)
library(ggtree)

# -------------------- USER SETTINGS --------------------
treefile <- "C:/Users/pvani/Downloads/aligned_filtered_unique_overall_WGS_hits_withref.fasta.treefile"
pctile   <- 0.75   # percentile cutoff for branch length
out_png  <- "clean_tree_length_filtered.png"
# -------------------------------------------------------

# 1. Load tree
tree <- read.tree(treefile)

# 2. Compute branch-length cutoff
bl <- tree$edge.length
bl_cutoff <- quantile(bl, pctile, na.rm = TRUE)
cat("Branch-length cutoff (", pctile*100, "th percentile) =", bl_cutoff, "\n")

# 3. Collapse branches shorter than cutoff
tree2 <- tree
tree2$edge.length[tree2$edge.length < bl_cutoff] <- 0
tree_clean <- di2multi(tree2, tol = 1e-12)

# 4. Plot tree WITHOUT labels first
p_nolabels <- ggtree(tree_clean)
print(p_nolabels)

# Save
ggsave(out_png, p_nolabels, width=12, height=12, dpi=300)
cat("Saved length-filtered tree as:", out_png, "\n")


library(ggtree)

# Plot circular tree WITHOUT labels
p_circular <- ggtree(tree_clean, layout = "circular")
print(p_circular)

# Optional: save the figure
ggsave("tree_circular_no_labels.png", p_circular, width = 12, height = 12, dpi = 300)
# tree_clean is your length-filtered tree
num_tips <- length(tree_clean$tip.label)
cat("Number of tips after branch-length filtering:", num_tips, "\n")
# Example: drop tips whose terminal branch < cutoff
short_tips <- tree_clean$tip.label[tree_clean$edge.length[1:length(tree_clean$tip.label)] < bl_cutoff]
tree_pruned <- drop.tip(tree_clean, short_tips)

# Number of tips remaining
num_tips_pruned <- length(tree_pruned$tip.label)
cat("Number of tips after dropping short terminal branches:", num_tips_pruned, "\n")

library(ape)
library(ggtree)

# -------------------- USER SETTINGS --------------------
treefile <- "C:/Users/pvani/Downloads/aligned_filtered_unique_overall_WGS_hits_withref.fasta.treefile"
pctile   <- 0.75   # branch-length percentile cutoff
out_png  <- "tree_862_rectangular.png"
# -------------------------------------------------------

# 1️⃣ Load tree
tree <- read.tree(treefile)

# 2️⃣ Compute branch-length cutoff (75th percentile)
bl <- tree$edge.length
bl_cutoff <- quantile(bl, pctile, na.rm = TRUE)
cat("Branch-length cutoff (", pctile*100, "th percentile) =", bl_cutoff, "\n")

# 3️⃣ Identify tips with terminal branch above cutoff
terminal_edges <- 1:length(tree$tip.label)
tips_keep <- tree$tip.label[tree$edge.length[terminal_edges] >= bl_cutoff]
cat("Number of tips above cutoff:", length(tips_keep), "\n")  # Should print 862

# 4️⃣ Prune tree to keep only those tips
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, tips_keep))

# 5️⃣ Collapse remaining short internal branches
tree_pruned$edge.length[tree_pruned$edge.length < bl_cutoff] <- 0
tree_clean <- di2multi(tree_pruned, tol = 1e-12)

# 6️⃣ Plot rectangular tree WITHOUT labels
p <- ggtree(tree_clean)
print(p)

# 7️⃣ Save figure
ggsave(out_png, p, width=12, height=12, dpi=300)
cat("Saved rectangular tree with 862 tips as:", out_png, "\n")

