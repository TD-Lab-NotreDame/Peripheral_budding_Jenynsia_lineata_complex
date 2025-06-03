#install.packages("ape")
library(ape)
#install.packages("phangorn")
library(phangorn)

# Read in your trees
trees <- read.nexus("data/Jenysia_RAD_SNAPP_c.trees")

#Check how many trees show onca as root, then lineataN spliting off, followed by luxata
onca_tips <- grep("onca", trees[[1]]$tip.label, value = TRUE)
lienataN_tips <- grep("AF_lineata", trees[[1]]$tip.label, value = TRUE)
luxata_tips <- grep("luxata", trees[[1]]$tip.label, value = TRUE)

# Function to check partial topology
check_partial_topology <- function(tree) {
  # Step 1: Check if onca is monophyletic and splits first
  if (!is.monophyletic(tree, onca_tips)) return(FALSE)
  
  root_node <- Ntip(tree) + 1
  root_children <- tree$edge[tree$edge[,1] == root_node, 2]
  
  # Get tips under each side of the root
  subtree1 <- Descendants(tree, root_children[1], "tips")[[1]]
  subtree2 <- Descendants(tree, root_children[2], "tips")[[1]]
  
  tips1 <- tree$tip.label[subtree1]
  tips2 <- tree$tip.label[subtree2]
  
  # One of the root splits should be just onca
  if (!(setequal(tips1, onca_tips) || setequal(tips2, onca_tips))) return(FALSE)
  
  # Step 2: Remove onca tips to analyze next split
  tree_reduced1 <- drop.tip(tree, onca_tips)
  
  # Step 3: Check if lienataN is monophyletic and splits next
  if (!is.monophyletic(tree_reduced1, lienataN_tips)) return(FALSE)
  
  root_children2 <- tree_reduced1$edge[tree_reduced1$edge[,1] == (Ntip(tree_reduced1) + 1), 2]
  subtree1b <- Descendants(tree_reduced1, root_children2[1], "tips")[[1]]
  subtree2b <- Descendants(tree_reduced1, root_children2[2], "tips")[[1]]
  
  tips1b <- tree_reduced1$tip.label[subtree1b]
  tips2b <- tree_reduced1$tip.label[subtree2b]
  
  if (!(setequal(tips1b, lienataN_tips) || setequal(tips2b, lienataN_tips))) return(FALSE)
  
  # Step 4: Remove lienataN tips to check for luxaqta
  tree_reduced2 <- drop.tip(tree_reduced1, lienataN_tips)
  
  if (!is.monophyletic(tree_reduced2, luxata_tips)) return(FALSE)
  
  # Check luxata splits next
  root_children3 <- tree_reduced2$edge[tree_reduced2$edge[,1] == (Ntip(tree_reduced2) + 1), 2]
  subtree1c <- Descendants(tree_reduced2, root_children3[1], "tips")[[1]]
  subtree2c <- Descendants(tree_reduced2, root_children3[2], "tips")[[1]]
  
  tips1c <- tree_reduced2$tip.label[subtree1c]
  tips2c <- tree_reduced2$tip.label[subtree2c]
  
  if (!(setequal(tips1c, luxata_tips) || setequal(tips2c, luxata_tips))) return(FALSE)
  
  # Passed all three steps
  return(TRUE)
}

# Apply to all trees
match_flags <- sapply(trees, check_partial_topology)

# Report results
percent <- mean(match_flags) * 100
percent

#alternative to check manually
# https://bioinformatics.stackexchange.com/questions/6954/calculate-the-percentage-of-each-unique-phylogenetic-tree-in-a-beast-output

install.packages('devtools')
library(devtools)
install_github('santiagosnchez/rBt')
library(rBt)

#import trees
trees <- read.annot.beast('data/Jenysia_RAD_SNAPP_c.trees')
#get all unique topologies
unique_topologies <- unique.multiPhylo(trees)
#
count <- function(item, list) {
  total = 0
  for (i in 1:length(list)) {
    if (all.equal.phylo(item, list[[i]], use.edge.length = FALSE)) {
      total = total + 1
    }
  }
  return(total)
}

result <- data.frame(unique_topology = rep(0, length(unique_topologies)),
                     count = rep(0, length(unique_topologies)))
for (i in 1:length(unique_topologies)) {
  result[i, ] <- c(i, count(unique_topologies[[i]], trees))
}

result$percentage <- ((result$count/length(trees))*100)


plot(unique_topologies[[1]])


# List of row numbers that follow the target tree pattern
rows_to_sum <- c(4,5,6,18,19,21,23,47,52,53,64)

# Use slice() to select rows, then summarise to sum column 3
result %>%
  slice(rows_to_sum) %>%
  summarise(total = sum(.[[3]]))




