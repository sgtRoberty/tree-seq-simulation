#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(FossilSim)
  library(ape)
  library(rjson)
})

# ---- Parse Named Command-Line Arguments ----
args <- commandArgs(trailingOnly = TRUE)

# Default values
SIMULATION_NUM <- 10
EXTANT_NUM <- 20 # Number of extant tips
MIN_FOSSIL_NUM <- 0 # Min number of fossil tips
MAX_FOSSIL_NUM <- Inf # Max number of fossil tips

# Helper function to extract named args
get_arg_value <- function(name, default) {
  flag <- paste0("--", name)
  if (flag %in% args) {
    idx <- match(flag, args)
    return(args[idx + 1])
  } else {
    return(default)
  }
}

# Parse named arguments
SIMULATION_NUM <- as.integer(get_arg_value("simulation-number", SIMULATION_NUM))
EXTANT_NUM <- as.integer(get_arg_value("extant-number", EXTANT_NUM))
MIN_FOSSIL_NUM <- as.integer(get_arg_value("min-fossil-number", MIN_FOSSIL_NUM))
MAX_FOSSIL_NUM <- as.numeric(get_arg_value("max-fossil-number", MAX_FOSSIL_NUM))

cat(sprintf("Using simulation-number = %s, extant-number = %s, min-fossil-number = %s, max-fossil-number = %s\n", SIMULATION_NUM, EXTANT_NUM, MIN_FOSSIL_NUM, MAX_FOSSIL_NUM))

# Initialize tree string lists
n <- SIMULATION_NUM
raw_time_trees <- character(n)
time_trees <- character(n)
distance_trees <- character(n)

build_tree <- function(tree, clockMean, clockSD) {
  node_heights <- node.depth.edgelength(tree)
  tree_height <- max(node_heights)
  origin_time <- tree_height + tree$root.edge
  parents <- tree$edge[, 1]
  children <- tree$edge[, 2]
  time_durations <- tree$edge.length
  num_edges <- nrow(tree$edge)

  branches_df <- data.frame(
    edge_index = 1:num_edges,
    parent = parents,
    child = children,
    origin_MYA = tree_height - node_heights[parents],
    end_MYA = tree_height - node_heights[children],
    time_duration = time_durations
  )

  # Add root edge if it exists
  has_root_edge <- !is.null(tree$root.edge)
  if (has_root_edge) {
    root_node <- Ntip(tree) + 1
    root_branch <- data.frame(
      edge_index = NA,
      parent = NA,
      child = root_node,
      origin_MYA = origin_time,
      end_MYA = tree_height,
      time_duration = tree$root.edge
    )
    branches_df <- rbind(branches_df, root_branch)
  }

  branches_df$rel_clock_rate <- rlnorm(nrow(branches_df), meanlog = -(clockSD)^2 / 2, sdlog = clockSD)
  branches_df$abs_clock_rate <- branches_df$rel_clock_rate * clockMean

  branches_df$branch_distance <- branches_df$time_duration * branches_df$abs_clock_rate

  build_newick <- function(node, distance_tree) {
    children <- which(parents == node)
    if (length(children) == 0) {
      tip_age <- tree_height - node_heights[node]

      if (tip_age < 1e-12) {
        label <- paste0("rhoSample_", tree$tip.label[node])
      } else {
        label <- paste0("psiSample_", tree$tip.label[node])
      }
    } else {
      subtrees <- sapply(
        tree$edge[tree$edge[, 1] == node, 2],
        function(x) build_newick(x, distance_tree)
      )
      label <- paste0("(", paste(subtrees, collapse = ","), ")")
    }
  
    # Root node number
    root_node <- Ntip(tree) + 1
  
  if (node == root_node && has_root_edge) {
    branch_data <- branches_df[nrow(branches_df), ]
  } else {
    edge_row <- match(node, tree$edge[, 2])
    branch_data <- branches_df[edge_row, ]
  }

  if (distance_tree) {
    return(sprintf("%s:%.16f", label, branch_data$branch_distance))
  } else {
    annotation <- sprintf("[&rel_clock_rate=%.16f,abs_clock_rate=%.16f]",
                          branch_data$rel_clock_rate,
                          branch_data$abs_clock_rate)
    return(sprintf("%s%s:%.16f", label, annotation, branch_data$time_duration))
  }

  }

  root_node <- Ntip(tree) + 1

  return(list(
  time_tree = paste0(build_newick(root_node, FALSE), ";"),
  distance_tree = paste0(build_newick(root_node, TRUE), ";")
  ))

}

# Initialize vector to store parameters
lambda_vec <- numeric(n)
mu_vec <- numeric(n)
psi_vec <- numeric(n)
rho_vec <- numeric(n)
clockMean_vec <- numeric(n)
clockSD_vec <- numeric(n)

num_taxa_vec <- numeric(n)
num_sa_vec <- numeric(n)
root_age_vec <- numeric(n)
origin_time_vec <- numeric(n)
oldest_sample_age_vec <- numeric(n)
tree_length_vec <- numeric(n)

# Main loop to simulate trees
for (i in 1:n) {

  lambda <- 0.1
  mu <- 0.01
  psi <- 0.01
  rho <- 1

  clockMean <- 0.01
  clockSD <- 0.5

  max_tips <- EXTANT_NUM + MAX_FOSSIL_NUM # Max number of tips
  min_tips <- EXTANT_NUM + MIN_FOSSIL_NUM # Min number of tips
  try <- 0

  cat("Simulating tree", i, "\n")
  cat(sprintf("Using lambda = %.4f, mu = %.4f, psi = %.4f, rho = %.4f\n", lambda, mu, psi, rho))

  while (TRUE) {
    try <- try + 1

    sim_tree <- sim.fbd.taxa(
      n = EXTANT_NUM,
      numbsim = 1,
      lambda = lambda,
      mu = mu,
      psi = psi,
      frac = rho,
      complete = FALSE
    )[[1]]

    ntips <- Ntip(sim_tree)

    if (ntips >= min_tips && ntips <= max_tips) {
      cat("  Accepted after", try, "tries (", ntips, "tips )\n")
      break
    } else {
      # Provide specific feedback on why it failed
      reason <- ifelse(ntips < min_tips, "too few", "too many")
      limit <- ifelse(ntips < min_tips, min_tips, max_tips)
    
      cat("  Try", try, "failed:", ntips, "tips is", reason, 
          "(limit:", limit, ") — simulating again...\n")
    }
  }

  built_trees <- build_tree(sim_tree, clockMean = clockMean, clockSD = clockSD)

  # Store parameters
  lambda_vec[i] <- lambda
  mu_vec[i] <- mu
  psi_vec[i] <- psi
  rho_vec[i] <- rho
  clockMean_vec[i] <- clockMean
  clockSD_vec[i] <- clockSD

  num_taxa_vec[i] <- ntips
  is_tip_edge <- sim_tree$edge[, 2] <= Ntip(sim_tree)
  num_sa_vec[i] <- sum(sim_tree$edge.length[is_tip_edge] == 0)
  node_heights <- node.depth.edgelength(sim_tree)
  root_age_vec[i] <- max(node_heights)
  origin_time_vec[i] <- root_age_vec[i] + sim_tree$root.edge

  tip_heights <- max(node_heights) - node_heights[1:Ntip(sim_tree)]
  tip_heights[abs(tip_heights) < 1e-12] <- 0
  oldest_sample_age_vec[i] <- max(tip_heights)

  tree_length_vec[i] <- sum(sim_tree$edge.length)

  # Store trees
  raw_time_trees[i] <- write.tree(sim_tree, digits = 16)
  time_trees[i] <- built_trees$time_tree
  distance_trees[i] <- built_trees$distance_tree

}

# Write NEXUS tree file
write_nexus_file <- function(filename, trees) {
  cat("#NEXUS\n\nBegin trees;\n", file = filename)
  
  con <- file(filename, open = "a")
  for (i in seq_along(trees)) {
    writeLines(paste0("tree STATE_", i - 1, " = ", trees[[i]]), con = con)
  }
  close(con)
  
  cat("End;\n", file = filename, append = TRUE)
}

truth_dir <- "truth"
dir.create(truth_dir, showWarnings = FALSE)

write_nexus_file(file.path(truth_dir, "truth_time_raw.trees"), raw_time_trees)
write_nexus_file(file.path(truth_dir, "truth_time.trees"), time_trees)
write_nexus_file(file.path(truth_dir, "truth_distance.trees"), distance_trees)

# Save truth.log
nTaxa <- num_taxa_vec
nSampledAncestors <- num_sa_vec
rootAge <- root_age_vec
originTime <- origin_time_vec
oldestSampleAge <- oldest_sample_age_vec
treeLength <- tree_length_vec

lambda <- lambda_vec
mu <- mu_vec
psi <- psi_vec
rho <- rho_vec

clockMean <- clockMean_vec
clockSD <- clockSD_vec

truth_df <- data.frame(nTaxa, nSampledAncestors, rootAge, originTime, oldestSampleAge, treeLength, lambda, mu, psi, rho, clockMean, clockSD)
write.table(truth_df, file = file.path(truth_dir, "truth.log"), sep = "\t", quote = FALSE, row.names = FALSE)

# Main loop to simulate sequences
for (i in seq_along(time_trees)) {

  rep_dir <- paste0("rep", i)

  dir.create(rep_dir, showWarnings = FALSE)

  time_file <- file.path(rep_dir, paste0("truth_time_rep", i, ".tre"))
  dist_file <- file.path(rep_dir, paste0("truth_distance_rep", i, ".tre"))

  writeLines(time_trees[[i]], con = time_file)
  writeLines(distance_trees[[i]], con = dist_file)

  JSON = list()
  for (x in colnames(truth_df)) {
    JSON[[x]] = truth_df[i,x]
  }
  JSON_str = as.character(rjson::toJSON(JSON, indent = 1))
  write(JSON_str, paste0(rep_dir, "/var.json"))

  old_wd <- getwd()
  
  setwd(rep_dir)

  args <- c(
    "--alisim", paste0("alignment_rep", i),
    "-t", paste0("truth_distance_rep", i, ".tre"),
    "-m", "GTR{1.0,1.0,1.0,1.0,1.0,1.0}+G4{0.5}+F{0.25,0.25,0.25,0.25}",
    "--out-format", "fasta"
  )

  system2("iqtree3", args = args, stdout = TRUE, stderr = TRUE, wait = TRUE)

  setwd(old_wd)

}