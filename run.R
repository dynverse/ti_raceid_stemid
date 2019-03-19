#!/usr/local/bin/Rscript

task <- dyncli::main()
# task = dyncli::main(
#   c("--dataset", "/code/example.h5", "--output", "/mnt/output"),
#   "/code/definition.yml"
# )

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(RaceID)

#   ____________________________________________________________________________
#   Load data                                                               ####

parameters <- task$parameters
counts <- as.matrix(task$counts)

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# initialize SCseq object with transcript expression
sc <- SCseq(data.frame(t(counts), check.names = FALSE))

# filtering of expression data
sc <- sc %>% filterdata(
  mintotal = 1,
  minexpr = 0,
  minnumber = 0,
  knn = parameters$knn,
  ccor = parameters$ccor
)

# compute pairwise distances
sc <- sc %>% compdist(
  metric = parameters$metric,
  FSelect = FALSE
)

# perform clustering
parameters$clustnr <- min(parameters$clustnr, ceiling(ncol(sc@expdata)/5))
sc <- sc %>% clustexp(
  sat = parameters$sat,
  samp = parameters$samp,
  cln = parameters$cln,
  clustnr = parameters$clustnr,
  bootnr = parameters$bootnr,
  FUNcluster = parameters$FUNcluster
)

# detect outliers and redefine clusters
sc <- sc %>% findoutliers(
  probthr = parameters$probthr,
  outminc = parameters$outminc,
  outlg = parameters$outlg,
  outdistquant = parameters$outdistquant
)

# compute t-SNE map
sc <- sc %>% comptsne(
  initial_cmd = parameters$initial_cmd,
  perplexity = parameters$perplexity
)

# initialization
ltr <- Ltree(sc)

# computation of the entropy
ltr <- ltr %>% compentropy()

# computation of the projections for all cells
ltr <- ltr %>% projcells(
  cthr = parameters$cthr,
  nmode = parameters$nmode,
  knn = parameters$projcells_knn,
  fr = parameters$fr
)

# computation of the projections for all cells after randomization
ltr <- ltr %>% projback(
  pdishuf = parameters$pdishuf,
  fast = parameters$fast
)

# assembly of the lineage tree
ltr <- ltr %>% lineagegraph()

# compute p-values for link significance
ltr <- ltr %>% comppvalue(
  pthr = parameters$pthr
)


# compute p value
ltr <- ltr %>% comppvalue()

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# collect information on dimreds
dimred_milestones <- ltr@ldata$cnl %>% as.matrix
rownames(dimred_milestones) <- paste0("M", ltr@ldata$m)
dimred <- ltr@ltcoord %>% na.omit
milestone_ids <- rownames(dimred_milestones)
grouping <- paste0("M", ltr@ldata$lp[rownames(dimred)])

# calculate distance between milestones
dist_milestones <- as.matrix(dist(dimred_milestones))

# fetch milestone network by filtering the linkscore
milestone_network <- ltr@cdata$linkscore %>%
  as.matrix() %>%
  reshape2::melt(varnames = c("from", "to"), value.name = "linkscore") %>%
  na.omit() %>%
  mutate_at(c("from", "to"), ~gsub("cl.", "M", ., fixed = TRUE)) %>%
  filter(linkscore >= parameters$scthr) %>%
  mutate(
    length =  dist_milestones[cbind(from, to)],
    directed = FALSE
  ) %>%
  dplyr::select(from, to, length, directed)

#   ____________________________________________________________________________
#   Save output                                                             ####

output <- dynwrap::wrap_data(cell_ids = rownames(dimred)) %>%
  dynwrap::add_dimred_projection(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    dimred = dimred,
    dimred_milestones = dimred_milestones
  ) %>%
  dynwrap::add_timings(checkpoints)

dyncli::write_output(output, task$output)
