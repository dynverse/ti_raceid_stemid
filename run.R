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

params <- task$params
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
  knn = params$knn,
  ccor = params$ccor
)

# compute pairwise distances
sc <- sc %>% compdist(
  metric = params$metric,
  FSelect = FALSE
)

# perform clustering
params$clustnr <- min(params$clustnr, ceiling(ncol(sc@expdata)/5))
sc <- sc %>% clustexp(
  sat = params$sat,
  samp = params$samp,
  cln = params$cln,
  clustnr = params$clustnr,
  bootnr = params$bootnr,
  FUNcluster = params$FUNcluster
)

# detect outliers and redefine clusters
sc <- sc %>% findoutliers(
  probthr = params$probthr,
  outminc = params$outminc,
  outlg = params$outlg,
  outdistquant = params$outdistquant
)

# compute t-SNE map
sc <- sc %>% comptsne(
  initial_cmd = params$initial_cmd,
  perplexity = params$perplexity
)

# initialization
ltr <- Ltree(sc)

# computation of the entropy
ltr <- ltr %>% compentropy()

# computation of the projections for all cells
ltr <- ltr %>% projcells(
  cthr = params$cthr,
  nmode = params$nmode,
  knn = params$projcells_knn,
  fr = params$fr
)

# computation of the projections for all cells after randomization
ltr <- ltr %>% projback(
  pdishuf = params$pdishuf,
  fast = params$fast
)

# assembly of the lineage tree
ltr <- ltr %>% lineagegraph()

# compute p-values for link significance
ltr <- ltr %>% comppvalue(
  pthr = params$pthr
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
  filter(linkscore >= params$scthr) %>%
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
