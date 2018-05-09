rm(list = ls())

library(dplyr)
library(ggplot2)
library(stringr)
library(corrplot)
library(igraph)

load("workspace.RData")
set.seed(24)

corr.plot <- F
dir.create("moa_plots")
quartz(width = 15, height = 12)

moas <- metadata %>%
  group_by(Metadata_moa) %>%
  tally() %>%
  filter(n >= 5) %>%
  select(Metadata_moa) %>%
  as.matrix() %>%
  as.vector()

moas <- "adrenergic receptor antagonist"

scale_sim <- function(sm) {
  sm <- sm %>%
    reshape2::melt() %>%
    mutate(value = 1*(ecdf(value)(value) - 0.0)) %>%
    reshape2::acast(Var1 ~ Var2, fun.aggregate = function(x) mean(x, na.rm = T))
  
  sm[sm >= 1] <- 1
  sm[sm <= 0] <- 0
  return(sm)
}

cr.median.mad.cov.nrm.main <- scale_sim(cr.median.mad.cov)
cr.mean.nrm.main <- scale_sim(cr.mean)
cr.mad.nrm.main <- scale_sim(cr.mad)
cr.cov.nrm.main <- scale_sim(cr.cov)

for (moa in setdiff(moas, NA)) {
  brd.sample <- metadata %>% 
    filter(Metadata_moa == moa) %>%
    select(Metadata_broad_sample) %>%
    unique %>%
    as.matrix() %>%
    as.vector()

  cr.median.mad.cov.nrm <- cr.median.mad.cov.nrm.main[brd.sample, brd.sample]
  cr.mean.nrm <- cr.mean.nrm.main[brd.sample, brd.sample]
  cr.mad.nrm <- cr.mad.nrm.main[brd.sample, brd.sample]
  cr.cov.nrm <- cr.cov.nrm.main[brd.sample, brd.sample]
  
  method <- "original" ## hclust or original
  hcl <- hclust(as.dist(1 - cr.median.mad.cov.nrm), method = "average")
  ord <- hcl$labels[hcl$order]
  cr.median.mad.cov.nrm <- cr.median.mad.cov.nrm[ord, ord]
  cr.mean.nrm <- cr.mean.nrm[ord, ord]
  cr.mad.nrm <- cr.mad.nrm[ord, ord]
  cr.cov.nrm <- cr.cov.nrm[ord, ord]
  
   rownames(cr.median.mad.cov.nrm) <- rep("", length(brd.sample))
   colnames(cr.median.mad.cov.nrm) <- rep("", length(brd.sample))
  rownames(cr.mean.nrm) <- rep("", length(brd.sample))
  colnames(cr.mean.nrm) <- rep("", length(brd.sample))
  rownames(cr.mad.nrm) <- rep("", length(brd.sample))
  colnames(cr.mad.nrm) <- rep("", length(brd.sample))
  rownames(cr.cov.nrm) <- rep("", length(brd.sample))
  colnames(cr.cov.nrm) <- rep("", length(brd.sample))

  n <- NROW(cr.median.mad.cov.nrm)

    #par(mfrow=c(1, 4))
    #corrplot::corrplot(cr.median.mad.cov.nrm[1:n, 1:n], method = "color", order = method, title = "median+MAD+cov.", mar=c(0,0,20,0))
    #corrplot::corrplot(cr.mean.nrm[1:n, 1:n], method = "color", order = method, title = "median", mar=c(0,0,20,0))
    #corrplot::corrplot(cr.mad.nrm[1:n, 1:n], method = "color", order = method, title = "MAD", mar=c(0,0,20,0))
    #corrplot::corrplot(cr.cov.nrm[1:n, 1:n], method = "color", order = method, title = "cov.", mar=c(0,0,20,0))
    #dev.print(device = pdf, sprintf("moa_plots/%s.pdf", moa))

    png(file=sprintf("moa_plots/graph_%s_top95.png", moa),width=15*27*5,height=12*27*5)
    par(mfrow=c(1, 4))
    A <- 1.0 * (cr.median.mad.cov.nrm[1:n, 1:n] > 0.95)
    diag(A) <- 0
    g <- graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected")
    l1 <- layout_nicely(g)
    plot(g, layout = l1, vertex.size=3, edge.width = 6, vertex.label.dist=1)
    
    A <- 1.0 * (cr.mean.nrm[1:n, 1:n] > 0.95)
    diag(A) <- 0
    g <- graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected")
    plot(g, layout = l1, vertex.size=3, edge.width = 6, vertex.label.dist=1)
    
    A <- 1.0 * (cr.mad.nrm[1:n, 1:n] > 0.95)
    diag(A) <- 0
    g <- graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected")
    plot(g, layout = l1, vertex.size=3, edge.width = 6, vertex.label.dist=1)
    
    A <- 1.0 * (cr.cov.nrm[1:n, 1:n] > 0.95)
    diag(A) <- 0
    g <- graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected")
    plot(g, layout = l1, vertex.size=3, edge.width = 6, vertex.label.dist=1)
    dev.off()
    
    png(file=sprintf("moa_plots/graph_%s_weighted.png", moa),width=15*27*5,height=12*27*5)
    sigma_weight <- 20
    par(mfrow=c(1, 4))
    A <- exp((-1 + cr.median.mad.cov.nrm[1:n, 1:n]) * sigma_weight) * 2
    diag(A) <- 0
    g <- graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected", weighted = TRUE)
    E(g)$width <- E(g)$weight * 3 
    plot(g, layout = l1, vertex.size=3, vertex.label.dist=1)
    
    A <- exp((-1 + cr.mean.nrm[1:n, 1:n]) * sigma_weight) * 2
    diag(A) <- 0
    g <- graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected", weighted = TRUE)
    E(g)$width <- E(g)$weight * 3 
    plot(g, layout = l1, vertex.size=3, vertex.label.dist=1)
    
    A <- exp((-1 + cr.mad.nrm[1:n, 1:n]) * sigma_weight) * 2
    diag(A) <- 0
    g <- graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected", weighted = TRUE)
    E(g)$width <- E(g)$weight * 3 
    plot(g, layout = l1, vertex.size=3, vertex.label.dist=1)
    
    A <- exp((-1 + cr.cov.nrm[1:n, 1:n]) * sigma_weight) * 2
    diag(A) <- 0
    g <- graph_from_adjacency_matrix(adjmatrix = A, mode = "undirected", weighted = TRUE)
    E(g)$width <- E(g)$weight * 3 
    plot(g, layout = l1, vertex.size=3, vertex.label.dist=1)
    
    dev.off()
}