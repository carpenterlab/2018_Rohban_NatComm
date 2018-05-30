library(dplyr)
library(magrittr)
library(cytominer)
library(foreach)
library(stringr)
library(readr)
library(doParallel)
library(ggplot2)

set.seed(24)

source("generate_component_matrix.R")

x <- readRDS("../tmp_sql/cov_profiles_24278.rds")
y <- readr::read_csv("../tmp_sql/24278_normalized.csv")
dmso.ids <- which(y$Metadata_broad_sample == "DMSO")
dmso.prf <- x[dmso.ids, ]
  
dmso.prf <- dmso.prf[1:(NCOL(dmso.prf) - 2)]
x <- x[1:(NCOL(x) - 2)]

mn <- apply(dmso.prf, 2, function(x) mean(x, na.rm = T))
sdv <- apply(dmso.prf, 2, function(x) sd(x, na.rm = T))

x <- scale(x, center = mn, scale = sdv)

ndim <- c(seq(from = 100, to = 500, by = 100), seq(from = 600, to = 1600, by = 200), seq(from = 1700, to = 3700, by = 400))
#ndim <- c(100, 200)

D <- foreach(ndim2 = ndim, .combine = rbind) %do% {
  cr.vals <- foreach (i = 1:20) %do% {
    A <- generate_component_matrix(NCOL(x), n_components = ndim2, density = 0.1)
    B <- generate_component_matrix(NCOL(x), n_components = ndim2, density = 0.1)
    
    x1 <- as.matrix(x) %*% A
    x2 <- as.matrix(x) %*% B
    
    x1 <- as.matrix(x1)
    x2 <- as.matrix(x2)
    
    x1 <- x1[complete.cases(x1), ]
    x2 <- x2[complete.cases(x2), ]
    
    c1 <- cor(t(x1)) 
    c2 <- cor(t(x2)) 
    
    d1 <- c1 %>% 
      reshape2::melt() %>%
      filter(Var1 > Var2) %>%
      filter(value > quantile(value, 0.99)) %>%
      dplyr::select(Var1, Var2)
    
    d2 <- c2 %>% 
      reshape2::melt() %>%
      filter(Var1 > Var2) %>%
      filter(value > quantile(value, 0.99)) %>%
      dplyr::select(Var1, Var2)
    
    n1 <- d1 %>%
      semi_join(d2, by = c("Var1", "Var2")) %>%
      NROW
    
    n2 <- d1 %>%
      NROW()
    
    (n1/n2)
  }
  
  v1 <- cr.vals %>% unlist() %>% mean()
  v2 <- cr.vals %>% unlist() %>% sd()
  #data.frame(ndim = ndim2, mean.overlap = v1, std.overlap = v2)
  data.frame(ndim = ndim2, overlap = cr.vals %>% unlist())
}

D2 <- D
D <- D2 %>% group_by(ndim) %>% summarise(mean.overlap = mean(overlap), std.overlap = sd(overlap)) %>% ungroup()
D3 <- D2 %>% left_join(D, by = "ndim")

g <- ggplot(D3 %>% mutate(mean.overlap = mean.overlap * 100, std.overlap = std.overlap * 100, overlap = overlap * 100), aes(x = ndim, y = mean.overlap)) + 
  geom_errorbar(aes(ymin=mean.overlap-std.overlap/(20^0.5)*2, ymax=mean.overlap+std.overlap/(20^0.5)*2), width=50) + geom_jitter(aes(x = ndim, y = overlap), size = 0.05) + 
  geom_line() + xlab("No. random projections") + ylab("Avg. overlap percentage \n of top 1% connections") + theme(axis.text = element_text(size=20), text = element_text(size=25)) + 
  theme(plot.title = element_text(hjust = 0.5)) 

ggsave("RP_stability.png", g)
ggsave("RP_stability.pdf", g)

