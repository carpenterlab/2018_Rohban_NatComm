rm(list = ls())

library(dplyr)
library(stringr)
library(ggplot2)

MOA <- "hdac inhibitor"
line.width <- 3
snf.on.mean <- F

sanitize <- function(a) {
  a.t <- c(a$Var1, a$Var2) %>% unique
  a.x <- a %>%
    select(Var1, Var2, value) %>%
    rbind(data.frame(Var1 = a.t, Var2 = a.t, value = 1)) %>% 
    reshape2::acast("Var1 ~ Var2") 
  a.x[lower.tri(a.x, diag = F)] <- t(a.x)[lower.tri(t(a.x), diag = F)]
  meta <- rbind(a %>% select(Var1, Metadata_moa.x) %>% rename(Var = Var1, Metadata_moa = Metadata_moa.x), a %>% select(Var2, Metadata_moa.y) %>% rename(Var = Var2, Metadata_moa = Metadata_moa.y)) %>% unique 
  a <- a.x %>%
    reshape2::melt() %>%
    left_join(meta, by = c("Var1" = "Var")) %>%
    left_join(meta, by = c("Var2" = "Var"))
  return(a)
}

a <- readRDS("sm_median_mad.rds")
b <- readRDS("sm_median_mad_cov.rds")

a <- sanitize(a)
b <- sanitize(b)

if (snf.on.mean) {
  cr.mean <- a %>%
    select(Var1, Var2, value) %>%
    group_by(Var1, Var2) %>%
    summarise(value = max(value)) %>%
    reshape2::acast("Var1 ~ Var2") 
  
  d <- apply(cr.mean, 1, function(x) !(sum(is.na(x)) >= (NROW(cr.mean) -1 )))
  cr.mean <- cr.mean[d, d]
  
  k.snf <- 7
  af.1 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = 0.5)
  af.2 <- SNFtool::affinityMatrix(Diff = 1 - cr.mean, K = k.snf, sigma = 0.5)
  af.snf <- SNFtool::SNF(list(af.1, af.2), K = k.snf, t = 10)
  rownames(af.snf) <- rownames(af.1)
  colnames(af.snf) <- colnames(af.1)
  cr.mean <- af.snf
  
  a <- cr.mean %>%
    reshape2::melt() %>%
    left_join(a %>% select(Metadata_moa.x, Metadata_moa.y, Var1, Var2), by = c("Var1", "Var2"))
}

a <- a %>% filter(!is.na(Metadata_moa.x) & !is.na(Metadata_moa.y) & Metadata_moa.x != "" & Metadata_moa.y != "" & Var1 != Var2)
b <- b %>% filter(!is.na(Metadata_moa.x) & !is.na(Metadata_moa.y) & Metadata_moa.x != "" & Metadata_moa.y != "" & Var1 != Var2)

n <- a$Var1 %>% unique %>% length
a.x <- a %>% filter(str_detect(Metadata_moa.x, MOA)) %>% arrange(-value) %>% group_by(Var1) %>% mutate(Var2.x = seq(n - 1, 1, -1)) %>% ungroup() %>%
  mutate(same.moa = str_detect(Metadata_moa.y, MOA))

b.x <- b %>% filter(str_detect(Metadata_moa.x, MOA)) %>% arrange(-value) %>% group_by(Var1) %>% mutate(Var2.x = seq(n - 1 , 1, -1)) %>% ungroup() %>%
  mutate(same.moa = str_detect(Metadata_moa.y, MOA))

a.y <- a.x %>% select(Var1, Var2.x, same.moa)
b.y <- b.x %>% select(Var1, Var2.x, same.moa)

expand <- function(x) {
  if (as.vector(as.matrix(x[1,"same.moa"]))) {
    return(data.frame(Var2.xx = seq(from = as.vector(as.matrix(x[,"Var2.x"])), to = as.vector(as.matrix(x[,"Var2.x"])) + line.width, by = 1)))
  } 
  return(data.frame(Var2.xx = as.vector(as.matrix(x[,"Var2.x"]))))
}

aa.y <- a.y %>%
  group_by(Var1, Var2.x, same.moa) %>%
  do(expand(.)) %>%
  ungroup()

a.y <- aa.y %>%
  select(-Var2.x) %>%
  rename(Var2.x = Var2.xx)

bb.y <- b.y %>%
  group_by(Var1, Var2.x, same.moa) %>%
  do(expand(.)) %>%
  ungroup()

b.y <- bb.y %>%
  select(-Var2.x) %>%
  rename(Var2.x = Var2.xx)


a.t <- a.y %>% 
  filter(same.moa) %>%
  group_by(Var1) %>%
  summarise(vl = max(Var2.x)) %>%
  arrange(-vl) %>%
  select(Var1) %>%
  as.matrix() %>% 
  as.vector()

b.t <- b.y %>% 
  filter(same.moa) %>%
  group_by(Var1) %>%
  summarise(vl = max(Var2.x)) %>%
  arrange(-vl) %>%
  select(Var1) %>%
  as.matrix() %>% 
  as.vector()

tit1 <- paste(MOA, "median+mad (SNF)", sep = "-", collapse = "")
tit2 <- paste(MOA, "median+mad+cov. (SNF)", sep = "-", collapse = "")

g1 <- ggplot(a.y, aes(x = Var1, y = Var2.x)) + 
  geom_raster(aes(fill=same.moa)) + 
  scale_fill_manual(values=c("white", "black")) + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) + 
  xlab("") + 
  ylab("") +
  theme(legend.position="none") + 
  ggtitle(tit1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = a.t) 

g2 <- ggplot(b.y, aes(x = Var1, y = Var2.x)) + 
  geom_raster(aes(fill=same.moa)) + 
  scale_fill_manual(values=c("white", "black")) + 
  scale_y_continuous(limits = c(0, NA)) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) + 
  xlab("") + 
  ylab("") +
  theme(legend.position="none") + 
  ggtitle(tit2) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = b.t) 

ggsave(sprintf("../figs/%s.png", tit1), g1, width = 10, height = 8)
ggsave(sprintf("../figs/%s.png", tit2), g2, width = 10, height = 8)
