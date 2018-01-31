library(dplyr)
library(stringr)
library(reshape2)
library(doParallel)

same.moa <- function(moa.list.1, moa.list.2) {
  if (is.na(moa.list.1) || is.na(moa.list.2) || moa.list.1 == "" || moa.list.2 == "") 
    return(FALSE)
  moas.1 <- str_split(moa.list.1, "\\|")[[1]]
  moas.2 <- str_split(moa.list.2, "\\|")[[1]]
  return(any(moas.1 %in% moas.2) | any(moas.2 %in% moas.1))
}

same.moa <- Vectorize(same.moa)

enrichment_top_conn <- function(sm, metadata, top.perc = 0.95) {
  sm <- sm %>% 
    reshape2::melt() %>% 
    filter(as.character(Var1) < as.character(Var2) & 
             Var1 != "DMSO" &
             Var2 != "DMSO") %>%
    left_join(., 
              metadata, 
              by = c("Var1" = "Metadata_broad_sample")) %>%
    left_join(., 
              metadata, 
              by = c("Var2" = "Metadata_broad_sample")) %>%
    mutate(same.moa = same.moa(Metadata_moa.x, Metadata_moa.y))
  
  thr <- quantile(sm$value, top.perc, na.rm = T)
  
  v11 <- sm %>%
    filter(value > thr & same.moa) %>%
    NROW
  
  v12 <- sm %>%
    filter(value > thr & !same.moa) %>%
    NROW
  
  v21 <- sm %>%
    filter(value < thr & same.moa) %>%
    NROW
  
  v22 <- sm %>%
    filter(value < thr & !same.moa) %>%
    NROW
  
  return(fisher.test(x = rbind(c(v11, v12), c(v21, v22)), 
                     alternative = "greater"))
}

contains.moa <- function(moa.list, moa) {
  moa %in% str_split(moa.list, "\\|")[[1]]  
}

contains.moa <- Vectorize(contains.moa)

moa_recall <- function(sm, metadata, n.cores = 1, N = 1000) {
  doParallel::registerDoParallel(cores = n.cores)
  
  moas <- unlist(lapply(metadata$Metadata_moa, function(x) str_split(x, "\\|")[[1]]))
  moas <- unique(moas)
  moas <- setdiff(moas, NA)
  moas <- setdiff(moas, "")

  group_recall <- function(sm, brds) {
    diag(sm[brds, brds]) <- NA
    
    median(apply(sm[brds, brds], 1, function(x) median(x, na.rm = T)), na.rm = T)
    #sm[brds, brds] %>%
    #  as.dist() %>%
    #  mean(., na.rm = T)
  }  
  
  brds.ref <- colnames(sm)
  
  res <- foreach (moa = moas, .combine = rbind) %dopar% {
    brds.moa <- metadata %>%
      filter(contains.moa(Metadata_moa, moa)) %>% 
      select(Metadata_broad_sample) %>% 
      as.matrix() %>%
      as.vector()
    
    brds.moa <- intersect(brds.moa, brds.ref)
    
    if (length(brds.moa) < 2) {
      return(data.frame(MOA = moa, p.value = NA))
    }
      
    x <- group_recall(sm = sm, brds = brds.moa)
    
    nulls <- foreach (j = 1:N, .combine = rbind) %do% {
      brds <- metadata %>%
        sample_n(NROW(brds.moa)) %>%
        select(Metadata_broad_sample) %>% 
        as.matrix() %>%
        as.vector()
      
      brds <- intersect(brds, brds.ref)
      
      group_recall(sm = sm, brds = brds)
    }
    
    data.frame(MOA = moa, p.value = 1 - ecdf(nulls)(x))
  }
  
  return(res)
}

moa_recall_fisher <- function(sm, metadata, thr) {
  moas <- unlist(lapply(metadata$Metadata_moa, function(x) str_split(x, "\\|")[[1]]))
  moas <- unique(moas)
  moas <- setdiff(moas, NA)
  moas <- setdiff(moas, "")
 
  brds.ref <- colnames(sm)
  
  res <- foreach (moa = moas, .combine = rbind) %dopar% {
    brds.moa <- metadata %>%
      filter(contains.moa(Metadata_moa, moa)) %>% 
      select(Metadata_broad_sample) %>% 
      as.matrix() %>%
      as.vector()
    
    brds.moa <- intersect(brds.moa, brds.ref)
    brds.other <- setdiff(brds.ref, brds.moa)
    
    if (length(brds.moa) < 6) {
      return(data.frame(MOA = moa, p.value = NA))
    }
    
    x <- sm[brds.moa, brds.moa] %>% as.dist() %>% as.vector()
    y <- sm[brds.moa, brds.other] %>% as.vector()
    
    z1 <- x
    z2 <- y
    
    v11 <- which(z1 > thr) %>% length
    v12 <- which(z1 < thr) %>% length
    v21 <- which(z2 > thr) %>% length
    v22 <- which(z2 < thr) %>% length
    
    V <- rbind(c(v11, v12), c(v21, v22))
    
    fsh <- fisher.test(V, alternative = "greater")
    data.frame(MOA = moa, p.value = fsh$p.value)
  }
  
  return(res)
}

signif.test <- function(x, k0) {
  moa.q <- str_split(x[1, "Metadata_moa.x"] %>% as.matrix() %>% as.vector(), "\\|") %>% unlist %>% unique 
  moas <- data.frame(Metadata_moa.y = (str_split(x[1:k0, "Metadata_moa.y"] %>% as.matrix() %>% as.vector(), "\\|") %>% unlist()))
  moas <- moas %>%
    group_by(Metadata_moa.y) %>%
    tally() %>%
    ungroup() %>%
    arrange(-n) 

  moas.left <- data.frame(Metadata_moa.y = (str_split(x[(k0+1):NROW(x), "Metadata_moa.y"] %>% as.matrix() %>% as.vector(), "\\|") %>% unlist())) %>%
    group_by(Metadata_moa.y) %>%
    tally() %>%
    ungroup() %>%
    arrange(-n) 
  
  df <- data.frame(MOA = c(), p.val = c())
  for (i in 1:NROW(moas)) {
    moa.x <- moas$Metadata_moa.y[i]
    v11 <- moas$n[i]
    v12 <- sum(moas$n) - v11
    
    v21 <- moas.left %>% 
      filter(Metadata_moa.y == as.character(moa.x)) %>%
      select(n) %>%
      as.matrix() %>%
      as.vector()
    
    if (length(v21) == 0) {
      v21 <- 0
    }
    v22 <- sum(moas.left$n) - v21
    V <- rbind(c(v11, v12), c(v21, v22))
    #print(v11)
    #print(v12)
    #print(v21)
    #print(v22)
    df.s <- data.frame(MOA = moa.x, p.val = fisher.test(V, "greater")$p.value)
    df <- rbind(df, df.s)
  }
  
  df <- df %>%
    mutate(p.val = stats::p.adjust(p.val)) %>%
    arrange(p.val) %>%
    slice(1) 
}

cmpd_classification <- function(sm, metadata, k0 = 5) {
  sm <- sm %>% 
    reshape2::melt() %>% 
    filter(Var1 != Var2 &
             Var1 != "DMSO" &
             Var2 != "DMSO") %>%
    left_join(., 
              metadata, 
              by = c("Var1" = "Metadata_broad_sample")) %>%
    left_join(., 
              metadata, 
              by = c("Var2" = "Metadata_broad_sample")) %>%
    filter(!is.na(Metadata_moa.x) & !is.na(Metadata_moa.y))
    #mutate(same.moa = same.moa(Metadata_moa.x, Metadata_moa.y))

  cmpd.true.pos <- sm %>% 
    arrange(-value) %>%
    group_by(Var1, Metadata_moa.x) %>%
    do(signif.test(.[], k0)) %>%
    #slice(1:k0) %>%
    #summarise(pass = sum(same.moa(Metadata_moa.x, Metadata_moa.y)) >= thr) %>%
    ungroup() %>%
    mutate(p.val = stats::p.adjust(p.val)) %>%
    mutate(pass = (p.val < 0.05) & same.moa(Metadata_moa.x, MOA)) 
  
  #return(cmpd.true.pos/(sm$Var1 %>% unique %>% length))
  return(cmpd.true.pos)
}