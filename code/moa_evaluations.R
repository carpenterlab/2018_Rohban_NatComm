library(dplyr)
library(stringr)
library(reshape2)
library(doParallel)
library(htmlTable)

sim_normalize_rect <- function(sim_mat) {
  # sim_mat_norm <- apply(sim_mat, 1, function(x) (ecdf(x)(x)))
  # sim_mat_norm <- (sim_mat_norm + t(sim_mat_norm))/2
  # rownames(sim_mat_norm) <- rownames(sim_mat)
  # colnames(sim_mat_norm) <- colnames(sim_mat)
  # return(sim_mat_norm)
  
  sm <- sim_mat[upper.tri(sim_mat)]
  sim_mat <- (sim_mat - median(sm))/mad(sm)
  sim_mat <- sim_mat/quantile(sim_mat, 0.999) * 0.999
  sim_mat[(sim_mat > 1)] <- 1
  sim_mat[(sim_mat < -1)] <- -1
  diag(sim_mat) <- 1
  
  # sm.melt <- sim_mat %>% 
  #   reshape2::melt()
  # 
  # sm.melt <- sm.melt %>% 
  #   mutate(value = ecdf(value)(value))
  # 
  # sim_mat <- sm.melt %>%
  #   reshape2::acast(Var1 ~ Var2)
  # 
  # diag(sim_mat) <- 1
  return(sim_mat)
}

sim_normalize <- function(sim_mat) {
  # sim_mat_norm <- apply(sim_mat, 1, function(x) (ecdf(x)(x)))
  # sim_mat_norm <- (sim_mat_norm + t(sim_mat_norm))/2
  # rownames(sim_mat_norm) <- rownames(sim_mat)
  # colnames(sim_mat_norm) <- colnames(sim_mat)
  # return(sim_mat_norm)
  # sm <- sim_mat[upper.tri(sim_mat)]
  # sim_mat <- (sim_mat - median(sm))/mad(sm)
  # sim_mat <- sim_mat/quantile(sim_mat, 0.999) * 0.999
  # diag(sim_mat) <- 1
  # sim_mat[(sim_mat > 1)] <- 1
  # sim_mat[(sim_mat < -1)] <- -1
  
  sm <- sim_mat[upper.tri(sim_mat)]
  sim_mat <- (sim_mat - median(sm))/mad(sm)
  sim_mat <- sim_mat/quantile(sim_mat, 0.999) * 0.999
  sim_mat[(sim_mat > 1)] <- 1
  sim_mat[(sim_mat < -1)] <- -1
  diag(sim_mat) <- 1
  
  # sm.melt <- sim_mat %>%
  #   reshape2::melt()
  # 
  # sm.melt <- sm.melt %>%
  #   mutate(value = ecdf(value)(value))
  # 
  # sim_mat <- sm.melt %>%
  #   reshape2::acast(Var1 ~ Var2)
  # 
  # diag(sim_mat) <- 1
  return(sim_mat)
}

same.moa <- function(moa.list.1, moa.list.2) {
  if (is.na(moa.list.1) || is.na(moa.list.2) || moa.list.1 == "" || moa.list.2 == "") 
    return(FALSE)
  moa.list.1 <- str_to_lower(moa.list.1)
  moa.list.2 <- str_to_lower(moa.list.2)
  
  moas.1 <- str_split(moa.list.1, "\\|")[[1]]
  moas.2 <- str_split(moa.list.2, "\\|")[[1]]
  return(any(moas.1 %in% moas.2) | any(moas.2 %in% moas.1))
}

same.moa <- Vectorize(same.moa)

perpare_sm <- function(sm, metadata) {
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
    filter(!is.na(Metadata_moa.x) & !is.na(Metadata_moa.y) & Metadata_moa.x != "" & Metadata_moa.y != "") %>%
    mutate(same.moa = same.moa(Metadata_moa.x, Metadata_moa.y))
  return(sm)
}

enrichment_top_conn <- function(sm, metadata, top.perc = 0.95, not.same.batch = F) {
  if (not.same.batch) {
    sm <- sm %>%
      filter((is.na(Metadata_Plate_Map_Name.x) & !is.na(Metadata_Plate_Map_Name.y))
             | (is.na(Metadata_Plate_Map_Name.y) & !is.na(Metadata_Plate_Map_Name.x))
             | (Metadata_Plate_Map_Name.x != Metadata_Plate_Map_Name.y))
  }
  
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

enrichment_top_conn_cross <- function(sm, metadata1, metadata2, top.perc = 0.95) {
  sm <- sm %>% 
    reshape2::melt() %>% 
    filter(Var1 != "DMSO" &
           Var2 != "DMSO") %>%
    left_join(., 
              metadata1, 
              by = c("Var1" = "Metadata_broad_sample")) %>%
    left_join(., 
              metadata2, 
              by = c("Var2" = "Metadata_broad_sample")) %>%
    filter(!is.na(Metadata_moa.x) & !is.na(Metadata_moa.y) & Metadata_moa.x != "" & Metadata_moa.y != "") %>%
    mutate(same.moa = same.moa(str_to_lower(Metadata_moa.x), str_to_lower(Metadata_moa.y)))
  
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
  
  V <- rbind(c(v11, v12), c(v21, v22))
  return(list(top.conn.valid = sm %>%
                filter(value > thr & same.moa), V = V, test = fisher.test(x = V, 
                     alternative = "greater")))
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
    #diag(sm[brds, brds]) <- NA

    #median(apply(sm[brds, brds], 1, function(x) median(x, na.rm = T)), na.rm = T)
    sm[brds, brds] %>%
      as.dist() %>%
      mean(., na.rm = T)
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
        sample_n(length(brds.moa)) %>%
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

moa_recall_fisher <- function(sm, metadata) {
  moas <- unlist(lapply(metadata$Metadata_moa, function(x) str_split(x, "\\|")[[1]]))
  moas <- unique(moas)
  moas <- setdiff(moas, NA)
  moas <- setdiff(moas, "")
 
  thr <- quantile(sm %>% as.dist(), 0.99, na.rm = T)
  brds.ref <- colnames(sm)
  
  res <- foreach (moa = moas, .combine = rbind) %dopar% {
    brds.moa <- metadata %>%
      filter(contains.moa(Metadata_moa, moa)) %>% 
      select(Metadata_broad_sample) %>% 
      as.matrix() %>%
      as.vector()
    
    brds.moa <- intersect(brds.moa, brds.ref)
    brds.other <- setdiff(brds.ref, brds.moa)
    
    if (length(brds.moa) < 2) {
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

signif.test <- function(x, thr) {
  sz <- x %>% filter(value > thr) %>% NROW()
  
  if (sz < 1) {
    return(data.frame(MOA = NA, p.val = 1))  
  }
  
  moa.q <- str_split(x[1, "Metadata_moa.x"] %>% as.matrix() %>% as.vector(), "\\|") %>% unlist %>% unique 
  moas <- data.frame(Metadata_moa.y = (str_split(x %>% filter(value > thr) %>% select(Metadata_moa.y) %>% as.matrix() %>% as.vector(), "\\|") %>% unlist()))
  moas <- moas %>%
    group_by(Metadata_moa.y) %>%
    tally() %>%
    ungroup() %>%
    arrange(-n) 

  moas.left <- data.frame(Metadata_moa.y = (str_split(x %>% filter(value < thr) %>% select(Metadata_moa.y) %>% as.matrix() %>% as.vector(), "\\|") %>% unlist())) %>%
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

cmpd_classification_curve <- function(sm, metadata, k0 = 5, not.same.batch = F) {
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
    filter(!is.na(Metadata_moa.x) & !is.na(Metadata_moa.y) & Metadata_moa.x != "" & Metadata_moa.y != "")
  
  thr <- quantile(sm$value, 0.99, na.rm = T)

  if (not.same.batch) {
    sm <- sm %>%
      filter((is.na(Metadata_Plate_Map_Name.x) & !is.na(Metadata_Plate_Map_Name.y))
             | (is.na(Metadata_Plate_Map_Name.y) & !is.na(Metadata_Plate_Map_Name.x))
             | (Metadata_Plate_Map_Name.x != Metadata_Plate_Map_Name.y))
  }
  
  return(
    sm %>% 
      arrange(-value) %>% 
      group_by(Var1, Metadata_moa.x) %>% 
      slice(1:k0) %>% 
      filter(value > thr) %>%
      summarise(num.same.moa.conn = sum(same.moa(str_to_lower(Metadata_moa.x), str_to_lower(Metadata_moa.y))))
  )
}

cmpd_knn_classification <- function(sm, metadata, k0 = 5, not.same.batch = F) {
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
    filter(!is.na(Metadata_moa.x) & !is.na(Metadata_moa.y) & Metadata_moa.x != "" & Metadata_moa.y != "")
  
  if (not.same.batch) {
    sm <- sm %>%
      filter((is.na(Metadata_Plate_Map_Name.x) & !is.na(Metadata_Plate_Map_Name.y))
             | (is.na(Metadata_Plate_Map_Name.y) & !is.na(Metadata_Plate_Map_Name.x))
             | (Metadata_Plate_Map_Name.x != Metadata_Plate_Map_Name.y))
  }
  
  thr <- quantile(sm$value, 0.99, na.rm = T)
  
  return(
    sm %>% 
      arrange(-value) %>% 
      group_by(Var1, Metadata_moa.x) %>% 
      slice(1:k0) %>% 
      filter(value > thr) %>%
      summarise(pass = ifelse(sum(same.moa(str_to_lower(Metadata_moa.x), str_to_lower(Metadata_moa.y))) >= 2, T, F)) %>% 
      filter(pass) 
  )
}

  
cmpd_classification <- function(sm, metadata, thr.perc, not.same.batch = F) {
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
    filter(!is.na(Metadata_moa.x) & !is.na(Metadata_moa.y) & Metadata_moa.x != "" & Metadata_moa.y != "")
    #mutate(same.moa = same.moa(Metadata_moa.x, Metadata_moa.y))

  if (not.same.batch) {
    sm <- sm %>%
      filter((is.na(Metadata_Plate_Map_Name.x) & !is.na(Metadata_Plate_Map_Name.y))
             | (is.na(Metadata_Plate_Map_Name.y) & !is.na(Metadata_Plate_Map_Name.x))
             | (Metadata_Plate_Map_Name.x != Metadata_Plate_Map_Name.y))
  }
  
  thr <- quantile(sm$value, thr.perc, na.rm = T)
  
  cmpd.true.pos <- sm %>% 
    arrange(-value) %>%
    group_by(Var1, Metadata_moa.x) %>%
    do(signif.test(.[], thr)) %>%
    ungroup() %>%
    filter(!is.na(MOA)) %>%
    mutate(p.val = stats::p.adjust(p.val)) %>%
    mutate(pass = (p.val < 0.05) & same.moa(Metadata_moa.x, MOA)) 
  
  return(cmpd.true.pos)
}

protein.interaction.comparison <- function(sm, interaction.data, extra.interaction.data, metadata, quant, suffix = "") {
  ## use gather_protein_interaction_data.R to obtain updated versions of the following PPI data
  ppi <- interaction.data
  ppi.detailed <- extra.interaction.data
  cr.org <- sm
  cr.org <- cr.org %>% 
    reshape2::melt() %>% 
    inner_join(metadata, by = c("Var1" = "Metadata_broad_sample")) %>%
    inner_join(metadata, by = c("Var2" = "Metadata_broad_sample")) %>%
    select(Metadata_Treatment.x, Metadata_Treatment.y, value) %>%
    reshape2::acast(Metadata_Treatment.x ~ Metadata_Treatment.y)

  cr <- cr.org
  
  genes <- rownames(cr) %>% lapply(., function(x) str_split(x, "_")[[1]][1]) %>% unlist
  alleles <- rownames(cr) %>% lapply(., function(x) str_split(x, "_")[[1]][2]) %>% unlist %>% 
    lapply(., function(x) str_split(x, "\\.")[[1]][1]) %>% unlist
  
  u <- outer(genes, genes, function(x1, x2) return(x1 == x2))
  v <- outer(alleles, alleles, function(x1, x2) return(x1 == "WT" & x2 == "WT"))
  z <- outer(rownames(cr), rownames(cr), function(x, y) return(x < y))
  
  w <- u & v & z
  rownames(w) <- rownames(cr)
  colnames(w) <- colnames(cr)
  cr.melt <- melt(cr)
  w.melt <- melt(w)
  wt.cr.melt <- cr.melt[which(w.melt$value),]
  
  trt.to.go <- c()
  for (i in 1:length(genes)) {
    if (alleles[i] == "WT" && (length(trt.to.go) == 0 || !(genes[i] %in% genes[trt.to.go]))) {
      trt.to.go <- c(trt.to.go, i)
    } 
  }
  cr <- cr[trt.to.go, trt.to.go]
  
  cr.melt <- cr %>% 
    reshape2::melt() 
    
  thr <- quantile(cr.melt$value, quant, na.rm = T)
  
  consistency.score <- function(cr, ppi, verbose = F, thr = 0.5, ppi.detailed = NULL) {
    cr.melt <- cr %>% melt
    cl <- colnames(cr.melt)
    cl[1] <- "Var1"
    cl[2] <- "Var2"
    colnames(cr.melt) <- cl
    
    cr.melt.subset <- cr.melt %>% dplyr::filter((Var1 %>% as.character()) < (Var2 %>% as.character()) & (value) >= thr)
    
    if (verbose) {
      print(min(cr.melt.subset$value))
      print(sprintf("The number of pairs suggested : %d", NROW(cr.melt.subset)))
    }
    genes1 <- lapply(cr.melt.subset[,1] %>% as.matrix %>% as.vector, function(x) str_split(x, "_")[[1]][1]) %>% unlist
    genes2 <- lapply(cr.melt.subset[,2] %>% as.matrix %>% as.vector, function(x) str_split(x, "_")[[1]][1]) %>% unlist
    gene.pairs <- cbind(genes1, genes2) %>% unique
    
    gene.ref <- rownames(ppi)
    gene.pairs <- gene.pairs[which(gene.pairs[,1] %in% gene.ref & gene.pairs[,2] %in% gene.ref),]
    sm <- apply(gene.pairs, 1, function(x) return(max(ppi[x[1], x[2]], ppi[x[1], x[2]]))) %>% sum
    u <- c()
    if (verbose) {
      for (i in 1:NROW(gene.pairs)) {
        if (ppi[gene.pairs[i, 1], gene.pairs[i, 2]] == 1) {
          if (!is.null(ppi.detailed)) {
            v <- ppi.detailed %>% dplyr::filter((Protein.1 == gene.pairs[i, 1] & Protein.2 == gene.pairs[i, 2]) | (Protein.2 == gene.pairs[i, 1] & Protein.1 == gene.pairs[i, 2])) 
            u <- rbind(u, v)
          }
        }
      }
      
      if (!is.null(u)) {
        cut.str <- function(x) {
          x <- as.character(x)
          if (str_length(x) > 100) {
            i <- str_locate(x, "\\)")
            return(str_sub(x, 1, i[1]))
          } else {
            return(x)
          }
        }
        
        u$Evidence <- lapply(u$Evidence, cut.str) %>% unlist
        u %>% htmlTable() %>% cat(., file = paste0("evidence_detailed_", suffix, ".html"))
        u %>% select(Protein.1, Protein.2) %>% unique %>% 
          mutate(Protein.1.x = ifelse(as.character(Protein.1) < as.character(Protein.2), as.character(Protein.1), as.character(Protein.2)), 
                 Protein.2.x = ifelse(as.character(Protein.1) < as.character(Protein.2), as.character(Protein.2), as.character(Protein.1))) %>%
          select(-Protein.1, -Protein.2) %>%
          rename(Protein.1 = Protein.1.x, Protein.2 = Protein.2.x) %>%
          unique %>% 
          htmlTable() %>% cat(., file = paste0("evidence_", suffix, ".html"))
      }
    }
    return(sm)
  }
  
  sig <- consistency.score(cr, ppi, T, thr, ppi.detailed)
  
  cr.melt <- cr %>% melt
  cl <- colnames(cr.melt)
  cl[1] <- "Var1"
  cl[2] <- "Var2"
  colnames(cr.melt) <- cl
  
  cr.melt.subset <- cr.melt %>% dplyr::filter((Var1 %>% as.character()) < (Var2 %>% as.character()) & (value) >= thr)
  
  genes <- rownames(cr) %>% lapply(., function(x) str_split(x, "_")[[1]][1]) %>% unlist
  genes.in <- intersect(colnames(ppi), genes)
  all.verified <- (ppi[genes.in, genes.in] %>% sum())/2
  x2 <- sig
  x1 <- NROW(cr.melt.subset) - x2
  
  data <- (rbind(c(x2, x1), c(all.verified-x2, NROW(cr) * (NROW(cr) - 1)/2 -all.verified-x1)))
  print(data)
  fisher.test(data, alternative = "greater") 
}
