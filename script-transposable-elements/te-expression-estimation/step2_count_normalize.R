
gene_count_normalizer <- function(gene_count, gene_lentable, out_dir){
  rm(list=ls())
  suppressMessages(library(dplyr))

  geneCount = read.table(gene_count, stringsAsFactors=T)
  geneLength = read.table(gene_lentable, stringsAsFactors=T)

  # Optional, remove header row
  #geneCount = geneCount[-1,]
  #geneLength = geneLength[-1,]

  colnames(geneCount) = c("gene_id","raw_read_count")
  colnames(geneLength) = c("gene_id", "length")

  geneCount$raw_read_count = as.numeric(geneCount$raw_read_count)
  geneLength$length = as.numeric(geneLength$length)

  # Remove duplicated feature (gene/TE) entries
  geneCount = unique(geneCount)
  geneLength = unique(geneLength)

  # Merge raw_counts df and length df into a large data frame
  tmp <- inner_join(geneCount, geneLength, by='gene_id')


  # count2fpkm
  fpkm <- data.frame(row.names = tmp$gene_id)
  for (i in 2:(dim(tmp)[2]-1)){
    col <- tmp[[i]]
    N <- sum(col) # count number of mapped reads
    FPKMi <- (col*1e9)/(N*tmp[[dim(tmp)[2]]]) # count FPKM
    FPKMi <- pmax(FPKMi,0) %>% as.data.frame() # drop minus values
    colnames(FPKMi) <- colnames(tmp)[i]
    fpkm <- cbind(fpkm,FPKMi)
  }
  fpkm_forvalidation = fpkm # for validation, convert to TPM
  colnames(fpkm) = c("fpkm")
  fpkm$gene_id = rownames(fpkm)
  fpkm = fpkm[,c(2,1)] # put gene_id in front of fpkm


  # count2tpm
  tpm <- data.frame(row.names = tmp$gene_id)
  for (i in 2:(dim(tmp)[2]-1)){
    col <- tmp[[i]] %>% as.numeric()
    len <- tmp[[dim(tmp)[2]]] %>% as.numeric()
    rate <- col/len
    N <- sum(col) # count number of mapped reads
    TPMi <- (rate*1e6)/(sum(rate)) # count TPM
    #print(sum(rate))
    TPMi <- pmax(TPMi,0) %>% as.data.frame() # drop minus values
    colnames(TPMi) <- colnames(tmp)[i]
    tpm <- cbind(tpm,TPMi)
  }
  colnames(tpm) = c("tpm")
  tpm$gene_id = rownames(tpm)
  tpm = tpm[,c(2,1)] # put gene_id in front of tpm


  # fpkm2tpm - Optional, check if algorithms are correct by comparing
  # tpm1 and tpm 0
  #
  tpm1 <- data.frame(row.names =tmp$gene_id)
  FPKMtoTPM <- function(x) {
    return(exp(log(x) - log(sum(x)) + log(1e6)))
  }
  for (i in 1:(dim(fpkm_forvalidation)[2])){
    tpmi1 <- FPKMtoTPM(fpkm_forvalidation[[i]]) %>% as.data.frame()
    colnames(tpmi1) <- colnames(fpkm_forvalidation)[i]
    tpm1 <- cbind(tpm1,tpmi1)
  }


  # Bind FPKM and TPM dataframe
  expMatrMain = inner_join(fpkm, tpm, by="gene_id")
  write.table(expMatrMain, file = paste0(out_dir, "gene.csv"), row.names=FALSE, col.names=TRUE, sep=",")
}


te_count_normalizer <- function(te_count, te_lentable, out_dir){
  rm(list=ls())
  suppressMessages(library(dplyr))

  teCount = read.table(te_count, stringsAsFactors=T)
  teLength = read.table(te_lentable, stringsAsFactors=T)

  # Optional, remove header row
  #teCount = teCount[-1,]
  #teLength = teLength[-1,]

  colnames(teCount) = c("te_id","raw_read_count","te_class")
  colnames(teLength) = c("te_id", "length")

  teCount$raw_read_count = as.numeric(teCount$raw_read_count)
  teLength$length = as.numeric(teLength$length)

  # Remove duplicated feature (gene/TE) entries
  teCount = unique(teCount)
  teLength = unique(teLength)

  # Merge raw_counts df and length df into a large data frame
  tmp <- inner_join(teCount, teLength, by='te_id')

  # Remove individual TEs, leaving only TE families (TE_xxxx)
  tmp = tmp[grepl("^TE_", tmp$te_id),]
  # Remove duplicated TE families
  tmp = tmp %>%
    group_by(te_id) %>%
    distinct(te_id, .keep_all = T)
  # Re-order columns
  tmp_TEclass = tmp[,c(1,3)]
  tmp = tmp[,-3]
  tmp$raw_read_count = as.numeric(tmp$raw_read_count)
  tmp$length = as.numeric(tmp$length)

  # count2fpkm
  fpkm <- data.frame(row.names = tmp$te_id)
  for (i in 2:(dim(tmp)[2]-1)){
    col <- tmp[[i]]
    N <- sum(col)
    FPKMi <- (col*1e9)/(N*tmp[[dim(tmp)[2]]])
    FPKMi <- pmax(FPKMi,0) %>% as.data.frame()
    colnames(FPKMi) <- colnames(tmp)[i]
    fpkm <- cbind(fpkm,FPKMi)
  }
  fpkm_forvalidation = fpkm # for validation, convert to TPM
  colnames(fpkm) = c("fpkm")
  fpkm$te_id = rownames(fpkm)
  fpkm = fpkm[,c(2,1)] # put gene_id in front of fpkm


  # count2tpm
  tpm <- data.frame(row.names = tmp$te_id)
  for (i in 2:(dim(tmp)[2]-1)){
    col <- tmp[[i]]
    len <- tmp[[dim(tmp)[2]]]
    rate <- col/len
    N <- sum(col)
    TPMi <- (rate*1e6)/(sum(rate))
    #print(sum(rate))
    TPMi <- pmax(TPMi,0) %>% as.data.frame()
    colnames(TPMi) <- colnames(tmp)[i]
    tpm <- cbind(tpm,TPMi)
  }
  colnames(tpm) = c("tpm")
  tpm$te_id = rownames(tpm)
  tpm = tpm[,c(2,1)] # put gene_id in front of tpm


  # fpkm2tpm - Optional, check if algorithms are correct by comparing
  # tpm1 and tpm 0
  #
  tpm1 <- data.frame(row.names =tmp$te_id)
  FPKMtoTPM <- function(x) {
    return(exp(log(x) - log(sum(x)) + log(1e6)))
  }
  for (i in 1:(dim(fpkm_forvalidation)[2])){
    tpmi1 <- FPKMtoTPM(fpkm_forvalidation[[i]]) %>% as.data.frame()
    colnames(tpmi1) <- colnames(fpkm_forvalidation)[i]
    tpm1 <- cbind(tpm1,tpmi1)
  }


  # Bind FPKM and TPM dataframe
  expMatrMain = inner_join(fpkm, tpm, by="te_id")
  expMatrMain = inner_join(expMatrMain, tmp_TEclass, by="te_id")

  te_total = expMatrMain
  # Cataloging TE types
  te_total$te_class_group = "TE"

  DNA_TE_list = c("DNA/DTA", "DNA/DTC", "DNA/DTH", "DNA/DTM", "DNA/DTT", "DNA/Helitron")
  LTR_TE_list = c("LTR/Copia", "LTR/Gypsy", "LTR/unknown")
  MITE_TE_list = c("MITE/DTA", "MITE/DTC", "MITE/DTH", "MITE/DTM", "MITE/DTT")

  te_total[which(te_total$te_class %in% DNA_TE_list),5] = "DNA"
  te_total[which(te_total$te_class %in% LTR_TE_list),5] = "LTR"
  te_total[which(te_total$te_class %in% MITE_TE_list),5] = "MITE"

  te_total$te_class_group = as.factor(te_total$te_class_group)

  write.table(te_total, file = paste0(out_dir, "TE.csv"), row.names=FALSE, col.names=TRUE, sep=",")
}