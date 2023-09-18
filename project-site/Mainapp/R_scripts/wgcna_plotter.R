
gene_wgcna_step1 <- function (sampleGeneData,
                              RcCutoff, samplePerc, rscut, datatype, anamethod, dirname){
  # Data cleaning and sft calculation
  # ====================================
  # 'datatype''Data type: count or FPKM. count means readcount, FPKM represent normalized count,eg: FPKM TPM CPM RPKM,
  # 'anamethod''Data transformat method: For count: varianceStabilizingTransformation or cpm; For FPKM, rawFPKM or logFPKM,
  # 'RcCutoff''Noise cutoff: count/normalized count lower than ... was thought to be noise. Default: count=>10, Normalized count=>1',
  # 'samplePerc''Sample percentage: parameter for noise remove. xx percent of all samples have readcount/FPKM > cutoff. Default: 0.3',
  # 'GeneNum''Gene number for WGCNA: how many gene you want retained for WGCNA after noise remove.',
  # 'cutmethod''MAD or SVR, Default: MAD',
  # 'rscut''sft Power cutoff: Power cutoff',
  
  #setwd(paste0("global-temp/", dirname))
  rm(list=ls())
  
  # Check if dependencies are fully installed
  # suppressMessages(suppressWarnings(if (!require('logr')) install.packages('logr')));
  suppressMessages(suppressWarnings(if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")))
  suppressMessages(suppressWarnings(if (!require('devtools')) install.packages('devtools')));
  suppressMessages(suppressWarnings(if (!require('WGCNA')) BiocManager::install('WGCNA',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('tidyverse')) install.packages('tidyverse')));
  suppressMessages(suppressWarnings(if (!require('ShinyWGCNA')) devtools::install_github("ShawnWx2019/WGCNAShinyFun",ref = "master")));
  suppressMessages(suppressWarnings(if (!require('DESeq2')) BiocManager::install('DESeq2',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('ggprism')) BiocManager::install('ggprism',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('patchwork')) BiocManager::install('patchwork',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('ggtree')) install.packages('ggtree')));
  suppressMessages(library(ape));
  # Enable multi-threading
  enableWGCNAThreads(nThreads = 2)
  # For macOS
  #allowWGCNAThreads()
  
  # Testing parameters
  # sampleGeneData = "WGCNA_fpkm_data.csv"
  # RcCutoff = 1
  # samplePerc = 0.9
  # GeneNum = 2000
  # cutmethod = "MAD"
  # rscut = 0.8
  # datatype = "FPKM"
  # anamethod = "rawFPKM"

  # Conversion from character (python-relay) to number
  RcCutoff = as.numeric(RcCutoff)
  samplePerc = as.numeric(samplePerc)
  rscut = as.numeric(rscut)

  # Define function for noise filtering
  getdatExpr_local = function(rawdata,RcCutoff,samplePerc,datatype,method){
    ##> Noise filtering
    x <-
      rawdata %>%
      setNames(c("id",colnames(.)[-1])) %>%
      mutate(id = as.character(id)) %>% ## aviod numeric columns.
      column_to_rownames("id") %>%
      filter(
        rowSums(. > RcCutoff) > (samplePerc*ncol(.))
      )
    if (
      datatype == "expected count"
    ) {
      x <- x %>%
        mutate_if(is.numeric,ceiling)
    }
    ##> Normalization
    if(method == "vst") {
      dx <- x %>%
        as.matrix() %>%
        varianceStabilizingTransformation(., blind = TRUE)
    } else if (method == "logFPKM") {
      if (datatype == "normalized count") {
        dx <- log2(x+1)
      } else {
        dx <- log2(x+1)
      }
    } else if (method == "rawFPKM") {
      dx <- x
    } else {
      return()
    }
    return(dx)
  }

  datExpr = read.table(sampleGeneData, header=T, sep = ",") # rowname-genes, colname-samples
  ngenes_beforenorm = nrow(datExpr) # Obtain total number of genes
  write.table(paste0( "[Input] Genes: ", ngenes_beforenorm, ", Samples: ", ncol(datExpr)-1),
          file = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-log.txt"),
          append = T, na ="", row.names=FALSE, quote=FALSE, col.names = FALSE
          )
  write.table(paste0( "[Noise removal and gene filtration] Genes before filtering: ", ngenes_beforenorm),
          file = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-log.txt"),
          append = T, na ="", row.names=FALSE, quote=FALSE, col.names = FALSE
          )

  # Dataframe format conversion, as.character -> as.numeric
  vec_datexpr_rownames = datExpr$X
  datExpr = datExpr[,-1]
  datExpr = as.data.frame(lapply(datExpr,as.numeric))
  datExpr = cbind(vec_datexpr_rownames, datExpr)
  #rownames(dx) = vec_datexpr_rownames

  # --NOISE REMOVAL--
  datExpr <- getdatExpr_local(datExpr, RcCutoff, samplePerc, datatype, anamethod)

  # --GENE FILTRATION--
  datExpr <- as.data.frame(t(datExpr))
  gsg = goodSamplesGenes(datExpr, verbose = 3)
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:",
                       paste(names(datExpr)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:",
                       paste(rownames(datExpr)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
  }
  ngenes_afternorm = ncol(datExpr) # Obtain total number of genes
  write.table(paste0( "[Noise removal and gene filtration] Genes after filtering: ", ngenes_afternorm),
          file = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-log.txt"),
          append = T, na ="", row.names=FALSE, quote=FALSE, col.names = FALSE
          )

  # Generate sample cluster dendrogram based on expression values
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  sampleTree = hclust(dist(datExpr), method = "average")
  treenew= as.phylo(sampleTree)
  anno = data.frame(lab = treenew$tip.label,
                    group = factor(gsub(".$","",treenew$tip.label),levels = unique(gsub(".$","",treenew$tip.label))))
  phylotree =ggtree(tr = treenew, layout = 'dendrogram', size=1)
  phylotree2 = phylotree %<+% anno +theme(legend.position = "") +
    geom_tiplab(hjust=0, vjust=1.5, angle=90,offset = -70) +
    geom_nodepoint(color="#a6e3e9", alpha=0.5, size=10) +
    geom_nodelab(hjust=-0.2, size=4, color="#b83b5e") + # Bootstrap values
    theme_tree()
  ggsave(paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-sample-dendro.png"), plot=phylotree2, dpi = 400)

  # --SFT PREDICTION--
  step1_sft = getpower(datExpr = datExpr, rscut = rscut)
  # Save scale-independence and mean connectivty plot
  png(filename = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-sft.png"), width=800, height=400, res=100)
  print(step1_sft$plot)
  dev.off()

  pdf(paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-sft.pdf"), width=8, height=4)
  print(step1_sft$plot)
  dev.off()
  # savePlot(filename = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-sft.png"),
  #          type = "png", device = dev.cur())

  write_csv(step1_sft$sft, file = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-sft-data.csv"))
  write.table(paste0( "[SFT prediction]: Recommended SFT soft power is: ", step1_sft$power),
            file = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-log.txt"),
            append = T, na ="", row.names=FALSE, quote=FALSE, col.names = FALSE
            )

  # Save workspace
  save(datExpr, RcCutoff, samplePerc, rscut, datatype, anamethod,
       file=paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-vars.RData"))
  write.table("[Output] All done, please proceed to section 2 of the analysis.",
          file = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-log.txt"),
          append = T, na ="", row.names=FALSE, quote=FALSE, col.names = FALSE
          )
}

gene_wgcna_step2 <- function (sftPower, minModuleSize, mergeCutHeight, dirname){
  # Gene expression network construction
  # ====================================
  # 'sftPower''SFT power.',
  # 'minModuleSize''Minimal module size: The gene number of minimal module. Default: 30',
  # 'mergeCutHeight''Merge cuttree height: tree height lower than this value will be merged. Default: 0.25',
  
  #setwd(paste0("global-temp/", dirname))
  rm(list=ls())
  
  # Check if dependencies are fully installed
  suppressMessages(suppressWarnings(library(WGCNA)));
  suppressMessages(suppressWarnings(library(ShinyWGCNA)));
  suppressMessages(suppressWarnings(library(tidyverse)));
  # Enable multi-threading
  enableWGCNAThreads(nThreads = 2)
  
  # Testing parameters
  # sftPower = 6
  # minModuleSize = 30
  # mergeCutHeight = 0.25

  # Conversion from character (python-relay) to number
  sftPower = as.numeric(sftPower)
  minModuleSize = as.numeric(minModuleSize)
  mergeCutHeight = as.numeric(mergeCutHeight)
  
  # Load data from previous steps
  load(paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step1-vars.RData"))
  
  # Set max block size
  maxBlocksize = ncol(datExpr)
  # 'maxBlocksize''max block size: For block-wised network construction method, the block size, which is set to be all filtered genes.'


  step2_network = getnetwork(datExpr = datExpr, power = sftPower,
                             minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight, maxBlocksize = maxBlocksize)
  # Visualization using the hierarchical clustering dendrogram
  gene_tree = step2_network$net$dendrograms[[1]]
  png(filename = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step2-network-dendro.png"), width=1000, height=600, res=100)
  plotDendroAndColors(
    gene_tree,
    step2_network$moduleColors[step2_network$net$blockGenes[[1]]],
    dendroLabels = F
  )
  dev.off()

  pdf(paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step2-network-dendro.pdf"), width=13, height=6)
  plotDendroAndColors(
    gene_tree,
    step2_network$moduleColors[step2_network$net$blockGenes[[1]]],
    dendroLabels = F
  )
  dev.off()

  ME_df = step2_network$MEs_col
  ME_dfnew = cbind(rownames(ME_df), ME_df)
  write_csv(ME_dfnew, file = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step2-module-eigengenes.csv"))
  write_csv(step2_network$Gene2module, file = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step2-module-allgenes.csv"))

  # Save workspace
  #save.image("step2.RData")
  save(step2_network, datExpr,
       RcCutoff, samplePerc, rscut, datatype, anamethod,
       sftPower, minModuleSize, mergeCutHeight, maxBlocksize,
       file=paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step2-vars.RData"))

}

gene_wgcna_step3 <- function (dirname){
  # Module-trait correlation analysis
  # ====================================
  # 'traitData''trait-sample information, trait represented with 0 or 1',
  
  #setwd(paste0("global-temp/", dirname))

  suppressMessages(suppressWarnings(library(ggplot2)));
  suppressMessages(suppressWarnings(library(WGCNA)));
  suppressMessages(suppressWarnings(library(ShinyWGCNA)));
  suppressMessages(suppressWarnings(library(tidyverse)));
  suppressMessages(suppressWarnings(library(ggprism)));
  suppressMessages(suppressWarnings(library(patchwork)));
  suppressMessages(suppressWarnings(library(ComplexHeatmap)));

  # Enable multi-threading
  enableWGCNAThreads(nThreads = 2)
  
  # Testing parameters
  # traitData = "WGCNA_data_trait.csv"
  
  # Load data from previous steps
  load(paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step2-vars.RData"))
  # set options for scientific number notation
  options(scipen = 6)
  
  #log_print("--MODULE-TRAIT CORRELATION ANALYSIS STARTED--")
  phen <- read.csv(paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/sample_trait_data.csv"))
  
  # Bind sample-id column to the dataframe index
  colnames(phen)[1] = "sample_id"
  phen <- left_join(
    data.frame(
      sample_id = rownames(datExpr)
    ),phen,"sample_id"
  ) %>% 
    column_to_rownames("sample_id")

  step3_mdule_trait <- getMt(
    phenotype = phen,
    nSamples = nrow(datExpr),
    moduleColors = step2_network$moduleColors,
    datExpr = datExpr
  )

  
  mod_color = gsub(pattern = "^..",replacement = "",rownames(step3_mdule_trait$modTraitCor))
  mod_color_anno = setNames(mod_color,rownames(step3_mdule_trait$modTraitCor))
  
  left_anno = rowAnnotation(
    Module = rownames(step3_mdule_trait$modTraitCor),
    col = list(Module = mod_color_anno),
    show_legend = F,
    show_annotation_name = F
  )
  
  # Generate data matrix, generate heatmap
  htmap_df = as.matrix(step3_mdule_trait$modTraitCor)
  
  # Save heatmap to static image
  png(filename = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step3-module-trait-heatmap.png"), width=2000, height=2000, res=200)
  print(pheatmap(htmap_df,
           display_numbers = TRUE,
           color = colorRampPalette(c("#00AFBB", "white", "#c93756"))(256),
           fontsize=10, border_color = "grey60",
           main = "Module-trait heatmap",
           treeheight_row=50, treeheight_col = 30,
           cellwidth = 50, cellheight = 40
  ))
  dev.off()

  pdf(paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step3-module-trait-heatmap.pdf"), width=20, height=20)
  print(pheatmap(htmap_df,
           display_numbers = TRUE,
           color = colorRampPalette(c("#00AFBB", "white", "#c93756"))(256),
           fontsize=10, border_color = "grey60",
           main = "Module-trait heatmap",
           treeheight_row=50, treeheight_col = 30,
           cellwidth = 50, cellheight = 40
  ))
  dev.off()

  htmap_dfnew = as.data.frame(htmap_df)
  htmap_dfnew = cbind(rownames(htmap_dfnew), htmap_dfnew)
  write_csv(as.data.frame(htmap_dfnew), file = paste0("Mainapp/static/Temp_R_wgcna/", dirname, "/step3-module-trait-correlation.csv"))

}
