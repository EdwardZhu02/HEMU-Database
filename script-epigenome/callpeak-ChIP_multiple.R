Peak_analysis_multiple <- function (){

  suppressMessages(suppressWarnings(if (!require('ChIPseeker')) BiocManager::install('ChIPseeker',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('GenomicFeatures')) BiocManager::install('GenomicFeatures',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('ggupset')) BiocManager::install('ggupset',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('ggplotify')) BiocManager::install('ggplotify',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('ggimage')) install.packages('ggimage')));
  suppressMessages(suppressWarnings(if (!require('viridis')) install.packages('viridis')));
  suppressMessages(suppressWarnings(if (!require('ggplot2')) install.packages('ggplot2')));
  
  suppressMessages(suppressWarnings(if (!require('graph')) BiocManager::install('graph',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('RBGL')) BiocManager::install('RBGL',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('Vennerable')) devtools::install_github("js229/Vennerable")));
  
  # Diag parameters
  narrowPeak_file = c("Mainapp/static/ChIP/SRR5748777V4.bam_deduplicate_bam_peaks.narrowPeak",
                      "Mainapp/static/ChIP/SRR5748778V4.bam_deduplicate_bam_peaks.narrowPeak",
                      "Mainapp/static/ChIP/SRR5748779V4.bam_deduplicate_bam_peaks.narrowPeak"
                      )
  tssRegion = "3000"
  species = "sorghum"
  output_folder_name = "test_multiple"
  options(ChIPseeker.ignore_1st_exon = F) # default F
  options(ChIPseeker.ignore_1st_intron = F)
  options(ChIPseeker.ignore_downstream = F)
  options(ChIPseeker.ignore_promoter_subcategory = F)
  # Diag parameters end
  
  files <- narrowPeak_file
  file_name_list <- list()
  for (file in files){
    file_name <- gsub(".*([SE]RR\\d+).*", "\\1", file)
    file_name_list[[file_name]] <- file
  }
  
  #tssRegion <- tssRegion
  #species = species
  tssRegion <- as.numeric(tssRegion)
  minus_tssRegion <<- 0 - as.numeric(tssRegion) # Define global variable
  
  if (species == "sorghum") {
    txdb <- loadDb("global-static/TxDb/TxDb.Sorghum_bicolor_NCBIv3.sqlite")
  } else if (species == "zea") {
    txdb <- loadDb("global-static/TxDb/TxDb.Zm-B73-REFERENCE-GRAMENE-4.0.sqlite")
  } else if (species == "saccharum") {
    txdb <- loadDb("global-static/TxDb/TxDb.Saccharum_spontaneum.sqlite")
  }
  
  # annotation options
  ##Change the catagories
  # options(ChIPseeker.ignore_1st_exon = ignore_1st_exon) # default F
  # options(ChIPseeker.ignore_1st_intron = ignore_1st_intron)
  # options(ChIPseeker.ignore_downstream = ignore_downstream)
  # options(ChIPseeker.ignore_promoter_subcategory = ignore_promoter_subcategory)
  
  # PeaK Annotation
  peakAnnoList <- lapply(file_name_list, annotatePeak, TxDb=txdb,tssRegion=c(minus_tssRegion, tssRegion))
  
  # Save the peak annotation to individual dataframes
  for(peakname in names(peakAnnoList)){
    df=as.data.frame(peakAnnoList[[peakname]])
    write.csv(df,paste0("Mainapp/static/Temp_R_ChIP_projects/", output_folder_name, "/", peakname, "_Peakannotation.csv"))
  }
  
  # promoter region definition
  promoter <- getPromoters(TxDb = txdb, upstream = tssRegion, downstream = tssRegion)
  tagMatrixList <- lapply(file_name_list, getTagMatrix, windows=promoter)
  
  #plot1: Average Profile of ChIP peaks binding to TSS region for multiple samples, showing in desperate windows (long time)
  AvgPro = as.ggplot(plotAvgProf(tagMatrixList, xlim=c(minus_tssRegion, tssRegion), conf=0.95,resample=500, facet="row"))
  plt_AvgPro = AvgPro+
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
    )
  ggsave(paste0("Mainapp/static/Temp_R_ChIP_projects/", output_folder_name, "/allbindingplots.png"),
         plot=plt_AvgPro, dpi = 400, width=5, height=3)
  
  #plot2: The reads aligning heatmap in the TSS region for multiple samples
  peakHeatmap = as.ggplot(~tagHeatmap(tagMatrixList, xlim = c(minus_tssRegion, tssRegion), xlab = "Distance to TSS (bp)", color = "#00AFBB"))
  plt_peakHeatmap = peakHeatmap+
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
    )
  ggsave(paste0("Mainapp/static/Temp_R_ChIP_projects/", output_folder_name, "/Alignplot.png"),
         plot=plt_peakHeatmap, dpi = 400, width=5, height=3)
  
  #plot3: Bar plot for Distance To TSS among multiple samples
  TssPlot = as.ggplot(plotDistToTSS(peakAnnoList))
  plt_TssPlot = TssPlot+
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
    )
  ggsave(paste0("Mainapp/static/Temp_R_ChIP_projects/", output_folder_name, "/TssPlot.png"),
         plot=plt_peakHeatmap, dpi = 400, width=5, height=3)
  
  #plot4: Bar plot for the peak regions among multiple samples
  peakAnnoList = as.ggplot(plotAnnoBar(peakAnnoList))
  plt_peakAnnoList = peakAnnoList+
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
    )
  ggsave(paste0("Mainapp/static/Temp_R_ChIP_projects/", output_folder_name, "/peakAnnoList.png"),
         plot=plt_peakAnnoList, dpi = 400, width=5, height=3)
  
  
  #plot5: Venn plot for the relating genes
  genes <- lapply(peakAnnoList, function(i) 
    as.data.frame(i)$geneId)
  data <- Vennerable::Venn(genes)
  VennPlot = as.ggplot(plot(data,doWeight=F))
  plt_VennPlot = VennPlot+
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
    )
  ggsave(paste0("Mainapp/static/Temp_R_ChIP_projects/", output_folder_name, "/Vennplot.png"),
         plot=plt_VennPlot, dpi = 400, width=5, height=3)
}
