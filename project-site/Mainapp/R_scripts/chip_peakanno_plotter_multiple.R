peak_analysis_multiple <- function (){

  rm(list=ls())

  suppressMessages(suppressWarnings(if (!require('ChIPseeker')) BiocManager::install('ChIPseeker',update = FALSE)));
  library(ChIPseeker)
  library(ggplot2)
  library(GenomicFeatures)
  library(ggupset)
  library(ggimage)
  library(ggplotify)
  library(viridis)
  library(Vennerable)

  files <- narrowPeak_file
  for (file in files){
    file_name <- gsub("([SE]RR\\d+).*", "\\1", file)
    file_name_list[[file_name]] <- file
  }
  tssRegion <- tssRegion
  species = species
  if (species == "sbi") {
    txdb <- loadDb("Mainapp/static/TxDb/TxDb.Sorghum_bicolor_NCBIv3.sqlite")
  } else if (species == "zma") {
    txdb <- loadDb("Mainapp/static/TxDb/TxDb.Zm-B73-REFERENCE-GRAMENE-4.0.sqlite")
  } else if (species == "ssp") {
    txdb <- loadDb("Mainapp/static/TxDb/TxDb.Saccharum_spontaneum.sqlite")
  }

  # annotation options
  ##Change the catagories
  options(ChIPseeker.ignore_1st_exon = ignore_1st_exon) # default F
  options(ChIPseeker.ignore_1st_intron = ignore_1st_intron)
  options(ChIPseeker.ignore_downstream = ignore_downstream)
  options(ChIPseeker.ignore_promoter_subcategory = ignore_promoter_subcategory)

  # PeaK Annotation
  peakAnnoList <- lapply(file_name_list, annotatePeak, TxDb=txdb,tssRegion=c(-tssRegion, tssRegion))
  df=as.data.frame(peakAnno)
  #save the dataframe
  write.csv(df,paste0("Mainapp/static/Temp_R_ChIP_projects/",output_folder_name,"/Peakannotation.csv"))

  # promoter region definition
  promoter <- getPromoters(TxDb = txdb, upstream = tssRegion, downstream = tssRegion)
  tagMatrixList <- lapply(file_name_list, getTagMatrix, windows=promoter)

  #plot1: Average Profile of ChIP peaks binding to TSS region for multiple samples, showing in desperate windows (long time)
  AvgPro = as.ggplot(plotAvgProf(tagMatrixList, xlim=c(-tssRegion, tssRegion), conf=0.95,resample=500, facet="row"))
  plt_AvgPro = AvgPro+
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
    )
  ggsave(paste0("Mainapp/static/Temp_R_ChIP_projects/", output_folder_name, "/allbindingplots.png"),
         plot=plt_AvgPro, dpi = 400, width=5, height=3)

  #plot2: The reads aligning heatmap in the TSS region for multiple samples
  peakHeatmap = as.ggplot(~tagHeatmap(tagMatrixList, xlim = c(-tssRegion, tssRegion), xlab = "Distance to TSS (bp)", color = "#00AFBB"))
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

  rm(list=ls())
}