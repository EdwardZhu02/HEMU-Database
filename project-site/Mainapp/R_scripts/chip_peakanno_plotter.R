peak_analysis_single <- function (narrowPeak_file, tssRegion, species,
                                  ignore_1st_exon, ignore_1st_intron, ignore_downstream, ignore_promoter_subcategory,
                                  output_folder_name){

  rm(list=ls())

  suppressMessages(suppressWarnings(if (!require('ChIPseeker')) BiocManager::install('ChIPseeker',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('GenomicFeatures')) BiocManager::install('GenomicFeatures',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('ggupset')) BiocManager::install('ggupset',update = FALSE)));
  #suppressMessages(suppressWarnings(if (!require('ggplotify')) BiocManager::install('ggplotify',update = FALSE)));
  suppressMessages(suppressWarnings(if (!require('ggimage')) install.packages('ggimage')));
  suppressMessages(suppressWarnings(if (!require('viridis')) install.packages('viridis')));
  suppressMessages(suppressWarnings(if (!require('ggplot2')) install.packages('ggplot2')));

  #print(paste(narrowPeak_file, tssRegion, species,
  #      ignore_1st_exon, ignore_1st_intron, ignore_downstream, ignore_promoter_subcategory,
  #      output_folder_name, sep=" "))

  f <- narrowPeak_file
  # file_name <- gsub("([SE]RR\\d+).*", "\\1", f)
  tssRegion <- as.numeric(tssRegion)
  minus_tssRegion <<- 0 - as.numeric(tssRegion) # Define global variable

  species = species
  if (species == "sorghum") {
    txdb <- loadDb("global-static/TxDb/TxDb.Sorghum_bicolor_NCBIv3.sqlite")
  } else if (species == "zea") {
    txdb <- loadDb("global-static/TxDb/TxDb.Zm-B73-REFERENCE-GRAMENE-4.0.sqlite")
  } else if (species == "saccharum") {
    txdb <- loadDb("global-static/TxDb/TxDb.Saccharum_spontaneum.sqlite")
  }
  
  # set additional parameters
  if (ignore_1st_exon == "T"){ options(ChIPseeker.ignore_1st_exon = T)  }
      else{ options(ChIPseeker.ignore_1st_exon = F) }

  if (ignore_1st_intron == "T"){ options(ChIPseeker.ignore_1st_intron = T) }
      else { options(ChIPseeker.ignore_1st_intron = F) }

  if (ignore_downstream == "T") { options(ChIPseeker.ignore_downstream = T) }
    else { options(ChIPseeker.ignore_downstream = F) }

  if (ignore_promoter_subcategory == "T") { options(ChIPseeker.ignore_promoter_subcategory = T) }
    else { options(ChIPseeker.ignore_promoter_subcategory = F) }


  # PeaK Annotation
  peakAnno <<- annotatePeak(f, tssRegion=c(minus_tssRegion, tssRegion), TxDb=txdb) # Define global variable
  df = as.data.frame(peakAnno)
  #save the dataframe
  write.csv(df,paste0("Mainapp/static/Temp_R_epigenome/",output_folder_name,"/Peakannotation.csv"))
  
  # promoter region definition
  promoter <- getPromoters(TxDb = txdb, upstream = tssRegion, downstream = tssRegion)
  tag_matrix <<- getTagMatrix(f, windows = promoter) # Define global variable

  #plot1: The reads aligning in the TSS region
  print("Generating plot1")
  #peakHeatmap = as.ggplot(~tagHeatmap(tag_matrix, xlim = c(minus_tssRegion, tssRegion), xlab = "Distance to TSS (bp)", color = "#00AFBB"))
  #plt_peakHeatmap = peakHeatmap+
  #  theme(panel.grid = element_blank(),
  #        panel.background = element_blank(),
  #  )
  #ggsave(paste0("Mainapp/static/Temp_R_epigenome/", output_folder_name, "/Alignplot.png"),
  #       plot=plt_peakHeatmap, dpi = 400, width=5, height=7)
  png(filename = paste0("Mainapp/static/Temp_R_epigenome/", output_folder_name, "/Alignplot.png"), width=2000, height=2800, res=400)
  tagHeatmap(tag_matrix, xlim = c(minus_tssRegion, tssRegion), xlab = "Distance to TSS (bp)", color = "#00AFBB")
  dev.off()

  print("Generating plot2")
  #plot2: Average Profile of ChIP peaks binding to TSS region
  bind_plotname = plotAvgProf(tag_matrix, xlim=c(minus_tssRegion, tssRegion),
              xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency",conf = 0.95)
  plt_bind_plot = bind_plotname +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
    )
  ggsave(paste0("Mainapp/static/Temp_R_epigenome/", output_folder_name, "/Bindplot.png"),
         plot=plt_bind_plot, dpi = 400, width=10, height=6)
  output_plt_bind_plot <- plotly::ggplotly(plt_bind_plot)
  htmlwidgets::saveWidget(output_plt_bind_plot, paste0("Mainapp/static/Temp_R_epigenome/", output_folder_name, "/Bindplot.html"),
                          selfcontained = TRUE)

  print("Generating plot3-pie")
  #plot3 Pie plot
  #region_pieplot <- as.ggplot(~plotAnnoPie(peakAnno))
  #plt_region_pie_plot = region_pieplot+
  #  theme(panel.grid = element_blank(),
  #        panel.background = element_blank(),
  #  )
  #ggsave(paste0("Mainapp/static/Temp_R_epigenome/", output_folder_name, "/Region_pieplot.png"),
  #       plot=plt_region_pie_plot, dpi = 400, width=10, height=6)
  png(filename = paste0("Mainapp/static/Temp_R_epigenome/", output_folder_name, "/Region_pieplot.png"), width=4000, height=2400, res=400)
  plotAnnoPie(peakAnno)
  dev.off()

  #print("Generating plot3-venn")
  #plot3: Pie plot combining Venn plot for the peak regions
  #region_vennpieplot <- upsetplot(peakAnno, vennpie=TRUE)
  #plt_region_vennpie_plot = region_vennpieplot+
  #  theme(panel.grid = element_blank(),
  #        panel.background = element_blank(),
  #  )
  #ggsave(paste0("Mainapp/static/Temp_R_epigenome/", output_folder_name, "/Region_vennpieplot.png"),
  #       plot=plt_region_vennpie_plot, dpi = 400, width=15, height=9)
  png(filename = paste0("Mainapp/static/Temp_R_epigenome/", output_folder_name, "/Region_vennpieplot.png"), width=6000, height=3600, res=400)
  print(upsetplot(peakAnno, vennpie=TRUE))
  dev.off()
  
  print("Generating plot4")
  #plot4: Distance To TSS
  plotDistToTSS <- plotDistToTSS(peakAnno)
  plt_DistToTSS_plot = plotDistToTSS +
    theme(panel.grid = element_blank(),
           panel.background = element_blank(),
           axis.line = element_line(colour = 'black'),
           axis.text.x = element_text(size = 10, hjust = 1),
           axis.text.y = element_text(size = 10, hjust = 1),
           plot.title = element_text(hjust = 0.5),
           )+
    labs(title="Distribution of TF binding loci relative to TSS (predicted)")
  ggsave(paste0("Mainapp/static/Temp_R_epigenome/", output_folder_name, "/DistToTSS.png"),
         plot=plt_DistToTSS_plot, dpi = 400, width=15, height=10)
  output_DistToTSS_plot <- plotly::ggplotly(plt_DistToTSS_plot)
  htmlwidgets::saveWidget(output_DistToTSS_plot, paste0("Mainapp/static/Temp_R_epigenome/", output_folder_name, "/DistToTSS.html"),
                          selfcontained = TRUE)

  rm(list=ls())
}