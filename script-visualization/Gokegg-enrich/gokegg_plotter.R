filename_id = 11231
gokegg_dataframe  =read.csv2("./test_gokegg.csv", sep=",", header=T)[,-1]

gokegg_plotter <- function (gokegg_dataframe, filename_id){
  
  rm(list=ls())
  library(stringr)
  library(dplyr)
  library(clusterProfiler)
  library(tidyverse)
  library(reticulate)
  library(plotly)
  
  source_python("./filter_extra_GO_enrichment_lines.py")
  raw_gokegg_data = gokegg_dataframe
  output_file_name = filename_id
  
  print("Data & lib load complete.")
  # GO enrichment
  gene_ids = raw_gokegg_data$query_name
  eggnog_lines_with_go = raw_gokegg_data$GOs!= ""
  eggnog_annoations_go = str_split(raw_gokegg_data[eggnog_lines_with_go,]$GOs, ",")
  gene_to_go = data.frame(gene = rep(gene_ids[eggnog_lines_with_go],
                                     times = sapply(eggnog_annoations_go, length)),
                          term = unlist(eggnog_annoations_go))
  
  # Select genes for downstream analysis
  gene_list = gene_to_go$gene[which(gene_to_go$term != "-")]
  term2gene = gene_to_go[,c(2,1)]
  # Perform GO enricher
  df = enricher(gene=gene_list,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                TERM2GENE = term2gene)
  gene_df = as.data.frame(df@result)
  
  print("GO enriching completed.")
  # GO annotation
  # Remove lines without proper GO term description
  gene_df_GOannot = filter_extra_GOlines(
    gene_df, 
    go2term(gene_df$ID))
  colnames(gene_df_GOannot) = colnames(gene_df)
  
  # Add GO and ontology information to each GO entry (BP/CC, etc)
  # Sort to match go_term and ontology order (by ID)
  gene_df_GOannot_sorted = gene_df_GOannot[order(gene_df_GOannot$ID),]
  tmp_go_df = go2term(gene_df$ID)
  tmp_ont_df = go2ont(gene_df_GOannot$ID)
  gene_df_GOannot_sorted$GOterm = tmp_go_df$Term
  gene_df_GOannot_sorted$ontology = tmp_ont_df$Ontology
  
  # Sort to display most enriched GO terms
  gene_df_GOannot_final = gene_df_GOannot_sorted[order(gene_df_GOannot_sorted$p.adjust),]
  
  # Prepare data used for plotting, separate different ontology annotations
  data_pltgoenrich = gene_df_GOannot_final[1:20,]
  data_pltgoenrich = data_pltgoenrich[order(data_pltgoenrich$ontology),]
  data_pltgoenrich$Count = as.numeric(data_pltgoenrich$Count)
  #data_pltgoenrich$GeneRatio = as.numeric(data_pltgoenrich$GeneRatio)
  
  plt_GO_enrich = ggplot(data_pltgoenrich,
                         aes(x=GOterm,y=-log10(p.adjust))) +
    geom_point(aes(color=ontology, size=Count)) +
    coord_flip() +
    labs(x="") +
    facet_grid(ontology~., scale = 'free_y', space = 'free_y') +
    theme_bw()
  
  print("GO analysis completed.")
  ###
  ######
  ########## KEGG enrichment ############
  ######
  ###
  
  gene_ids = raw_gokegg_data$query_name
  
  eggnog_lines_with_ko = raw_gokegg_data$KEGG_ko!= ""
  eggnog_annoations_ko = str_split(raw_gokegg_data[eggnog_lines_with_ko,]$KEGG_ko, ",")
  gene_to_ko = data.frame(gene = rep(gene_ids[eggnog_lines_with_ko],
                                     times = sapply(eggnog_annoations_ko, length)),
                          ko = unlist(eggnog_annoations_ko)) %>%
    na.omit()
  
  gene_to_ko[,2] = gsub("ko:","",gene_to_ko[,2])
  
  # Data Cleaning
  gene_to_ko_cleaned = gene_to_ko[which(gene_to_ko$ko != "-"),]
  ko_to_gene = gene_to_ko_cleaned[,c(2,1)]
  
  
  ko_pathway_name = read.table("./KEGG_pathway_ko_uniq.txt", header=T, sep="\t")
  ko_to_enzyme = ko_pathway_name %>%
    dplyr::select(ko = ko, level3name = level3_pathway_name) %>%
    na.omit()
  ko_to_pathway = ko_pathway_name %>%
    dplyr::select(ko = ko, level2name = level2_pathway_name) %>%
    na.omit()
  
  pathway_to_gene = gene_to_ko %>% left_join(ko_to_enzyme, by = "ko") %>%
    dplyr::select(pathway=level3name,gene=gene) %>%
    na.omit()
  
  # Select genes for downstream analysis
  gene_list = gene_to_ko_cleaned$gene
  
  df_enrich_kegg = enricher(gene=gene_list,
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            TERM2GENE = pathway_to_gene,
                            # TERM2GENE = ko_to_gene,
                            # TERM2NAME = ko_to_pathway
  )
  print("KEGG enriching completed.")
  
  gene_df_kegg = as.data.frame(df_enrich_kegg@result)
  gene_df_kegg_cleaned = gene_df_kegg[!is.na(gene_df_kegg$Description),]
  
  # Sort to display most enriched GO terms
  #gene_df_kegg_final = gene_df_kegg_cleaned[rev(order(gene_df_kegg_cleaned$Count)),]
  gene_df_kegg_final = gene_df_kegg_cleaned[order(gene_df_kegg_cleaned$p.adjust),]
  
  # Prepare enrichment plot data
  data_pltkeggenrich = gene_df_kegg_final[1:20,]
  data_pltkeggenrich$Count = as.numeric(data_pltkeggenrich$Count)
  #data_pltkeggenrich$GeneRatio = as.numeric(data_pltkeggenrich$GeneRatio)
  
  
  plt_KEGG_enrich = ggplot(data_pltkeggenrich,
                           aes(x=Description,y=GeneRatio)) +
    geom_point(aes(color=-log10(p.adjust), size=Count)) +
    coord_flip() +
    labs(x="") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #plt_KEGG_enrich
  # df_enrich_kegg_significant = filter(df_enrich_kegg,
  #                                     pvalue<.05, qvalue<0.2)
  # 
  # plt_KEGG_dot = dotplot(df_enrich_kegg_significant,
  #                        showCategory=20,
  #                        color = "p.adjust") +
  #   #scale_size(range=c(2,12)) +
  #   scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE))
  
  #plt_KEGG_dot
  print("depicting KEGG plots completed.")
  
  ggsave(paste0("Mainapp/static/Temp_R_gokegg/", output_file_name, "_GOenrichment.png"), 
         plot=plt_GO_enrich, dpi = 400, width=10, height=10)
  output_go_enrichment = plotly::ggplotly(plt_GO_enrich)
  htmlwidgets::saveWidget(output_go_enrichment, 
                          paste0("Mainapp/static/Temp_R_gokegg/", output_file_name, "_GOenrichment.html"),
                          selfcontained = TRUE)
  print("Saved GO plot.")
  ggsave(paste0("Mainapp/static/Temp_R_gokegg/", output_file_name, "_KEGGenrichment.png"), 
         plot=plt_KEGG_enrich, dpi = 400, width=10, height=10)
  output_kegg_enrichment = plotly::ggplotly(plt_KEGG_enrich)
  htmlwidgets::saveWidget(output_kegg_enrichment, 
                          paste0("Mainapp/static/Temp_R_gokegg/", output_file_name, "_KEGGenrichment.html"),
                          selfcontained = TRUE)
  print("Saved KEGG plot 1.")
  # ggsave(paste0("Mainapp/static/Temp_R_gokegg/", output_file_name, "_KEGGenrichment2.png"), 
  #        plot=plt_KEGG_dot, dpi = 400, width=10, height=10)
  # output_kegg_enrichment2 = plotly::ggplotly(plt_KEGG_dot)
  # htmlwidgets::saveWidget(output_kegg_enrichment2, 
  #                         paste0("Mainapp/static/Temp_R_gokegg/", output_file_name, "_KEGGenrichment2.html"),
  #                         selfcontained = TRUE)
  # print("Saved KEGG plot 2.")
}

gokegg_plotter(gokegg_dataframe, filename_id)