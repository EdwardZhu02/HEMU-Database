go_enricher <- function (gokegg_dataframe, filename_id, bg_annot){
    suppressMessages(library(stringr))
    suppressMessages(library(dplyr))
    suppressMessages(library(clusterProfiler))
    suppressMessages(library(tidyverse))
    suppressMessages(library(plotly))

    raw_gokegg_data = gokegg_dataframe  # Dataframe generaated from curator
    # ['query_name', 'Description', 'GOs', 'KEGG_ko', 'KEGG_Pathway']

    output_file_name = filename_id  # Unique task filename identifier
    # str(random.randint(int(1e7), int(1e8)))

    ###
    ######
    ########## Preparation of background annotation ############
    ######
    ###
    background_gokegg_data = read.table(bg_annot,sep="\t",header=F, fill = T)
    background_gokegg_data[background_gokegg_data==""]<-NA
    colnames(background_gokegg_data) = c("query_name", "seed_ortholog", "evalue", "score",	"eggNOG_OGs",	"max_annot_lvl",	"COG_category",	"Description",	"Preferred_name",	"GOs",	"EC",	"KEGG_ko",	"KEGG_Pathway",	"KEGG_Module",	"KEGG_Reaction",	"KEGG_rclass",	"BRITE",	"KEGG_TC",	"CAZy",	"BiGG_Reaction",	"PFAMs")
    background_gene_ids = background_gokegg_data$query_name

    eggnog_lines_with_go = background_gokegg_data$GOs!= ""  # Filter background annotation lines with a GO annotation
    eggnog_annoations_go = str_split(background_gokegg_data[eggnog_lines_with_go,]$GOs, ",")
    background_gene_to_go = data.frame(gene = rep(background_gene_ids[eggnog_lines_with_go],
                                       times = sapply(eggnog_annoations_go, length)),
                                       term = unlist(eggnog_annoations_go))
    background_go_to_gene = background_gene_to_go[,c(2,1)]
    background_go_to_gene = background_go_to_gene[which(background_go_to_gene$term != "-"),]
    background_go_to_gene$gene<-gsub("transcript:","",background_go_to_gene$gene)
    # Result: 'background_go_to_gene' dataframe

    ###
    ######
    ########## GO enrich ############
    ######
    ###
    gene_list = raw_gokegg_data$query_name

    # Perform GO enricher
    df_enrich_go<-enricher(gene=gene_list,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 TERM2GENE = background_go_to_gene)
    gene_df = as.data.frame(df_enrich_go@result)

    # Remove lines without proper GO term description
    go2term_ID = go2term(gene_df$ID)
    rownames(go2term_ID) = go2term_ID$go_id
    gene_df_filtered = gene_df[-which(!gene_df$ID %in% rownames(go2term_ID)),]
    gene_df_GOannot = gene_df_filtered[order(gene_df_filtered$ID),]


    # Add GO and ontology information to each GO entry (BP/CC, etc)
    # Sort to match go_term and ontology order (by ID)
    gene_df_GOannot_sorted = gene_df_GOannot[order(gene_df_GOannot$ID),]
    tmp_go_df = go2term(gene_df$ID)
    tmp_ont_df = go2ont(gene_df_GOannot$ID)
    gene_df_GOannot_sorted$GOterm = tmp_go_df$Term
    gene_df_GOannot_sorted$ontology = tmp_ont_df$Ontology

    # Sort to display most enriched GO terms (top 20)
    gene_df_GOannot_final = gene_df_GOannot_sorted[order(gene_df_GOannot_sorted$pvalue),]
    if (nrow(gene_df_GOannot_final) > 20) {
        data_pltgoenrich = gene_df_GOannot_final[1:20,]
    } else {
        data_pltgoenrich = gene_df_GOannot_final
    }


    # Prepare data used for plotting, separate different ontology annotations
    data_pltgoenrich = data_pltgoenrich[order(data_pltgoenrich$ontology),]
    data_pltgoenrich$Count = as.numeric(data_pltgoenrich$Count)
    data_pltgoenrich$pvalue = as.numeric(data_pltgoenrich$pvalue)
    # Extract numeric gene ratio
    data_pltgoenrich$Gene_Ratio = 0
    for (i in 1:nrow(data_pltgoenrich)){
      tmp_linesplit = unlist(strsplit(data_pltgoenrich[i,3], "/")) # GeneRatio
      data_pltgoenrich[i,ncol(data_pltgoenrich)] = as.numeric(tmp_linesplit[1]) / as.numeric(tmp_linesplit[2])
    }

    # Output enrichment result into a single .csv file.
    write.csv(gene_df_GOannot_final, paste0("Mainapp/static/Temp_R_gokegg/", output_file_name, "_enrichment.csv"), row.names = F)

    #plt_GO_enrich = ggplot(data_pltgoenrich, aes(x=GOterm,y=-log10(p.adjust))) +
    plt_GO_enrich = ggplot(data_pltgoenrich, aes(x=GOterm,y=Gene_Ratio)) +
      geom_point(aes(color=pvalue, size=Count)) +
      coord_flip() +
      labs(x="") +
      facet_grid(ontology~., scale = 'free_y', space = 'free_y') +
      theme_bw() +
      #theme(legend.position = 'none') +
      #scale_x_discrete(labels=function(x) str_wrap(x, width=40))
      scale_x_discrete(labels = function(x)  str_sub(x, 1, 70))


    ggsave(paste0("Mainapp/static/Temp_R_gokegg/", output_file_name, "_enrichment.png"),
           plot=plt_GO_enrich, dpi = 400, width=10, height=10)
    output_go_enrichment = plotly::ggplotly(plt_GO_enrich)
    htmlwidgets::saveWidget(output_go_enrichment,
                            paste0("Mainapp/static/Temp_R_gokegg/", output_file_name, "_enrichment.html"),
                            selfcontained = TRUE)

    #df_enrich_go_significant = filter(df_enrich_go, pvalue<0.05)
    #plt_GO_cnet = cnetplot(df_enrich_go,categorySize="p.adjust")
    #ggsave(paste0("Mainapp/static/Temp_R_gokegg/", output_file_name, "_GOcnet.png"),
               #plot=plt_GO_cnet, dpi = 150, width=10, height=10)
    rm(list=ls())
}


kegg_enricher <- function (gokegg_dataframe, filename_id, bg_annot){
    library(stringr)
    library(dplyr)
    library(clusterProfiler)
    library(tidyverse)
    library(plotly)

    raw_gokegg_data = gokegg_dataframe  # Dataframe generaated from curator
    # ['query_name', 'Description', 'GOs', 'KEGG_ko', 'KEGG_Pathway']

    output_file_name = filename_id  # Unique task filename identifier
    # str(random.randint(int(1e7), int(1e8)))

    ###
    ######
    ########## Preparation of background annotation ############
    ######
    ###
    background_gokegg_data = read.table(bg_annot,sep="\t",header=F, fill = T)
    background_gokegg_data[background_gokegg_data==""]<-NA
    colnames(background_gokegg_data) = c("query_name", "seed_ortholog", "evalue", "score",	"eggNOG_OGs",	"max_annot_lvl",	"COG_category",	"Description",	"Preferred_name",	"GOs",	"EC",	"KEGG_ko",	"KEGG_Pathway",	"KEGG_Module",	"KEGG_Reaction",	"KEGG_rclass",	"BRITE",	"KEGG_TC",	"CAZy",	"BiGG_Reaction",	"PFAMs")
    background_gene_ids = background_gokegg_data$query_name

    eggnog_lines_with_ko <- background_gokegg_data$KEGG_ko != ""
    eggnog_annoations_ko <- str_split(background_gokegg_data[eggnog_lines_with_ko,]$KEGG_ko, ",")
    background_gene_to_ko <- data.frame(gene = rep(background_gene_ids[eggnog_lines_with_ko],
                                        times = sapply(eggnog_annoations_ko, length)),
                                        ko = unlist(eggnog_annoations_ko)) %>% na.omit()

    background_gene_to_ko[,2]<-gsub("ko:","",background_gene_to_ko[,2])

    # Data transformation
    #gene_to_ko_cleaned = background_gene_to_ko[which(background_gene_to_ko$ko != "-"),]
    #ko_to_gene = gene_to_ko_cleaned[,c(2,1)]

    ko_pathway_name = read.table("global-static/gokegg-bgannot/KEGG_pathway_ko_uniq.txt", header=T, sep="\t")
    ko_to_enzyme = ko_pathway_name %>% dplyr::select(ko = ko, level3name = level3_pathway_name) %>% na.omit()

    #ko_to_pathway = ko_pathway_name %>%
    #   dplyr::select(ko = ko, level2name = level2_pathway_name) %>%
    #   na.omit()

    pathway_to_gene <- background_gene_to_ko %>% left_join(ko_to_enzyme, by = "ko") %>%
      dplyr::select(pathway=level3name,gene=gene) %>%
      na.omit()
    pathway_to_gene$gene<-gsub("transcript:","",pathway_to_gene$gene)
    # Result: 'pathway_to_gene' dataframe

    ###
    ######
    ########## KEGG enrich ############
    ######
    ###
    gene_list = raw_gokegg_data$query_name

    df_enrich_kegg = enricher(gene=gene_list,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 TERM2GENE = pathway_to_gene,
                 # TERM2GENE = ko_to_gene,
                 # TERM2NAME = ko_to_pathway
                 )

    df_enrich_kegg_result = df_enrich_kegg@result
    df_enrich_kegg_result$Gene_Ratio = 0
    df_enrich_kegg_result$pvalue = as.numeric(df_enrich_kegg_result$pvalue)

    for (i in seq_len(nrow(df_enrich_kegg_result))){
      tmp_linesplit = unlist(strsplit(df_enrich_kegg_result[i,3], "/")) # GeneRatio
      df_enrich_kegg_result[i,ncol(df_enrich_kegg_result)] = as.numeric(tmp_linesplit[1]) / as.numeric(tmp_linesplit[2])
    }

    # Sort to display most enriched KEGG terms (top 20)
    df_enrich_kegg_result = df_enrich_kegg_result[order(df_enrich_kegg_result$pvalue),]
    if (nrow(df_enrich_kegg_result) > 20) {
        df_enrich_kegg_result = df_enrich_kegg_result[1:20,]


    # Output enrichment result into a single .csv file.
    write.csv(df_enrich_kegg_result, paste0("Mainapp/static/Temp_R_gokegg/", output_file_name, "_enrichment.csv"), row.names = F)

    plt_KEGG_enrich = ggplot(df_enrich_kegg_result, aes(x=Description,y=Gene_Ratio)) +
      geom_point(aes(color=pvalue, size=Count)) +
      coord_flip() +
      labs(x="") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)
      )

    ggsave(paste0("Mainapp/static/Temp_R_gokegg/", output_file_name, "_enrichment.png"),
               plot=plt_KEGG_enrich, dpi = 400, width=10, height=10)
    output_kegg_enrichment = plotly::ggplotly(plt_KEGG_enrich)
    htmlwidgets::saveWidget(output_kegg_enrichment,
                            paste0("Mainapp/static/Temp_R_gokegg/", output_file_name, "_enrichment.html"),
                            selfcontained = TRUE)
    rm(list=ls())
}