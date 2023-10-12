
exp_heatmap_plotter <- function(wkdir_name, exp_df, query_format){
  # Enforcing multiple sequence alignment
  # Params:
  # wkdir_name: "htmap" + 8-digit random number, ex,. htmap12345678
  # exp_df: dataframe with col names ["gene_id", "sample_id", "fpkm", "tissue"] or ["gene_id", "sample_id", "tpm", "tissue"]
  # query_format: FPKM / TPM, used for plot labeling
  
  suppressMessages(library(tidyverse))
  suppressMessages(library(ggtree)) # sample clustering
  suppressMessages(library(reshape2)) # dataframe reformatting, melt/cast
  suppressMessages(library(aplot)) # concatenating plots
  suppressMessages(library(RColorBrewer)) # concatenating plots
  rm(list=ls())
  
  # Diagnosis, read matrix from file
  #exp_df = read.csv2("genefam_tester.csv", sep=",")[,-1]
  query_format = "TPM"
  
  # Prepare plot labels
  if (query_format == "FPKM") {
    str_exp_level = 'FPKM'
  }else if (query_format == "TPM"){
    str_exp_level = 'TPM'
  }
  
  exp_df$exp_level = as.numeric(exp_df$exp_level)
  value_level_max = max(exp_df$exp_level)
  value_level_min = min(exp_df$exp_level)
  
  exp_df_cast = acast(exp_df, sample_id~gene_id, value.var='exp_level')
  
  exp_df_dist = dist(exp_df_cast, method="euclidean")
  exp_df_hclust=hclust(exp_df_dist, method="complete")
  
  # Generate heatmap
  exp_heatmap = ggplot(data=exp_df, aes(x=gene_id, y=sample_id)) +
    geom_tile(alpha=0.6, aes(fill=exp_level), color='black', size=0.1, linetype='dashed') +
    
    theme_classic() +
    scale_fill_gradient2(low="#FFFFFF", mid="#f0c27b", high="#990033", midpoint = value_level_max/2) +
    #scale_fill_gradient(low="white", high="#3a1c71") +
    theme_minimal() +
    theme(#axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.x =element_text(angle=40,hjust=1,vjust=1.1),
    ) +
    geom_text(aes(label=round(exp_level, 1), color=exp_level), size=2)+
    scale_color_gradient2(low="white", mid="#eef2f3", high="black", midpoint = value_level_max/2) +
    guides(color='none',
           fill=guide_legend(title=str_exp_level)
           ) # hide legend for gradient text +
  
  
  tree_dendrogram = ggtree(exp_df_hclust, layout="rectangular",branch.length="none") +
    #geom_tiplab(hjust=-0.1) +
    geom_nodepoint(color="#a6e3e9", alpha=0.8, size=3) +
    theme_tree2() +
    xlim(NA, 6) # x axis limit
  
  heatmap_all = exp_heatmap %>%
    insert_left(tree_dendrogram, width = 0.1)
  
  ggsave("test2.png", plot=heatmap_all, dpi = 300, width = 9, height = 5)
  exp_df_cast_new = cbind(rownames(as.matrix(exp_df_cast)), as.matrix(exp_df_cast))
  write.csv(as.matrix(exp_df_cast_new), "gene_member_exp.csv", row.names = F)
  
}