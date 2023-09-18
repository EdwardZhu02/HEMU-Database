
genefam_phylogenetic_analysis <- function(wkdir_name, msa_method, pairwisedist_method, phylotree_layout, bootstrap_rep) {

  rm(list=ls())
  # Enforcing multiple sequence alignment
  # Params:
  # wkdir_name: "phylo" + 8-digit random number, ex,. phylo12345678
  # msa_method: ClustalW, CluatalOmega, Muscle
  # pairwisedist_method: K80, K81, T92
  # phylotree_layout: rectangular, slanted, circular

  # Testing
  # msa_method = "ClustalW"
  # pairwisedist_method = "K80"
  # phylotree_layout = "rectangular"
  # bootstrap_rep = 1000
  
  # Load packages
  # Hierarchical clustering
  suppressMessages(library(ape))
  #suppressMessages(library(msa))
  suppressMessages(suppressWarnings(if (!require('msa')) BiocManager::install('msa',update = FALSE)));
  # Common data processing
  suppressMessages(library(Biostrings))
  suppressMessages(library(reticulate))
  # Data visualization
  suppressMessages(library(ggtree))
  suppressMessages(library(pheatmap))
  
  # Set work directory and read sequence (.fasta format) - Abandoned, may lead to system bug
  # setwd(paste0("Mainapp/static/Temp_R_genefam/", wkdir_name))
  seq_cluster_data = readAAStringSet(paste0("Mainapp/static/Temp_R_genefam/", wkdir_name, "/seq_original.fasta"))
  
  # MSA, normalize sequence lengths for dendrogram construction
  seq_cluster_msa = msa(seq_cluster_data, method = msa_method)
  seq_cluster_dnabin = as.AAbin(seq_cluster_msa)
  
  # Construct neighbour-joining tree
  seq_cluster_dist = as.matrix(dist.aa(seq_cluster_dnabin))
  # If the sequences are very different, most evolutionary distances are undefined
  # and a non-finite value (Inf or NaN) is returned.
  # You may do dist.dna(, model = "raw") to check whether some values are higher than 0.75.

  # Sequence distance heatmap construction
  png(filename = paste0("Mainapp/static/Temp_R_genefam/", wkdir_name, "/dist_heatmap.png"), width=3500, height=3500, res=330)
  print(pheatmap(seq_cluster_dist, color = colorRampPalette(c("#00AFBB", "white", "#c93756"))(256), border_color = "grey60"))
  dev.off()

  # Normalize genomic distances
  seq_cluster_dist[sapply(seq_cluster_dist,is.infinite)] = NaN
  seq_cluster_nj = bionjs(seq_cluster_dist)

  # ---Bootstrap validation of neighbour-joining tree---
  # Define bootstrap function, independent to NJ dendrogram construction.
  func_bootstrap = function(x) {
    tmpmatr1 = dist.dna(x, model=pairwisedist_method, gamma=T, as.matrix=T)
    # Normalize genomic distances that are significantly bigger or smaller than
    # the pre-set limit, which will influence data visualization.
    tmpmatr1[sapply(tmpmatr1,is.infinite)] = NaN
    #print(tmpmatr1)
    # Construct neighbour-joining tree
    bionjs(tmpmatr1)
  }

  # Tree Bipartition and Bootstrapping Phylogenies
  seq_cluster_nj_bootstrap = boot.phylo(seq_cluster_nj, seq_cluster_dnabin, 
                                        func_bootstrap, quiet=T, B=as.numeric(bootstrap_rep), trees = T)
  
  # ---Process previously-constructed neighbour-joining tree---
  #seq_cluster_nj_bootstrap_backup = seq_cluster_nj_bootstrap
  # Add bootstrap values to the constructed tree
  # Conducted 1000 bootstrap tests, converting to hit in 100 tests
  seq_cluster_nj$node.label = round(seq_cluster_nj_bootstrap$BP/as.numeric(bootstrap_rep)*100, 0)
  
  options(ignore.negative.edge=T)
  #seq_cluster_nj_withroot = root(seq_cluster_nj, "stem_sequence")
  
  plt_njtree = ggtree(seq_cluster_nj, size=1, layout=phylotree_layout) +
    geom_tiplab(hjust=-0.1, fontface="italic") +
    geom_nodepoint(color="#a6e3e9", alpha=0.5, size=10) +
    geom_nodelab(hjust=-0.2, size=4, color="#b83b5e") + # Bootstrap values
    #geom_text(aes(label=node)) + # Node number
    #geom_highlight(node=21, fill="red", alpha=0.3)+
    # geom_strip(6,11,label = "AAA",offset = 0.1,offset.text = 0.02,
    #            color = "green",barsize = 3,fontsize = 5,angle = 90,
    #            hjust = 0.5) +
    theme_tree2() # Theme to display x axis
    #xlim(NA, 10) + # x axis lmit
    #geom_treescale(x=9, y=2, color="#3f72af") # Scale bar
  
  #plt_njtree
  # Save phylogenetic tree to static image
  ggsave(paste0("Mainapp/static/Temp_R_genefam/", wkdir_name, "/phylotree_nj.png"), plot=plt_njtree,
         dpi = 300, width = 10, height = 7, units = "in", scale = 0.7)
  # Save editable tree format (.nwk)
  write.tree(seq_cluster_nj, file = paste0("Mainapp/static/Temp_R_genefam/", wkdir_name, "/phylotree_nj.nwk"))

  # CLear workspace
  rm(list=ls())
}


