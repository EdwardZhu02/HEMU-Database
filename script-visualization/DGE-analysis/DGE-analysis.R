rm(list=ls())
DE_data_raw = read.csv("./testdata.csv")

DE_group_list = c("B73","B73","B73","Fas","Fas","Fas")
DE_group_color_list = c("#2878B5", "#2878B5", "#2878B5", "#9AC9DB", "#9AC9DB", "#9AC9DB") # For normalization plot coloring

logfc_threshold = 2
pvalue_threshold = 0.05
heatmap_gene_count = 50
group1_name = 'B73'
group2_name = 'Fas'
output_folder_name = 'folder1'

#DE_analysis <- function (DE_data_raw, DE_group_list, logfc_threshold, pvalue_threshold, heatmap_gene_count,
                         #group1_name, group2_name, output_folder_name){

# Load packages
suppressMessages(library(limma)) # DE analysis
suppressMessages(library(tidyverse)) # Plot generate
suppressMessages(library(plotly)) # Interactive plot convert
suppressMessages(library(reshape2)) # Dataframe manipulation, melt raw expression dataframe
suppressMessages(library(FactoMineR)) # PCA
suppressMessages(library(factoextra)) # PCA
suppressMessages(library(pheatmap)) # Heatmap
suppressMessages(library(ggcorrplot)) # Correlation plot
suppressMessages(library(viridis)) # Theme customization

# Variable format conversion
logfc_threshold = as.numeric(logfc_threshold)
pvalue_threshold = as.numeric(pvalue_threshold)
heatmap_top_gene_count = as.numeric(heatmap_gene_count)
DE_group_list = as.character(DE_group_list)
group_list = factor(DE_group_list)
group_color_list = as.character(DE_group_color_list)

# Dataframe header removal
rownames(DE_data_raw) = DE_data_raw[,1]
DE_data_raw = DE_data_raw[,-1]
exprSet <- DE_data_raw

# Normalization of TPM data
exprSet_boxplt_data = melt(data=exprSet)
colnames(exprSet_boxplt_data) = c("sample_id", "Raw_TPM")

plt_boxplot_beforenorm = ggplot(data=exprSet_boxplt_data) + 
  geom_boxplot(aes(x=sample_id, y=Raw_TPM), outlier.shape = NA, alpha=0.7,
               outlier.size=0.1,
               fill = group_color_list) +
  coord_cartesian(ylim = c(0, 60)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
        axis.text.y = element_text(size = 10, hjust = 1),
        axis.title.x = element_blank()
  )

# Normalize Between Arrays
exprSet=normalizeBetweenArrays(exprSet)

exprSet_boxplt_data = melt(data=exprSet)
colnames(exprSet_boxplt_data) = c("sample_id", "Normalized_TPM")

plt_boxplot_afternorm = ggplot(data=exprSet_boxplt_data) + 
  geom_boxplot(aes(x=sample_id, y=Normalized_TPM), outlier.shape = NA, alpha=0.7,
               outlier.size=0.1,
               fill = group_color_list) +
  coord_cartesian(ylim = c(0, 60)) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
        axis.text.y = element_text(size = 10, hjust = 1),
        axis.title.x = element_blank()
  )

print("Normalization of TPM data complete.")

exprSet <- log2(exprSet+1)

# Execution of DE analysis
dat <- exprSet
design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
deg=topTable(fit,coef=2,adjust='BH',number = Inf)

if(T){
  logFC_t = logfc_threshold
  pvalue_t = pvalue_threshold
  deg$status=ifelse(deg$P.Value > pvalue_t,'NS',
                    ifelse( deg$logFC > logFC_t,'UP',
                            ifelse( deg$logFC < -logFC_t,'DOWN','NS') )
  )
  deg$symbol=rownames(deg)
}

### Output DEG results into a single .csv file.
write.csv(deg, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/differential_gene.csv"), row.names = FALSE)

### Density plot of expression ###
DEraw_density_data = melt(data=DE_data_raw)
colnames(DEraw_density_data) = c("sample_id", "Raw_TPM")
plt_sample_exp_density = ggplot(data=DEraw_density_data) +
  geom_density(aes(x=Raw_TPM, fill=sample_id), adjust=1.5) +
  coord_cartesian(xlim = c(0, 300)) +
  scale_fill_viridis(option = 'cividis', discrete = T) +
  theme_minimal() +
  facet_wrap(~sample_id) +
  theme(#panel.grid = element_blank(),
    #panel.background = element_blank(),
    axis.line = element_line(colour = 'black'),
    axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
    axis.text.y = element_text(size = 10, hjust = 1),
    legend.position = 'none',
  ) +
  labs(y="Density", x="TPM (Before Normalization)")

### Sample Correlation Plot ###
# Use TPM values after both limma normalization and log2(N+1) transformation
exp_correlation = cor(exprSet)
exp_correlation_p <- cor_pmat(exprSet)

plt_sample_correlation = ggcorrplot(
  exp_correlation,
  outline.color = "black",
  #method = "circle",
  ggtheme = ggplot2::theme_bw,
  lab_size = 4, p.mat = exp_correlation_p, sig.level = 0.05, insig = "blank", pch = 4, pch.col = "black", pch.cex = 5, tl.cex = 12,
  colors = c("#839EDB", "white", "#FF8D8D"),lab = T, title=paste0("Correlation: ", group1_name, " vs ", group2_name)) +
  theme(
    axis.text.x = element_text(size = 10, hjust = 1),
    axis.text.y = element_text(size = 10, hjust = 1),
  )


plt_volcano = ggplot(data=deg, aes(x=logFC, y=-log10(P.Value), color=status)) +
  geom_point(size=0.8, alpha=0.6) +
  scale_color_manual(values =c("DOWN"='#00AFBB',"NS"='#eaeaea',"UP"="#f38181")) +
  
  geom_vline(aes(xintercept=logFC_t),color="#40514e",size=0.5,linetype="dashed") +
  geom_vline(aes(xintercept=-logFC_t),color="#40514e",size=0.5,linetype="dashed") +
  geom_hline(aes(yintercept=-log10(pvalue_threshold)),color="#40514e",size=0.5,linetype="dashed") + # -log10(0.05)
  
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text.x = element_text(size = 10, hjust = 1),
        axis.text.y = element_text(size = 10, hjust = 1),
        #axis.title.x = element_blank(),
        legend.position = 'none',
  ) +
  labs(title=paste0("Volcano: ", group1_name, " vs ", group2_name))

### PCA Analysis ###

raw_PCA_data = t(exprSet)
dat.pca = PCA(raw_PCA_data, graph = F)

plt_pca_samples = fviz_pca_ind(dat.pca,
                               #geom.ind = "point", # show points only (but not "text")
                               col.ind = DE_group_list, # color by groups
                               labelsize = 3,
                               palette = c("#00AFBB", "#F8AC8C"),
                               addEllipses = TRUE, # Concentration ellipses 
                               legend.title = "Groups") +
  
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'),
        axis.text.x = element_text(size = 10, hjust = 1),
        axis.text.y = element_text(size = 10, hjust = 1),
        #legend.position = 'none',
  ) +
  labs(title=paste0("PCA-samples: ", group1_name, " vs ", group2_name))

### Heatmap of top most DE genes ###
print("Plotting heatmap")

DE_data_heatmap = deg[which(deg$status %in% c("UP", "DOWN")),]
DE_data_heatmap$logFC_abs = abs(DE_data_heatmap$logFC)
DE_data_heatmap = DE_data_heatmap[order(-DE_data_heatmap$logFC_abs),]

DE_data_heatmap_toplist = DE_data_heatmap[1:heatmap_top_gene_count,8]

TPM_data_raw_heatmap = DE_data_raw
TPM_data_raw_heatmap$symbol = rownames(TPM_data_raw_heatmap)
top_DE_rownames = rownames(TPM_data_raw_heatmap)


TPM_data_top_heatmap = TPM_data_raw_heatmap[which(TPM_data_raw_heatmap$symbol %in% DE_data_heatmap_toplist),-ncol(TPM_data_raw_heatmap)]
TPM_data_top_heatmap = as.matrix(TPM_data_top_heatmap)
#TPM_data_top_heatmap = apply(TPM_data_top_heatmap,2,as.numeric)

TPM_data_top_heatmap[TPM_data_top_heatmap==0] = NA
TPM_data_top_heatmap[is.na(TPM_data_top_heatmap)] = min(TPM_data_top_heatmap,na.rm = T)*0.01
TPM_data_top_heatmap = TPM_data_top_heatmap[apply(TPM_data_top_heatmap, 1, function(x) sd(x)!=0),]
# Eliminate pheatmap error: NA/NaN/Inf in foreign function call

rownames(TPM_data_raw_heatmap) = top_DE_rownames

plt_heatmap_topDE = pheatmap(log10(TPM_data_top_heatmap),
                             display_numbers = TRUE,
                             color = colorRampPalette(c("#00AFBB", "white", "#c93756"))(256),
                             fontsize=7, border_color = "grey60",
                             main = paste0("Heatmap - top DE genes: ", group1_name, " vs ", group2_name, "\n(log10 transformed)"),
                             treeheight_row=50, treeheight_col = 30,
                             cellwidth = 30, cellheight = 6,
                             silent=T,
)

#ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/boxplot_beforenorm.png"),
#       plot=plt_boxplot_beforenorm, dpi = 400, width=5, height=5)
#output_box_beforenorm <- plotly::ggplotly(plt_boxplot_beforenorm)
#htmlwidgets::saveWidget(output_box_beforenorm, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/boxplot_beforenorm.html"),
#                        selfcontained = TRUE)

#ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/boxplot_afternorm.png"),
#       plot=plt_boxplot_afternorm, dpi = 400, width=5, height=5)
#output_box_afternorm <- plotly::ggplotly(plt_boxplot_afternorm)
#htmlwidgets::saveWidget(output_box_afternorm, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/boxplot_afternorm.html"),
#                        selfcontained = TRUE)


ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/exp_density.png"),
       plot=plt_sample_exp_density, dpi = 400, width=5, height=5)
output_sample_exp_density <- plotly::ggplotly(plt_sample_exp_density)
htmlwidgets::saveWidget(output_sample_exp_density, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/exp_density.html"),
                        selfcontained = TRUE)

ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/exp_correlation.png"),
       plot=plt_sample_correlation, dpi = 400, width=5, height=5)
output_sample_correlation <- plotly::ggplotly(plt_sample_correlation)
htmlwidgets::saveWidget(output_sample_correlation, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/exp_correlation.html"),
                        selfcontained = TRUE)

output_volcano <- plotly::ggplotly(plt_volcano)
htmlwidgets::saveWidget(output_volcano, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/volcano.html"))
ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/volcano.png"), 
       plot=plt_volcano, dpi = 400, width=6, height=5)

output_pca_samples <- plotly::ggplotly(plt_pca_samples)
htmlwidgets::saveWidget(output_pca_samples, paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/pca_samples.html"))
ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/pca_samples.png"), 
       plot=plt_pca_samples, dpi = 400, width=6, height=5)

ggsave(paste0("Mainapp/static/Temp_R_DEprojects/", output_folder_name, "/top_de_genes_heatmap.png"), 
       plot=plt_heatmap_topDE, dpi = 400,
       width = length(DE_group_list), height = heatmap_top_gene_count/7)

