TE_sample_plotter_tpm <- function (te_dataframe, filename_id){

    rm(list=ls())

    suppressMessages(library(tidyverse))
    suppressMessages(library(plotly))

    te_total_expressed = te_dataframe
    output_id = filename_id

    te_total_expressed$TE_class = as.factor(te_total_expressed$TE_class)
    te_total_expressed$TE_class_group = as.factor(te_total_expressed$TE_class_group)
    te_total_expressed$fpkm = as.numeric(te_total_expressed$fpkm)
    te_total_expressed$tpm = as.numeric(te_total_expressed$tpm)

    write.csv(te_total_expressed, file = paste0("Mainapp/static/Temp_R_TE/", output_id, "_TEexpr_tpm.csv"))

    plt_TE_supfam = ggplot(data=te_total_expressed, aes(x=TE_class, y=tpm)) +
      geom_jitter(aes(color=TE_class_group), width=0.2, size=0.5, alpha=0.4) +
      geom_boxplot(aes(fill=TE_class_group), alpha=0.6, notch = T, outlier.shape = NA) +

      theme(axis.text.x = element_text(angle = 0, size = 8),) +
      theme_bw() +
      scale_fill_brewer(palette = "BrBG") +
      scale_color_brewer(palette = "BrBG") +
      coord_flip() +
      theme(
            #panel.grid = element_blank(),
            #panel.background = element_blank(),
            axis.line = element_line(colour = 'black'),
            axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
            axis.text.y = element_text(size = 10, hjust = 1),
        ) +
      ylim(1,200) +
      #facet_grid(TE_class_group~., drop = T) +
      xlab("TE Superfamily") + ylab("TPM") +
      labs(color = "Group", fill = "Group")


    plt_TE_group = ggplot(data=te_total_expressed, aes(x=TE_class_group, y=tpm)) +
      geom_jitter(aes(color=TE_class_group), width=0.1, size=0.5, alpha=0.2) +
      geom_boxplot(aes(fill=TE_class_group), alpha=0.6, notch = T, outlier.shape = NA) +
      scale_fill_brewer(palette = "BrBG") +
      scale_color_brewer(palette = "BrBG") +
      #coord_flip() +

      theme(
            panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black'),
            axis.text.x = element_text(size = 10, hjust = 1, angle = 60),
            axis.text.y = element_text(size = 10, hjust = 1),
            legend.position = 'none'
        ) +
      ylim(1,200) +
      #facet_grid(TE_class_group~., drop = T) +
      xlab("TE Group") + ylab("TPM")


    output_plt_TE_supfam <- plotly::ggplotly(plt_TE_supfam)
    htmlwidgets::saveWidget(output_plt_TE_supfam,
                            paste0("Mainapp/static/Temp_R_TE/", output_id,
                                   "_plt_TE_superfamily.html"),
                            selfcontained = TRUE)
    ggsave(paste0("Mainapp/static/Temp_R_TE/", output_id, "_plt_TE_superfamily.png"), plot=plt_TE_supfam)

    output_plt_TE_group <- plotly::ggplotly(plt_TE_group)
    htmlwidgets::saveWidget(output_plt_TE_group,
                            paste0("Mainapp/static/Temp_R_TE/", output_id,
                                   "_plt_TE_group.html"),
                            selfcontained = TRUE)
    ggsave(paste0("Mainapp/static/Temp_R_TE/", output_id, "_plt_TE_group.png"), plot=plt_TE_group)
    rm(list=ls())
}

TE_family_plotter_tpm <- function (te_dataframe, filename_id){

    suppressMessages(library(tidyverse))
    suppressMessages(library(patchwork))
    suppressMessages(library(plotly))

    te_raw_data = te_dataframe
    output_id = filename_id

    # colnames(te_raw_data) = c('sample_id', 'fpkm', 'tpm', 'sample_tissue')

    te_raw_data$fpkm = as.numeric(te_raw_data$fpkm)
    te_raw_data$tpm = as.numeric(te_raw_data$tpm)
    te_raw_data$sample_tissue = as.factor(te_raw_data$sample_tissue)
    te_raw_data$id = seq_len(nrow(te_raw_data))

    write.csv(te_raw_data, file = paste0("Mainapp/static/Temp_R_TE/", output_id, "_TEexpr.csv"))

    min_tpm = min(te_raw_data$tpm)
    mean_tpm = mean(te_raw_data$tpm)
    max_tpm = max(te_raw_data$tpm)
    sample_number = nrow(te_raw_data)
    sample_number_expressed = nrow(te_raw_data[which(te_raw_data$tpm > 1),])

    plt_TE_sample = ggplot(data=te_raw_data) +
      geom_segment(aes(x=id, xend=id, y=0, yend=tpm), color="#3366CC", alpha=0.4, size=1) +
      geom_point(aes(x=id, y=tpm, shape=sample_id), color="#66CCCC", size=0.4) +
      geom_hline(yintercept = mean_tpm, color="#f38181", linetype = "dashed", size=0.2) +
      geom_hline(yintercept = max_tpm, color="#f38181", linetype = "dashed", size=0.2) +
      #geom_text(aes(x=sample_number*2/3, y=mean_tpm*1.5),
      #          label=paste0("Mean: ", round(mean_tpm,2))) +
      #geom_text(aes(x=sample_number*2/3, y=max_tpm*0.9),
      #          label=paste0("Max: ", round(max_tpm,2))) +
      theme_bw() +
      theme(
        panel.border = element_blank(),
        legend.position="none",
      ) +
      labs( x=paste0("samples (n=", sample_number, ")"), y="TPM")

    plt_TE_tissue = ggplot(data = te_raw_data, mapping = aes(x=sample_tissue, y=tpm)) +
      geom_boxplot(alpha=0.5, fill="gray", outlier.fill="gray", outlier.size=0.3) +
      geom_jitter(color="#66CCCC", alpha=0.2, size=0.5, width=0.1) +
      theme_bw() +
      scale_fill_brewer() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none",
        panel.border = element_blank(),
      ) +
      labs(x = "Tissue types", y = "TPM")

    output_plt_TE_supfam <- plotly::ggplotly(plt_TE_sample)
    htmlwidgets::saveWidget(output_plt_TE_supfam,
                            paste0("Mainapp/static/Temp_R_TE/", output_id,
                                   "_plt_TE_famsample.html"),
                            selfcontained = TRUE)

    output_plt_TE_group <- plotly::ggplotly(plt_TE_tissue)
    htmlwidgets::saveWidget(output_plt_TE_group,
                            paste0("Mainapp/static/Temp_R_TE/", output_id,
                                   "_plt_TE_famtissue.html"),
                            selfcontained = TRUE)

    return(list(min_tpm, mean_tpm, max_tpm, sample_number, sample_number_expressed))
    rm(list=ls())
}

TE_sample_plotter_fpkm <- function (te_dataframe, filename_id){

    suppressMessages(library(tidyverse))
    suppressMessages(library(plotly))

    te_total_expressed = te_dataframe
    output_id = filename_id

    te_total_expressed$TE_class = as.factor(te_total_expressed$TE_class)
    te_total_expressed$TE_class_group = as.factor(te_total_expressed$TE_class_group)
    te_total_expressed$fpkm = as.numeric(te_total_expressed$fpkm)
    te_total_expressed$tpm = as.numeric(te_total_expressed$tpm)

    write.csv(te_total_expressed, file = paste0("Mainapp/static/Temp_R_TE/", output_id, "_TEexpr.csv"))

    plt_TE_supfam = ggplot(data=te_total_expressed, aes(x=TE_class, y=fpkm)) +
      geom_jitter(aes(color=TE_class_group), width=0.2, size=0.5, alpha=0.4) +
      geom_boxplot(aes(fill=TE_class_group), alpha=0.6, notch = T, outlier.shape = NA) +

      theme(axis.text.x = element_text(angle = 0, size = 8),) +
      theme_bw() +
      scale_fill_brewer(palette = "BrBG") +
      scale_color_brewer(palette = "BrBG") +
      coord_flip() +
      theme(
            #panel.grid = element_blank(),
            #panel.background = element_blank(),
            axis.line = element_line(colour = 'black'),
            axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
            axis.text.y = element_text(size = 10, hjust = 1),
        ) +
      ylim(1,200) +
      #facet_grid(TE_class_group~., drop = T) +
      xlab("TE Superfamily") + ylab("FPKM") +
      labs(color = "Group", fill = "Group")


    plt_TE_group = ggplot(data=te_total_expressed, aes(x=TE_class_group, y=fpkm)) +
      geom_jitter(aes(color=TE_class_group), width=0.1, size=0.5, alpha=0.2) +
      geom_boxplot(aes(fill=TE_class_group), alpha=0.6, notch = T, outlier.shape = NA) +
      scale_fill_brewer(palette = "BrBG") +
      scale_color_brewer(palette = "BrBG") +
      #coord_flip() +

      theme(
            panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black'),
            axis.text.x = element_text(size = 10, hjust = 1, angle = 60),
            axis.text.y = element_text(size = 10, hjust = 1),
            legend.position = 'none'
        ) +
      ylim(1,200) +
      #facet_grid(TE_class_group~., drop = T) +
      xlab("TE Group") + ylab("FPKM")


    output_plt_TE_supfam <- plotly::ggplotly(plt_TE_supfam)
    htmlwidgets::saveWidget(output_plt_TE_supfam,
                            paste0("Mainapp/static/Temp_R_TE/", output_id,
                                   "_plt_TE_superfamily.html"),
                            selfcontained = TRUE)
    ggsave(paste0("Mainapp/static/Temp_R_TE/", output_id, "_plt_TE_superfamily.png"), plot=plt_TE_supfam)

    output_plt_TE_group <- plotly::ggplotly(plt_TE_group)
    htmlwidgets::saveWidget(output_plt_TE_group,
                            paste0("Mainapp/static/Temp_R_TE/", output_id,
                                   "_plt_TE_group.html"),
                            selfcontained = TRUE)
    ggsave(paste0("Mainapp/static/Temp_R_TE/", output_id, "_plt_TE_group.png"), plot=plt_TE_group)
    rm(list=ls())
}

TE_family_plotter_fpkm <- function (te_dataframe, filename_id){

    suppressMessages(library(tidyverse))
    suppressMessages(library(patchwork))
    suppressMessages(library(plotly))

    te_raw_data = te_dataframe
    output_id = filename_id

    # colnames(te_raw_data) = c('sample_id', 'fpkm', 'tpm', 'sample_tissue')

    te_raw_data$fpkm = as.numeric(te_raw_data$fpkm)
    te_raw_data$tpm = as.numeric(te_raw_data$tpm)
    te_raw_data$sample_tissue = as.factor(te_raw_data$sample_tissue)
    te_raw_data$id = seq_len(nrow(te_raw_data))

    write.csv(te_raw_data, file = paste0("Mainapp/static/Temp_R_TE/", output_id, "_TEexpr.csv"))

    min_fpkm = min(te_raw_data$fpkm)
    mean_fpkm = mean(te_raw_data$fpkm)
    max_fpkm = max(te_raw_data$fpkm)
    sample_number = nrow(te_raw_data)
    sample_number_expressed = nrow(te_raw_data[which(te_raw_data$fpkm > 1),])

    plt_TE_sample = ggplot(data=te_raw_data) +
      geom_segment(aes(x=id, xend=id, y=0, yend=fpkm), color="#3366CC", alpha=0.4, size=1) +
      geom_point(aes(x=id, y=fpkm, shape=sample_id), color="#66CCCC", size=0.4) +
      geom_hline(yintercept = mean_fpkm, color="#f38181", linetype = "dashed", size=0.2) +
      geom_hline(yintercept = max_fpkm, color="#f38181", linetype = "dashed", size=0.2) +
      #geom_text(aes(x=sample_number*2/3, y=mean_fpkm*1.5),
      #          label=paste0("Mean: ", round(mean_fpkm,2))) +
      #geom_text(aes(x=sample_number*2/3, y=max_fpkm*0.9),
      #          label=paste0("Max: ", round(max_fpkm,2))) +
      theme_bw() +
      theme(
        panel.border = element_blank(),
        legend.position="none",
      ) +
      labs( x=paste0("samples (n=", sample_number, ")"), y="FPKM")

    plt_TE_tissue = ggplot(data = te_raw_data, mapping = aes(x=sample_tissue, y=fpkm)) +
      geom_boxplot(alpha=0.5, fill="gray", outlier.fill="gray", outlier.size=0.3) +
      geom_jitter(color="#66CCCC", alpha=0.2, size=0.5, width=0.1) +
      theme_bw() +
      scale_fill_brewer() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="none",
        panel.border = element_blank(),
      ) +
      labs(x = "Tissue types", y = "FPKM")

    output_plt_TE_supfam <- plotly::ggplotly(plt_TE_sample)
    htmlwidgets::saveWidget(output_plt_TE_supfam,
                            paste0("Mainapp/static/Temp_R_TE/", output_id,
                                   "_plt_TE_famsample.html"),
                            selfcontained = TRUE)

    output_plt_TE_group <- plotly::ggplotly(plt_TE_tissue)
    htmlwidgets::saveWidget(output_plt_TE_group,
                            paste0("Mainapp/static/Temp_R_TE/", output_id,
                                   "_plt_TE_famtissue.html"),
                            selfcontained = TRUE)

    return(list(min_fpkm, mean_fpkm, max_fpkm, sample_number, sample_number_expressed))
    rm(list=ls())
}