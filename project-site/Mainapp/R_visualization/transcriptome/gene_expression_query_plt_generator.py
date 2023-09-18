import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

# File deletion
import shutil
from pathlib import Path

# Used for passing pandas dataframe variable to R
pandas2ri.activate()  # Used for variable conversion (pd.DataFrame)


def overview_barplot(exp_data, gene_id, query_format, query_species):
    # Detect and delete existing plots
    tmp_plot_path = "Mainapp/static/Temp_R_html"
    for elm in Path(tmp_plot_path).glob(gene_id + '*'):
        elm.unlink() if elm.is_file() else shutil.rmtree(elm)  # delete folder and file

    robjects.globalenv['exp_data'] = exp_data
    plt_rscript_tpm = '''

        overview_barplot_T <- function (exp_data, gene_id){
          suppressMessages(library(ggplot2))
          suppressMessages(library(ggthemes))
          suppressMessages(library(plotly))
          suppressMessages(library(ggrepel))

          # Generate overview barplot using expression data from a single gene
          # sample_id fpkm tpm tissue_type

          exp_data$tissue_type <- as.factor(exp_data$tissue_type);
          exp_data$fpkm <- as.numeric(as.character(exp_data$fpkm));
          exp_data$tpm <- as.numeric(as.character(exp_data$tpm));
          exp_data <- exp_data %>% arrange(tissue_type, sample_id)
          # exp_data = exp_data %>% arrange(tissue_type, tpm)
          exp_data$id <- seq_len(nrow(exp_data))


          overview_barplot1 <- ggplot(data=exp_data) +
              geom_segment(aes(x=id, xend=id, y=0, yend=tpm),
                            color="#3366CC", alpha=0.4, size=1) +
              geom_point(aes(x=id, y=tpm, shape=sample_id), color="#66CCCC", size=0.4) +

              labs( x="Samples", y="TPM",caption = gene_id) +
              theme_bw() +
              theme(
                panel.border = element_blank(),
                legend.position="none",
              )
          output <- plotly::ggplotly(overview_barplot1)
          HTML_name <- paste0("Mainapp/static/Temp_R_html/", gene_id, "_plt1_TPM.html")
          htmlwidgets::saveWidget(output, HTML_name, selfcontained = TRUE)
          return(TRUE)
        }
        '''
    plt_rscript_fpkm = '''

            overview_barplot_F <- function (exp_data, gene_id){
              suppressMessages(library(ggplot2))
              suppressMessages(library(ggthemes))
              suppressMessages(library(plotly))
              suppressMessages(library(ggrepel))

              # Generate overview barplot using expression data from a single gene
              # sample_id fpkm tpm tissue_type

              exp_data$tissue_type <- as.factor(exp_data$tissue_type);
              exp_data$fpkm <- as.numeric(as.character(exp_data$fpkm));
              exp_data$tpm <- as.numeric(as.character(exp_data$tpm));
              exp_data <- exp_data %>% arrange(tissue_type, sample_id)
              # exp_data = exp_data %>% arrange(tissue_type, tpm)
              exp_data$id <- seq_len(nrow(exp_data))


              overview_barplot1 <- ggplot(data=exp_data) +
                  geom_segment(aes(x=id, xend=id, y=0, yend=fpkm),
                                color="#3366CC", alpha=0.4, size=1) +
                  geom_point(aes(x=id, y=fpkm, shape=sample_id), color="#66CCCC", size=0.4) +

                  labs( x="Samples", y="FPKM",caption = gene_id) +
                  theme_bw() +
                  theme(
                    panel.border = element_blank(),
                    legend.position="none",
                  )
              output <- plotly::ggplotly(overview_barplot1)
              HTML_name <- paste0("Mainapp/static/Temp_R_html/", gene_id, "_plt1_FPKM.html")
              htmlwidgets::saveWidget(output, HTML_name, selfcontained = TRUE)
              return(TRUE)
            }
            '''

    try:
        if query_format == 'TPM':
            robjects.r(plt_rscript_tpm)
            robjects.r['overview_barplot_T'](exp_data, gene_id)
        elif query_format == 'FPKM':
            robjects.r(plt_rscript_fpkm)
            robjects.r['overview_barplot_F'](exp_data, gene_id)
        else:
            # Frontend query type error
            raise SystemError
    except SystemError:
        return None


def tissue_specific_barplot(exp_data, gene_id, query_format, query_species):

    # Existing plots already deleted in the above function, considering these two functions are used simultaneously.
    # If used individually, additional cleaning mechanisms should be added.

    robjects.globalenv['exp_data'] = exp_data

    plt2_rscript_tpm = '''

        tissue_specific_barplot_T <- function (exp_data, gene_id){

          suppressMessages(library(ggplot2))
          suppressMessages(library(ggthemes))
          suppressMessages(library(plotly))
          suppressMessages(library(ggrepel))

          # Generate tissue-specific barplot using expression data from a single gene
          # sample_id fpkm tpm tissue_type

          exp_data$tissue_type <- as.factor(exp_data$tissue_type);
          exp_data$fpkm <- as.numeric(as.character(exp_data$fpkm));
          exp_data$tpm <- as.numeric(as.character(exp_data$tpm));
          exp_data <- exp_data %>% arrange(tissue_type, sample_id)
          # exp_data = exp_data %>% arrange(tissue_type, tpm)
          exp_data$id <- seq_len(nrow(exp_data))


        tissue_specific_barplot_t_big <- ggplot(data = exp_data, mapping = aes(x=tissue_type, y=tpm)) + 
                
          geom_jitter(
            color="#66CCCC",
            alpha=0.5,
            size=0.2,
            ) +
          geom_boxplot(
            alpha=0.5,
            fill="gray",
            outlier.fill="gray",
            outlier.size=0.3
          ) +
          theme_bw() +
          labs(
            x = "Tissue types",
            y = "TPM",
            caption = gene_id,
          ) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position="none",
            panel.border = element_blank(),
          ) 

          output <- plotly::ggplotly(tissue_specific_barplot_t_big)
          HTML_name <- paste0("Mainapp/static/Temp_R_html/", gene_id, "_plt2_TPM.html")
          htmlwidgets::saveWidget(output, HTML_name, selfcontained = TRUE)
          return(TRUE)
        }
        '''
    plt2_rscript_fpkm = '''

        tissue_specific_barplot_F <- function (exp_data, gene_id){

          suppressMessages(library(ggplot2))
          suppressMessages(library(ggthemes))
          suppressMessages(library(plotly))
          suppressMessages(library(ggrepel))

          # Generate tissue-specific barplot using expression data from a single gene
          # sample_id fpkm tpm tissue_type

          exp_data$tissue_type <- as.factor(exp_data$tissue_type);
          exp_data$fpkm <- as.numeric(as.character(exp_data$fpkm));
          exp_data$tpm <- as.numeric(as.character(exp_data$tpm));
          exp_data <- exp_data %>% arrange(tissue_type, sample_id)
          # exp_data = exp_data %>% arrange(tissue_type, tpm)
          exp_data$id <- seq_len(nrow(exp_data))


          tissue_specific_barplot_f_big <- ggplot(data = exp_data, mapping = aes(x=tissue_type, y=fpkm)) + 
                  geom_jitter(
                    color="#66CCCC",
                    alpha=0.5,
                    size=0.5,
                    ) +
                  geom_boxplot(
                    alpha=0.5,
                    fill="gray",
                    outlier.fill="gray",
                    outlier.size=0.3
                  ) +
                  theme_bw() +
                  labs(
                    x = "Tissue types",
                    y = "FPKM",
                    caption = gene_id,
                  ) +
                  theme(
                    axis.text.x = element_text(angle = 45, hjust = 1),
                    legend.position="none",
                    panel.border = element_blank(),
                  )


          output <- plotly::ggplotly(tissue_specific_barplot_f_big)
          HTML_name <- paste0("Mainapp/static/Temp_R_html/", gene_id, "_plt2_FPKM.html")
          htmlwidgets::saveWidget(output, HTML_name, selfcontained = TRUE)
          return(TRUE)
        }
        '''

    try:
        if query_format == 'TPM':
            robjects.r(plt2_rscript_tpm)
            robjects.r['tissue_specific_barplot_T'](exp_data, gene_id)
        elif query_format == 'FPKM':
            robjects.r(plt2_rscript_fpkm)
            robjects.r['tissue_specific_barplot_F'](exp_data, gene_id)
        else:
            # Frontend query type error
            return SystemError
    except SystemError:
        return None
