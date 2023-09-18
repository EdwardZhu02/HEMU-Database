#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# App: Gene structure visualization

# Load essential packages
suppressMessages(suppressWarnings(if (!require('shiny')) install.packages('shiny')))
suppressMessages(suppressWarnings(if (!require('ggplot2')) install.packages('ggplot2')))
suppressMessages(suppressWarnings(if (!require('dplyr')) install.packages('dplyr')))
suppressMessages(suppressWarnings(if (!require('stringr')) install.packages('stringr'))) # str_split
suppressMessages(suppressWarnings(if (!require('ggh4x')) install.packages('ggh4x')))
suppressMessages(suppressWarnings(if (!require('ggbio')) BiocManager::install('ggbio',update = FALSE)))
suppressMessages(suppressWarnings(if (!require('GenomicRanges')) BiocManager::install('GenomicRanges',update = FALSE)))

# Define server logic
shinyServer(function(input, output) {
  
  df_filtered = reactive({
    
    # Read genome annotation
    load(file = paste0("./Resources/", as.character(input$species_select), ".RData"))
    # genome_gff = paste0("./Resources/", as.character(input$species_select), "_filtered.gff")
    # df = read.table(genome_gff, header=F, sep="\t", quote = "")
    
    return(df)
  })
  
  # Define function to generate plot
  generatePlot <- function(df_filtered, gene_id, selected_species) {
    
    # Clear error messages
    output$error_message = renderText(" ")

    pltgen_progress_name = c(
      "Filter gene information",
      "Extract transcripts",
      "Render plot"
    )
    withProgress(message = 'Plot update', value = 0, expr = {
      for (i in 1:3) {
        
        # Update progress bar and display task information
        incProgress(1/3, detail = pltgen_progress_name[i])
        
        if (i == 1) { # Filter gene information
          
          gene_name_final = gene_id
          
          # Filter target gene (Exact match)
          df_containing_targetgene = df_filtered[which(df_filtered$V10 == input$gene_name),]
          
          # Validate if the program can obtain gene information from the given GFF
          if (nrow(df_containing_targetgene) == 0){
            output$error_message = renderUI({
              tagList(
                HTML("<h4>Gene ID not found in the selected species!</h4>"),
                HTML("No entry detected in genomic annotation of the selected species, please <strong style='color: red'>CHECK GENE ID QUERY</strong> <strong>and submit again.</strong>")
              )
            })
            output$df_targetgene = renderDataTable(as.data.frame(matrix(nrow=0,ncol=7)))
            return()
          }
          
          # Identify unique transcript using the first unit of attribute column (Example: ID) in simplified GFF
          # The split rule might be different among annotations.
          # Result example: Zm00001d027233_T001
          # ===========
          if (as.character(selected_species) == "Clacr"){
            df_targetgene_attrs = lapply(df_containing_targetgene$V9, function(x){
              str_split(
                str_split(
                  str_split(x, 
                    ";", simplify = T)[[1]], # ID column
                  ":", simplify = T)[[1]],
                "=", simplify = T)[[2]]
            })
          } else if (as.character(selected_species) == "Sbico"){
            df_targetgene_attrs = lapply(df_containing_targetgene$V9, function(x){
              str_split(
                str_split(x,
                  ";", simplify = T)[[1]], # ID column (CDS), Parent column (3UTR, 5UTR)
                "=", simplify = T)[[2]]
            })
          } else if (as.character(selected_species) == "Sspon"){
            df_targetgene_attrs = lapply(df_containing_targetgene$V9, function(x){
              str_split(
                str_split(
                  str_split(x, 
                    ";", simplify = T)[[1]], # ID column
                  ":", simplify = T)[[1]],
                "=", simplify = T)[[2]]
            })
          } else if (as.character(selected_species) == "Mluta"){
            df_targetgene_attrs = lapply(df_containing_targetgene$V9, function(x){
              str_split(
                str_split(x,
                  ";", simplify = T)[[1]], # ID column (CDS), Parent column (3UTR, 5UTR)
                "=", simplify = T)[[2]]
            })
          } else if (as.character(selected_species) == "Msine"){
            df_targetgene_attrs = lapply(df_containing_targetgene$V9, function(x){
              str_split(
                str_split(x,
                  ";", simplify = T)[[1]], # ID column (CDS), Parent column (3UTR, 5UTR)
                "=", simplify = T)[[2]]
            })
          } else if (as.character(selected_species) %in% c("Cserr", "Hdipl", "Ttria")){
            df_targetgene_attrs = lapply(df_containing_targetgene$V9, function(x){
              str_split(
                str_split(
                  str_split(x, 
                    ";", simplify = T)[[1]], # ID column
                  ":", simplify = T)[[1]],
                "=", simplify = T)[[2]]
            })
          } else if (as.character(selected_species) %in% 
                     c("Zmaysv4", "Zmaysv5", "Zmaysv3", "ZmaysPE0075", 
                       "ZmaysDK105", "ZmaysOh7B", "ZmaysP39", "ZmaysCML247", 
                       "ZmaysCML277", "ZmaysM162W", "ZmaysNC350", "ZmaysOH43", 
                       "ZmaysCML333", "ZmaysNC358", "ZmaysKy21", "ZmaysCML103", 
                       "ZmaysIl14H", "ZmaysB97", "ZmaysHP301", "ZmaysCML322", 
                       "ZmaysKi3", "ZmaysKi11", "ZmaysCML69", "ZmaysCML228", 
                       "ZmaysMo18W", "ZmaysM37W", "ZmaysTzi8", "ZmaysTx303", 
                       "ZmaysMs71", "ZmaysLH244", "ZmaysW22", "ZmaysEP1", 
                       "ZmaysF7")){
            df_targetgene_attrs = lapply(df_containing_targetgene$V9, function(x){
              str_split(
                str_split(x,
                          ";", simplify = T)[[1]], # ID column (CDS), Parent column (3UTR, 5UTR)
                "=", simplify = T)[[2]]
            })
          } else {
            output$error_message = renderUI({
              tagList(
                HTML("<h4>Annotation filter rule error!</h4>")
              )
            })
            output$df_targetgene = renderDataTable(as.data.frame(matrix(nrow=0,ncol=7)))
            return()
          }
          
          df_targetgene_attrs = t(as.data.frame(df_targetgene_attrs))
          colnames(df_targetgene_attrs) = c("transcript_id")
          # Merge with previously filtered dataframes
          df_containing_targetgene = cbind(df_containing_targetgene, df_targetgene_attrs)
          # ===========
          
          # Extract the 'gene' label row
          dftmp_generow = df_containing_targetgene[which(df_containing_targetgene$V3 == "gene"),]
          gene_start = as.numeric(dftmp_generow$V4)
          gene_end = as.numeric(dftmp_generow$V5)
          
          # Removing the 'gene' label row from the main df
          df_containing_targetgene = df_containing_targetgene[-which(df_containing_targetgene$V3 %in% c("gene","mRNA","exon")),]
          
        } else if (i == 2){
          
          # Sort row by transcript ID and element start position
          df_containing_targetgene = df_containing_targetgene[order(df_containing_targetgene$transcript_id, df_containing_targetgene$V4),]
          
          df_output = df_containing_targetgene[,-c(6,8,10,11)]
          colnames(df_output) = c("seqid", "source", "type", "start", "end", "strand", "attributes")
          
          # Validating potential non-translating gene
          if (nrow(df_output) == 0){
            # Probably a non-translating gene
            output$error_message = renderUI({
              tagList(
                HTML("<h4>Non-translating gene!</h4>"),
                HTML("The query gene ID represents a <strong style='color: red'>non-translating gene</strong>. Thus, no corresponding transcript is produced or identified.")
              )
            })
            output$df_targetgene = renderDataTable(as.data.frame(matrix(nrow=0,ncol=7)))
            return()
          }
          
          # Render frontend table
          output$df_targetgene = renderDataTable(df_output, options = list(pageLength = 10))
          
          # Create a GRanges object from the dataframe
          gene_range = GRanges(seqnames = dftmp_generow$V1[1], 
                               ranges = IRanges(start=df_containing_targetgene$V4, end=df_containing_targetgene$V5, group=df_containing_targetgene$V3),
                               strand = dftmp_generow$V7[1],
                               mcols = DataFrame(group = df_containing_targetgene$V3, transcript_id = df_containing_targetgene$transcript_id)
          )
        } else{
          
          # Generate the main plot
          pltmain = autoplot(gene_range,aes(group=mcols.transcript_id, fill=mcols.group),geom="alignment") +
            theme_bw() +
            scale_fill_brewer(palette = "Paired") +
            labs(fill = paste0("Feature:")) + 
            theme(legend.position = "bottom")
          return(pltmain)
        }
      }
    })
    
  }
  
  observeEvent(
    input$action_generate_plt,
    {
      output$plot <- renderPlot({
        generatePlot(
          df_filtered(),
          input$gene_name,
          input$species_select
        )
      })
    })
  
}) # shinyServer
