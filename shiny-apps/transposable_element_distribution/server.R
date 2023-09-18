#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# App: TE distribution

# Load essential packages
suppressMessages(suppressWarnings(if (!require('ggplot2')) install.packages('ggplot2')));
suppressMessages(suppressWarnings(if (!require('shiny')) install.packages('shiny')));
suppressMessages(suppressWarnings(if (!require('karyoploteR')) BiocManager::install('karyoploteR',update = FALSE)));
suppressMessages(suppressWarnings(if (!require('rtracklayer')) BiocManager::install('rtracklayer',update = FALSE)));

# Define server logic
shinyServer(function(input, output) {

  # Define function to generate plot for a specific genome
  generatePlot <- function(genome, pltheight, slider_windowsize, panel_number) {
    
    pltgen_progress_name = c(
      "Render plot"
    )
    withProgress(message = 'Plot update', value = 0, expr = {
      for (i in 1:2) {
      
      # Update progress bar and display task information
      incProgress(1/2, detail = pltgen_progress_name[i])
      
      if (i == 1) { # Load TE distribution data
        
        # Define color used for plotting
        color1 = '#A5B6C5' # LTR
        color2 = '#F0988C' # Helitron
        color3 = '#496C88' # DNA
        
        plt_id = as.numeric(panel_number) # Obtain plot ID for output
        load(paste0("./Resources/", genome, ".RData"))
        
      } else { # Render plot
        # Render plot ----
        output[[paste0("chr_dist_plt", plt_id)]] <<- renderPlot({
          
          # Curate sequences
          selected_sequences = input[[paste0("select_chrid", plt_id)]]
          seq_ranges <- toGRanges(sequence_fai_chr[selected_sequences,c(1,2,3)])
          
          # specify plot parameters
          plot_parameters = getDefaultPlotParams(plot.type=2)
          plot_parameters$data1height = input[[paste0("data1_height", plt_id)]]
          plot_parameters$data2height = input[[paste0("data2_height", plt_id)]]
          plot_parameters$data1inmargin = input[[paste0("data1_inmargin", plt_id)]]
          plot_parameters$data2inmargin = input[[paste0("data2_inmargin", plt_id)]]
          plot_parameters$data1outmargin = 100
          plot_parameters$data2outmargin = 100
          
          # Raw chr sequence rectangles ----
          kp = plotKaryotype(genome=seq_ranges, plot.type=2,
                             cex=0.7, # Font sizes for sequence names
                             plot.params=plot_parameters)
          kpAddBaseNumbers(kp)
          kpAddCytobandsAsLine(kp)
          
          # LTR-RT ----
          # Copia
          kpPlotDensity(kp, feature_copiaLTR, window.size = as.numeric(input[[paste0("slider_windowsize", plt_id)]])*10e5, 
                        col=color1, border=color1,
                        data.panel=2, # upside
                        r0=0, r1=0.2)
          kpAddLabels(kp, "RLC        ", data.panel=2, cex=0.7, col="#888888",
                      r0=0, r1=0.2)
          kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, data.panel=2, r0=0, r1=0.2, cex=0.4)
          
          # Gypsy
          kpPlotDensity(kp, feature_gypsyLTR, window.size = as.numeric(input[[paste0("slider_windowsize", plt_id)]])*10e5,
                        col=color1, border=color1,
                        data.panel=2, # upside
                        r0=0.4, r1=0.6)
          kpAddLabels(kp, "RLG        ", data.panel=2, cex=0.7, col="#888888",
                      r0=0.4, r1=0.6)
          kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, data.panel=2, r0=0.4, r1=0.6, cex=0.4)
          
          
          # Helitron ----
          kpPlotDensity(kp, feature_helitron, window.size = as.numeric(input[[paste0("slider_windowsize", plt_id)]])*10e5, 
                        col=color2, border=color2,
                        data.panel=2, # upside
                        r0=0.8, r1=1.0)
          kpAddLabels(kp, "helitron      ", data.panel=2, cex=0.7, col="#888888",
                      r0=0.8, r1=1.0)
          kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, data.panel=2, r0=0.8, r1=1.0, cex=0.4)
          
          
          # DNA-TIR ----
          # hAT/DTA
          kpPlotDensity(kp, feature_DTATIR, window.size = as.numeric(input[[paste0("slider_windowsize", plt_id)]])*10e5, 
                        col=color3, border=color3,
                        data.panel=1, # upside
                        r0=0, r1=0.2)
          kpAddLabels(kp, "DTA        ", data.panel=1, cex=0.7, col="#888888",
                      r0=0, r1=0.2)
          kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, data.panel=1, r0=0, r1=0.2, cex=0.4)
          
          # PIF/DTH
          kpPlotDensity(kp, feature_DTHTIR, window.size = as.numeric(input[[paste0("slider_windowsize", plt_id)]])*10e5, 
                        col=color3, border=color3,
                        data.panel=1, # upside
                        r0=0.4, r1=0.6)
          kpAddLabels(kp, "DTH        ", data.panel=1, cex=0.7, col="#888888",
                      r0=0.4, r1=0.6)
          kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, data.panel=1, r0=0.4, r1=0.6, cex=0.4)
          
          # Mutator/DTM
          kpPlotDensity(kp, feature_DTMTIR, window.size = as.numeric(input[[paste0("slider_windowsize", plt_id)]])*10e5, 
                        col=color3, border=color3,
                        data.panel=1, # upside
                        r0=0.8, r1=1)
          kpAddLabels(kp, "DTM        ", data.panel=1, cex=0.7, col="#888888",
                      r0=0.8, r1=1)
          kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, data.panel=1, r0=0.8, r1=1, cex=0.4)
          
          # CACTA/DTC
          kpPlotDensity(kp, feature_DTCTIR, window.size = as.numeric(input[[paste0("slider_windowsize", plt_id)]])*10e5, 
                        col=color3, border=color3,
                        data.panel=1, # upside
                        r0=1.2, r1=1.4)
          kpAddLabels(kp, "DTC        ", data.panel=1, cex=0.7, col="#888888",
                      r0=1.2, r1=1.4)
          kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, data.panel=1, r0=1.2, r1=1.4, cex=0.4)
          
          # Tc1/DTT
          kpPlotDensity(kp, feature_DTTTIR, window.size = as.numeric(input[[paste0("slider_windowsize", plt_id)]])*10e5, 
                        col=color3, border=color3,
                        data.panel=1, # upside
                        r0=1.6, r1=1.8)
          kpAddLabels(kp, "DTT        ", data.panel=1, cex=0.7, col="#888888",
                      r0=1.6, r1=1.8)
          kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, data.panel=1, r0=1.6, r1=1.8, cex=0.4)
          
        }, height = pltheight ) # renderPlot
      } # else
    } # for
    })
  } # function generatePlot
  
  # Zea mays B73v4 ----
  observeEvent(
    input$action_generate_plt1,
    {
      generatePlot("GCF_000005005.2", reactive(input$plot_height1), input$slider_windowsize1, 1)
    })
  
  # Sorghum bicolor BTx623 ----
  observeEvent(
    input$action_generate_plt2,
    {
      generatePlot("GCF_000003195.3", reactive(input$plot_height2), input$slider_windowsize2, 2)
    })
  
  # Coix lacryma-jobi var lacryma-jobi ----
  observeEvent(
    input$action_generate_plt3,
    {
      generatePlot("Clacr", reactive(input$plot_height3), input$slider_windowsize3, 3)
    })
  
  # Miscanthus lutarioriparius ----
  observeEvent(
    input$action_generate_plt4,
    {
      generatePlot("GCA_904845875.1", reactive(input$plot_height4), input$slider_windowsize4, 4)
    })
  
  # Saccharum spontaneum ----
  observeEvent(
    input$action_generate_plt5,
    {
      generatePlot("GCA_003544955.1", reactive(input$plot_height5), input$slider_windowsize5, 5)
    })
  
  # Miscanthus sinensis ----
  observeEvent(
    input$action_generate_plt6,
    {
      generatePlot("Msine", reactive(input$plot_height6), input$slider_windowsize6, 6)
    })

})
