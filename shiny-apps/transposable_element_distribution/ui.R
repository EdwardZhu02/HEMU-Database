#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# App: TE distribution

# Load essential packages
suppressMessages(suppressWarnings(if (!require('bslib')) install.packages('bslib')))
suppressMessages(suppressWarnings(if (!require('shiny')) install.packages('shiny')))

# Define a function that generate UI layout for each panel
generate_panel <- function(panel_name, panel_number, choices, selected) {
  chrid1 <- as.numeric(panel_number)
  tabPanel(em(panel_name),
           sidebarLayout(
             sidebarPanel(
               # Sidebar configuration
               width = 4,
               # Sidebar content
               HTML("<h5>Hit the button below <strong style='color: darkred'>(Initialize Visualization Panel)</strong> to load TE annotation data and draw example plot.</h5>"),
               actionButton(paste0("action_generate_plt", chrid1), label = "Initialize Visualization Panel", class="btn-primary"),
               br(),
               br(),
               checkboxGroupInput(paste0("select_chrid", chrid1), label = "Sequences", 
                                  choices = choices,
                                  selected = selected),
               br(),
               sliderInput(paste0("slider_windowsize", chrid1), label = p("Window size (×10⁵bp)"), min = 0.5, 
                           max = 10, value = 3),
               br(),
               numericInput(paste0("plot_height", chrid1), label = "Plot height (px)", value = 400)
             ),
             mainPanel(
               fluidRow(
                 br(),
                 column(6,
                        numericInput(paste0("data1_height", chrid1), label = "Upper track height", value = 100),
                        numericInput(paste0("data2_height", chrid1), label = "Lower track height", value = 100),
                 ),
                 column(6,
                        numericInput(paste0("data1_inmargin", chrid1), label = "Upper track margin", value = 30),
                        numericInput(paste0("data2_inmargin", chrid1), label = "Lower track margin", value = 40),
                 ),
                 hr(),
                 br(),
                 p("LTR-RTs are classified into two superfamilies: RLG (Gypsy) and RLC (Copia)"),
                 p("TIR TEs are classified into five major superfamilies: DTA (hAT), DTC (CACTA), DTH (PIF/Harbinger), DTM (Mutator) and DTT(Tc1/Mariner)"),
               ),
               fluidRow(
                 hr(),
                 plotOutput(paste0("chr_dist_plt", chrid1), width="100%"),
               )
             )
           ) # sidebarLayout
  ) # tabPanel
}

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    theme = bs_theme(version = 4, bootswatch = "flatly"),
    tags$head(
      # Note the wrapping of the string in HTML()
      # Disable internal scroll bar when being embedded into an iframe
      tags$style(HTML("
      ::-webkit-scrollbar {
          display: none;
      }"))
    ),
    titlePanel(h4("TE Toolkit: Chromosomal TE Insertion Density Visualizer")),
    tabsetPanel(
      # Zea mays B73v4 ----
      generate_panel("Zea mays B73v4", 1, 
                     c( "chromosome 1" = 2, 
                        "chromosome 2" = 3, 
                        "chromosome 3" = 4,
                        "chromosome 4" = 5,
                        "chromosome 5" = 6,
                        "chromosome 6" = 7,
                        "chromosome 7" = 8,
                        "chromosome 8" = 9,
                        "chromosome 9" = 10,
                        "chromosome 10" = 1), 
                     c(2,3)
      ), # generate_panel
      # Sorghum bicolor BTx623 ----
      generate_panel("Sorghum bicolor BTx623", 2, 
                     c( "chromosome 1" = 2, 
                        "chromosome 2" = 3, 
                        "chromosome 3" = 4,
                        "chromosome 4" = 5,
                        "chromosome 5" = 6,
                        "chromosome 6" = 7,
                        "chromosome 7" = 8,
                        "chromosome 8" = 9,
                        "chromosome 9" = 10,
                        "chromosome 10" = 1), 
                     c(2,3)
                    ), # generate_panel
      # Coix lacryma-jobi var lacryma-jobi ----
      generate_panel("Coix lacryma-jobi var lacryma-jobi", 3, 
                     c( "chromosome 1" = 1, 
                        "chromosome 2" = 2, 
                        "chromosome 3" = 3,
                        "chromosome 4" = 4,
                        "chromosome 5" = 5,
                        "chromosome 6" = 6,
                        "chromosome 7" = 7,
                        "chromosome 8" = 8,
                        "chromosome 9" = 9,
                        "chromosome 10" = 10),
                     c(2,3)
                    ), # generate_panel
      # Miscanthus lutarioriparius ----
      generate_panel("Miscanthus lutarioriparius", 4, 
                     c( "chromosome 1" = 1, 
                        "chromosome 2" = 2, 
                        "chromosome 3" = 3,
                        "chromosome 4" = 4,
                        "chromosome 5" = 5,
                        "chromosome 6" = 6,
                        "chromosome 7" = 7,
                        "chromosome 8" = 8,
                        "chromosome 9" = 9,
                        "chromosome 10" = 10,
                        "chromosome 11" = 11, 
                        "chromosome 12" = 12, 
                        "chromosome 13" = 13,
                        "chromosome 14" = 14,
                        "chromosome 15" = 15,
                        "chromosome 16" = 16,
                        "chromosome 17" = 17,
                        "chromosome 18" = 18,
                        "chromosome 19" = 19
                     ),
                     c(2,3)
                    ), # generate_panel
      # Saccharum spontaneum ----
      generate_panel("Saccharum spontaneum", 5, 
                     c( "chromosome 1A" = 1, "chromosome 1B" = 2, "chromosome 1C" = 3, "chromosome 1D" = 4,
                        "chromosome 2A" = 5, "chromosome 2B" = 6, "chromosome 2C" = 7, "chromosome 2D" = 8,
                        "chromosome 3A" = 9, "chromosome 3B" = 10, "chromosome 3C" = 11, "chromosome 3D" = 12,
                        "chromosome 4A" = 13, "chromosome 4B" = 14, "chromosome 4C" = 15, "chromosome 4D" = 16,
                        "chromosome 5A" = 17, "chromosome 5B" = 18, "chromosome 5C" = 19, "chromosome 5D" = 20,
                        "chromosome 6A" = 21, "chromosome 6B" = 22, "chromosome 6C" = 23, "chromosome 6D" = 24,
                        "chromosome 7A" = 25, "chromosome 7B" = 26 , "chromosome 7C" = 27, "chromosome 7D" = 28,
                        "chromosome 8A" = 29, "chromosome 8B" = 30, "chromosome 8C" = 31, "chromosome 8D" = 32
                        ),
                     c(2,3)
      ), # generate_panel
      # Miscanthus sinensis ----
      generate_panel("Miscanthus sinesis", 6, 
                     c( "chromosome 1" = 1, 
                        "chromosome 2" = 2, 
                        "chromosome 3" = 3,
                        "chromosome 4" = 4,
                        "chromosome 5" = 5,
                        "chromosome 6" = 6,
                        "chromosome 7" = 7,
                        "chromosome 8" = 8,
                        "chromosome 9" = 9,
                        "chromosome 10" = 10,
                        "chromosome 11" = 11, 
                        "chromosome 12" = 12, 
                        "chromosome 13" = 13,
                        "chromosome 14" = 14,
                        "chromosome 15" = 15,
                        "chromosome 16" = 16,
                        "chromosome 17" = 17,
                        "chromosome 18" = 18,
                        "chromosome 19" = 19
                     ),
                     c(2,3)
      ) # generate_panel
      
    )
  )
)
