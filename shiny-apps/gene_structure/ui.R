#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# App: Gene structure visualization

# Load essential packages
suppressMessages(suppressWarnings(if (!require('bslib')) install.packages('bslib')))
suppressMessages(suppressWarnings(if (!require('shiny')) install.packages('shiny')))


# Define UI
shinyUI(
  fluidPage(
    theme = bs_theme(version = 4, bootswatch = "flatly"),
    tags$head(
      # Note the wrapping of the string in HTML()
      # Disable internal scroll bar when being embedded into an iframe
      tags$style(HTML("
      ::-webkit-scrollbar {
          display: none;
      }
      .selectize-input {
        height: 35px;
        padding-top: 5px;
      }"))
    ),
    # titlePanel(h4("Gene family analysis toolkit: Gene structure visualizer")),
  
    fluidRow(
      column(12, align="center",
        br(),
        # Species to load genome annotation from
        selectInput("species_select",
                    label=("Select species gene dataset"), width = "500px",
                    list("Zea mays B73v4 (GCF_000005005.2)" = "Zmaysv4",
                         "Zea mays B73v5 (GCF_902167145.1)" = "Zmaysv5",
                         "Coix lacryma-jobi var.lacryma-jobi (GCA_009763385.1)" = "Clacr",
                         "Saccharum spontaneum (GCA_003544955.1)" = "Sspon",
                         "Sorghum bicolor BTx623 (GCF_000003195.3)" = "Sbico",
                         "Saccharum spontaneum (GCA_003544955.1)" = "Sspon",
                         "Miscanthus lutarioriparius (GCA_904845875.1)" = "Mluta",
                         "Miscanthus sinensis (Phytozome v7.0)" = "Msine",
                         "Chrysopogon serrulatus (GCA_015844335.1)" = "Cserr",
                         "Hyparrhenia diplandra (GCA_015847255.1)" = "Hdipl",
                         "Themeda triandra (GCA_018135685.1)" = "Ttria",
                         "Zea mays B73v3 (GCF_000005005.1)" = "Zmaysv3",
                         "Zea mays PE0075 (GCA_003704525.1)" = "ZmaysPE0075",
                         "Zea mays DK105 (GCA_003709335.1)" = "ZmaysDK105",
                         "Zea mays Oh7B (GCA_902166955.1)" = "ZmaysOh7B",
                         "Zea mays P39 (GCA_902166965.1)" = "ZmaysP39",
                         "Zea mays CML247 (GCA_902166975.1)" = "ZmaysCML247",
                         "Zea mays CML277 (GCA_902166985.1)" = "ZmaysCML277",
                         "Zea mays M162W (GCA_902166995.1)" = "ZmaysM162W",
                         "Zea mays NC350 (GCA_902167005.1)" = "ZmaysNC350",
                         "Zea mays OH43 (GCA_902167015.1)" = "ZmaysOH43",
                         "Zea mays CML333 (GCA_902167025.1)" = "ZmaysCML333",
                         "Zea mays NC358 (GCA_902167035.1)" = "ZmaysNC358",
                         "Zea mays Ky21 (GCA_902167045.1)" = "ZmaysKy21",
                         "Zea mays CML103 (GCA_902167055.1)" = "ZmaysCML103",
                         "Zea mays Il14H (GCA_902167065.1)" = "ZmaysIl14H",
                         "Zea mays B97 (GCA_902167075.1)" = "ZmaysB97",
                         "Zea mays HP301 (GCA_902167085.1)" = "ZmaysHP301",
                         "Zea mays CML322 (GCA_902167095.1)" = "ZmaysCML322",
                         "Zea mays Ki3 (GCA_902167105.1)" = "ZmaysKi3",
                         "Zea mays Ki11 (GCA_902167115.1)" = "ZmaysKi11",
                         "Zea mays CML69 (GCA_902167135.1)" = "ZmaysCML69",
                         "Zea mays CML228 (GCA_902167155.1)" = "ZmaysCML228",
                         "Zea mays Mo18W (GCA_902167165.1)" = "ZmaysMo18W",
                         "Zea mays M37W (GCA_902167175.1)" = "ZmaysM37W",
                         "Zea mays Tzi8 (GCA_902167185.1)" = "ZmaysTzi8",
                         "Zea mays Tx303 (GCA_902167205.1)" = "ZmaysTx303",
                         "Zea mays Ms71 (GCA_902167375.1)" = "ZmaysMs71",
                         "Zea mays LH244 (GCA_905067065.1)" = "ZmaysLH244",
                         "Zea mays W22 (GCA_001644905.2)" = "ZmaysW22",
                         "Zea mays EP1 (GCA_001984235.2)" = "ZmaysEP1",
                         "Zea mays F7 (GCA_001990705.1)" = "ZmaysF7"
                         )),
        textInput("gene_name", "Enter gene ID (only single query is supported)", value = "Zm00001d023520", placeholder = "Gene ID:", width = "400px"),
        br(),
        # HTML("<h5>Hit the button below <strong style='color: darkred'>(Initialize Visualization Panel)</strong> to load annotation data and draw example plot.</h5>"),
        actionButton("action_generate_plt", label = "Search Gene Information and Structure", style='height:40px; width:400px', class="btn-primary"),
        br(),
        br(),
        HTML("<p>The <strong>initial rendering progress may take longer</strong> due to the loading of R packages.<br> Please wait patiently and <strong style='color:darkred'>AVOID REPEATED SUBMISSION</strong>.</p>"),
        hr()
      )
    ),
      
    # Main panel for output plot
    fluidRow(
      br(),
      uiOutput("error_message"),
      br()
    ),
    plotOutput("plot"), # The main plot
    dataTableOutput("df_targetgene") # Simplified feature information dataframe
  )
) # shinyUI
