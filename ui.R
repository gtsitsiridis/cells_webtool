# Define UI for application that plots random distributions
shinyUI(tagList(
  includeCSS("www/style.css"),
  navbarPage(
    header= fluidRow(class = "control-bar",
                      
                     column(
                       4,
                       selectInput("gene", "Query gene/protein:", genes), offset = 2
                     ),
                     column(
                       4,
                       selectInput("cell_type", "Query cell type:", cell_types)
                     )
                     ),
    
    # Application title
    title = "Mouse lung aging atlas",
    
    tabPanel(
      "Overview",
      # Show a plot of the generated distribution
      column(6,offset = 3, plotOutput("tab1_distPlot", height = 250))
    ),
    
    tabPanel(
      "Solubility",
        fluidRow(# Show a plot of the generated distribution
          column(6, plotOutput("tab3_dotplot")),
          column(
            6, plotOutput("tab3_solubilityPlot", height = 250)
          ))
      ),
      tabPanel(
        "Age DE RNA",
        fluidRow(
          # Show a plot of the generated distribution
          column(4, plotOutput("tab2_dotplot")),
          column(4, plotlyOutput("tab2_volcanoPlot")),
          column(4, plotOutput("tab2_boxplot"))
        )
      )
    )
  )
)
