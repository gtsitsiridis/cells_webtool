# Define UI for application that plots random distributions
spinner <- function(ui_element){
  withSpinner(ui_element, type=2, color.background="white")
}

shinyUI(tagList(
  includeCSS("www/style.css"),
  dashboardPage(
    dashboardHeader(title = "Mouse lung aging atlas"
                    , tags$li(class="dropdown",

                              conditionalPanel(condition= "input.tabs != 'overview'", downloadButton(label ="Download plots", class='btn-primary',outputId = "download_plots_button")))
                    ),
    dashboardSidebar(
      sidebarMenu(
        id = "tabs",
        # menuItem("Overview",
        #          tabName = "overview"),
        menuItem("Lung cell type signatures", tabName = "celltype_tab"),
        menuItem("Lung aging protein", tabName = "solubility_tab"),
        menuItem("Lung aging mRNA", tabName = "mRNA_tab")
        
      ),
      
      conditionalPanel("input.tabs != 'overview'",
                       uiOutput("gene_selector")),
      conditionalPanel("input.tabs == 'mRNA_tab' | input.tabs == 'celltype_tab'",
                       uiOutput("cell_type_selector"))
    ),
    dashboardBody(
      htmlOutput("help"),
      tabItems(
      # First tab content
      # tabItem(tabName = "overview",
      #        HTML("<center><h1>Add description of the tool</h1></center>")),
      tabItem(tabName = "celltype_tab",
              fluidRow(
                box(
                  collapsible = TRUE,
                  width = 4,
                  spinner(plotOutput("dotplot_1", height = "600px"))
                )
                ,
                box(
                  collapsible = TRUE,
                  width = 4,
                  spinner( plotOutput("distplot", height = "600px"))
                ),
                box(
                  collapsible = TRUE,
                  width = 4,
                  spinner( DT::dataTableOutput("markers_table", height = "600px"))
                )
              )),
      tabItem(tabName = "solubility_tab",
              fluidRow(
                box(
                  collapsible = TRUE,
                  width = 4,
                  spinner( plotOutput("dotplot_2", height = "600px"))
                ),
                box(
                  collapsible = TRUE,
                  width = 3,
                  spinner( plotOutput("protein_violinplot", height = "600px"))
                ),
                box(
                  collapsible = TRUE,
                  width = 5,
                  spinner( plotOutput("solubility", height = "600px"))
                )
              )),
      
      tabItem(tabName = "mRNA_tab",
              fluidRow(
                box(
                  collapsible = TRUE,
                  spinner( plotOutput("dotplot_3", height = "600px")),
                  width = 4
                ),
                box(
                  collapsible = TRUE,
                  spinner( plotlyOutput("gene_volcano", height = "600px")),
                  width = 5
                ),
                box(
                  collapsible = TRUE,
                  spinner( plotOutput("gene_violinplot", height = "600px")),
                  width = 3
                )
              ))
    ))
  )
))