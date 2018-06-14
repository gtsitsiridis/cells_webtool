# Define UI for application that plots random distributions
shinyUI(tagList(
  includeCSS("www/style.css"),
  dashboardPage(
    dashboardHeader(title = "Mouse lung aging atlas"
                    # , tags$li(class="dropdown", downloadButton(label ="Download", outputId = "download_plots_button"))
                    ),
    dashboardSidebar(
      sidebarMenu(
        id = "tabs",
        # menuItem("Overview",
        #          tabName = "overview"),
        menuItem("Lung cell atlas", tabName = "celltype_tab"),
        menuItem("Lung aging protein", tabName = "solubility_tab"),
        menuItem("Lung aging mRNA", tabName = "mRNA_tab")
        
      ),
      conditionalPanel("input.tabs != 'overview'",
                       uiOutput("gene_selector")),
      conditionalPanel("input.tabs == 'mRNA_tab'",
                       uiOutput("cell_type_selector"))
    ),
    dashboardBody(
      htmlOutput("help"),
      tabItems(
      # First tab content
      # tabItem(tabName = "overview",
      #         htmlOutput("description")),
      tabItem(tabName = "celltype_tab",
              fluidRow(
                box(
                  collapsible = TRUE,
                  width = 4,
                  plotOutput("dotplot_1", height = "600px")
                ),
                box(
                  collapsible = TRUE,
                  width = 4,
                  plotOutput("distplot", height = "600px")
                ),
                box(
                  collapsible = TRUE,
                  width = 4,
                  DT::dataTableOutput("markers_table", height = "600px")
                )
              )),
      tabItem(tabName = "solubility_tab",
              fluidRow(
                box(
                  collapsible = TRUE,
                  width = 4,
                  plotOutput("dotplot_2", height = "600px")
                ),
                box(
                  collapsible = TRUE,
                  width = 3,
                  plotOutput("protein_violinplot", height = "600px")
                ),
                box(
                  collapsible = TRUE,
                  width = 5,
                  plotOutput("solubility", height = "600px")
                )
              )),
      
      tabItem(tabName = "mRNA_tab",
              fluidRow(
                box(
                  collapsible = TRUE,
                  plotOutput("dotplot_3", height = "600px"),
                  width = 4
                ),
                box(
                  collapsible = TRUE,
                  plotlyOutput("gene_volcano", height = "600px"),
                  width = 5
                ),
                box(
                  collapsible = TRUE,
                  plotOutput("gene_violinplot", height = "600px"),
                  width = 3
                )
              ))
    ))
  )
))