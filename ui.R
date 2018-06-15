# Define UI for application that plots random distributions
spinner <- function(ui_element) {
  withSpinner(ui_element, type = 2, color.background = "white")
}

shinyUI(tagList(
  includeCSS("www/style.css"),
  dashboardPage(
    skin = "black",
    dashboardHeader(
      titleWidth = 407,
      title = HTML("Lung Aging Atlas – Schiller and Theis labs"),
      tags$li(
        class = "dropdown",
        HTML(
          "<img src='Overview_logos.png' style='padding-top:10px;padding-right:10px;' height='70'/>"
        )
      )
    ),
    dashboardSidebar(
      width = 250,
      sidebarMenu(
        id = "tabs",
        # menuItem("Overview",
        #          tabName = "overview"),
        menuItem("Lung cell type signatures", tabName = "celltype_tab"),
        menuItem("Lung aging protein", tabName = "solubility_tab"),
        menuItem("Lung aging mRNA", tabName = "mRNA_tab"),
        menuItem("Lung aging – annotation enrichments", tabName = "enrichment_tab")
        
      ),
      
      conditionalPanel(
        "input.tabs == 'mRNA_tab' | input.tabs == 'celltype_tab'|input.tabs=='solubility_tab'",
        withSpinner(
          uiOutput("gene_selector"),
          type = 2,
          color.background = "#222d32"
        )
      ),
      conditionalPanel(
        "input.tabs == 'mRNA_tab' | input.tabs == 'celltype_tab'|input.tabs=='enrichment_tab'",
        uiOutput("cell_type_selector")
      ),
      conditionalPanel(
        "input.tabs=='enrichment_tab'",
        uiOutput("enrichment_type_selector")
      )
    ),
    dashboardBody(
      fluidRow(column(10, htmlOutput("help")), column(
        2,
        conditionalPanel(
          condition = "input.tabs == 'mRNA_tab' | input.tabs == 'celltype_tab'|input.tabs=='solubility_tab'|input.tabs=='enrichment_tab'",
          downloadButton(
            label = "Download plots",
            # class = 'btn-primary',
            outputId = "download_plots_button"
          )
        )
      )),
      
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
                    spinner(plotOutput("distplot", height = "600px"))
                  ),
                  box(
                    collapsible = TRUE,
                    width = 4,
                    spinner(DT::dataTableOutput("markers_table", height = "600px"))
                  )
                )),
        tabItem(tabName = "solubility_tab",
                fluidRow(
                  box(
                    collapsible = TRUE,
                    width = 4,
                    spinner(plotOutput("dotplot_2", height = "600px"))
                  ),
                  box(
                    collapsible = TRUE,
                    width = 3,
                    spinner(plotOutput("protein_violinplot", height = "600px"))
                  ),
                  box(
                    collapsible = TRUE,
                    width = 5,
                    spinner(plotOutput("solubility", height = "600px"))
                  )
                )),
        
        tabItem(tabName = "mRNA_tab",
                fluidRow(
                  box(
                    collapsible = TRUE,
                    spinner(plotOutput("dotplot_3", height = "600px")),
                    width = 4
                  ),
                  box(
                    collapsible = TRUE,
                    spinner(plotlyOutput("gene_volcano", height = "600px")),
                    width = 5
                  ),
                  box(
                    collapsible = TRUE,
                    spinner(plotOutput("gene_violinplot", height = "600px")),
                    width = 3
                  )
                )),
        tabItem(tabName = "enrichment_tab",
                fluidRow(
                  collapsible = TRUE,
                  box(width = 6, spinner(
                    DT::dataTableOutput("enrichment_table", height = "700px")
                  )),
                  box(width = 6, spinner(
                    plotOutput("enrichment_barplot", height = "700px")
                  ))
                ))
        
      )
    )
  )
))