  # Define UI for application that plots random distributions
  shinyUI(tagList(
    includeCSS("www/style.css"),
    dashboardPage(
      dashboardHeader(title = "Mouse lung aging atlas"),
      dashboardSidebar(
        sidebarMenu(
          menuItem("Overview",
                   tabName = "overview"),
          menuItem("Lung cell atlas", tabName = "lca"),
          menuItem("Lung aging protein", tabName = "lap"),
          menuItem("Lung aging mRNA", tabName = "lam"),
          selectInput("gene", "Query gene/protein:", genes),
          selectInput("cell_type", "Query cell type:", cell_types)
        )
      ),
      dashboardBody(tabItems(
        # First tab content
        tabItem(tabName = "overview",
                htmlOutput("description")),
        tabItem(tabName = "lca",
                fluidRow(
                  box(
                    collapsible = TRUE,
                    width = 8,
                    plotOutput("distplot", height = "600px")
                  ),
                  box(
                    collapsible = TRUE,
                    width = 4,
                    plotOutput("dotplot_1", height = "600px")
                  )
                )),
        tabItem(tabName = "lap",
                fluidRow(
                  box(
                    collapsible = TRUE,
                    width = 8,
                    plotOutput("protein_volcano", height = "600px")
                  ),
                  box(
                    collapsible = TRUE,
                    width = 4,
                    plotOutput("dotplot_2", height = "600px")
                  )
                ),
                
                fluidRow(
                  box(
                    collapsible = TRUE,
                    width = 8,
                    plotOutput("solubility", height = "600px")
                  )
                )),
        tabItem(tabName = "lam",
                fluidRow(
                  box(
                    collapsible = TRUE,
                    plotlyOutput("gene_volcano", height = "600px"),
                    width = 8
                  ),
                  box(
                    collapsible = TRUE,
                    plotOutput("dotplot_3", height = "600px"),
                    width = 4
                  )
                  
                ),
                fluidRow(
                  box(
                    collapsible = TRUE,
                    plotOutput("violinplot", height = "600px"),
                    width = 8
                  )
                ))
      ))
    )
    
    # navbarPage(
    #   id = "atlas",
    #   header =
    #     fluidRow(
    #       class = "control-bar",
    #       column(4,
    #              selectInput("gene", "Query gene/protein:", genes), offset = 2),
    #       column(4,
    #              selectInput(
    #                "cell_type", "Query cell type:", cell_types
    #              ))
    #     ),
    #
    #   # Application title
    #   title = "Mouse lung aging atlas",
    #
    #   tabPanel("Overview",
    #            htmlOutput("description")),
    #   tabPanel("Lung cell atlas signatures",
    #            fluidRow(
    #              # Show a plot of the generated distribution
    #              column(6, plotOutput("dotplot_1", height = 250)),
    #              column(6, plotOutput("distplot"))
    #            )),
    #
    #   tabPanel("Lung aging protein",
    #            fluidRow(
    #              # Show a plot of the generated distribution
    #              column(4, plotOutput("dotplot_2")),
    #              column(4, plotOutput("protein_volcano", height = 250)),
    #              column(4, plotOutput("solubility"))
    #            )),
    #   tabPanel("Lung aging mRNA",
    #            fluidRow(
    #              # Show a plot of the generated distribution
    #              column(4, plotOutput("dotplot_3")),
    #              column(4, plotlyOutput("gene_volcano")),
    #              column(4, plotOutput("violinplot"))
    #            ))
    # )
  ))
