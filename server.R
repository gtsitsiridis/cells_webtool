shinyServer(function(input, output, session) {
  ### Descriptions
  help <- list()
  
  help[["celltype_tab"]] <-
    HTML('<center><strong><p>celltype_tab help</strong></p></center>')
  help[["solubility_tab"]] <-
    HTML('<center><strong><p>solubility_tab help</strong></p></center>')
  help[["mRNA_tab"]] <-
    HTML('<center><strong><p>mRNA_tab help</strong></p></center>')
  help[["overview"]] <-
    HTML('')
  help[["enrichment_tab"]] <-
    HTML('<center><strong><p>enrichment_tab help</strong></p></center>')
  
  output$description <-
    renderUI(HTML("<center><h1>Add description of the tool</h1></center>"))
  
  output$help <- renderUI({
    tab <- input$tabs
    help[[tab]]
  })
  
  ### Helper functions
  check_save <- function(plot) {
    # Check if exists
    if (class(plot) == "try-error") {
      plot  <-
        emptyPlot()
    } else{
      # Save plot
      cl <- class(plot)[3]
      plots[[cl]] <- plot
    }
    plot
  }
  
  ### Reactive values
  plots <- reactiveValues(
    distplot = NULL,
    dotplot = NULL,
    protein_volcano = NULL,
    gene_volcano = NULL,
    solubility = NULL,
    gene_violinplot = NULL,
    protein_violinplot = NULL
  )
  
  values <- reactiveValues(gene = genes[1],
                           cell_type = cell_types[1])
  
  ### Pass input to values
  # observeEvent(values$gene, {
  #   new_gene_name <- values$gene
  #   updateSelectInput(session,
  #                     "gene", "Query gene/protein:", genes, selected = new_gene_name)
  # })
  observeEvent(input$gene, {
    new_gene_name <- input$gene
    values$gene <- new_gene_name
  })
  # observeEvent(values$cell_type, {
  #   new_cell_type <- values$cell_type
  #   updateSelectInput(session,
  #                     "cell_type",
  #                     "Query cell type:",
  #                     cell_types,
  #                     selected = new_cell_type)
  # })
  observeEvent(input$cell_type, {
    new_cell_type <- input$cell_type
    values$cell_type <- new_cell_type
  })

  ### Define gene and cell type selectors
  output$cell_type_selector <- renderUI({
    selectInput("cell_type", "Query cell type:", cell_types)
  })
  output$gene_selector <- renderUI({
    selectInput("gene", "Query gene/protein:", genes)
  }) 
  output$enrichment_type_selector <- renderUI({
    selectInput("enrichment_type", "Query enrichment type:", enrichment_types)
  })
  
  ### Create plots
  output$distplot <- renderPlot({
    gene_name <- values$gene
    p <- try(genTSNEplot(gene_name),
             silent = T)
    class(p)[3] <- "distplot"
    check_save(p)
  })
  
  output$gene_volcano <- renderPlotly({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      cell_type <- values$cell_type
      gene_name <- values$gene
      p <-
        try(p <- plot_volcano(de_table = gene_de_table,
                              cell_type = cell_type,
                              gene_name = gene_name)
            ,
            silent = T)
      class(p)[3] <- "gene_volcano"
      ggplotly(
        check_save(p),
        tooltip = c("Gene", "-log10(pvalue)", "log2FoldChange"),
        source = "gene_volcano"
      )
    })
  })
  
  output$gene_violinplot <- renderPlot({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      gene_name <- values$gene
      cell_type <- values$cell_type
      p <- try(genBoxplot(gene_name, cell_type), silent = T)
      class(p)[3] <- "gene_violinplot"
      check_save(p)
    })
  })
  output$protein_violinplot <- renderPlot({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      gene_name <- values$gene
      p <- try(genBoxplot_protein(gene_name), silent = T)
      class(p)[3] <- "protein_violinplot"
      check_save(p)
    })
  })
  output$dotplot_1 <- renderPlot({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      gene_name <- values$gene
      cell_type <- values$cell_type
      p <- try(dotPlot(gene_name), silent = T)
      class(p)[3] <- "dotplot"
      check_save(p)
    })
  })
  output$dotplot_2 <- renderPlot({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      gene_name <- values$gene
      cell_type <- values$cell_type
      p <- try(dotPlot(gene_name), silent = T)
      class(p)[3] <- "dotplot"
      check_save(p)
    })
  })
  output$dotplot_3 <- renderPlot({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      gene_name <- values$gene
      cell_type <- values$cell_type
      p <- try(dotPlot(gene_name), silent = T)
      class(p)[3] <- "dotplot"
      check_save(p)
    })
  })
  output$solubility <- renderPlot({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      protein <- values$gene
      cell_type <- values$cell_type
      p <- try(genLinePlot(protein), silent = T)
      class(p)[3] <- "solubility"
      check_save(p)
    })
  })
  
  output$enrichment_barplot <- renderPlot({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      cell_type <- values$cell_type
      enrichment_type <- input$enrichment_type
      if(is.null(enrichment_type)){
        return()
      }
      p <- try(enrichmentBarPlot(cell_type, enrichment_type), silent = F)
      print(p)
      class(p)[3] <- "enrichment_barplot"
      check_save(p)
    })
  })
  
  output$markers_table <- DT::renderDataTable({
    cell_type <- values$cell_type
    # gene <- values$gene
    dt <- getMarkersTable(cell_type)
    DT::datatable(
      dt,
      extensions = 'Buttons',
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        scrollY = "400px",
        searchHighlight = T,
        dom = '<"top"Bf>rt<"bottom"lip><"clear">',
        buttons = list(
          'print',
          list(
            extend =  "csv",
            title = "file"
          ),
          list(
            extend =  "pdf",
            title = "file"
          )
        )
      ),
      rownames = FALSE,
      selection = list(
        mode = 'single',
        target = 'row'
        # selected = which(dt$gene == gene)
      )
    )
  })
  output$enrichment_table <- DT::renderDataTable({
    cell_type <- values$cell_type
    enrichment_type <- input$enrichment_type
    if(is.null(enrichment_type)){
      return()
    }
    dt <- getEnrichmentTable(cell_type, enrichment_type)
    DT::datatable(
      dt,
      extensions = 'Buttons',
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        scrollY = "500px",
        searchHighlight = T,
        dom = '<"top"Bf>rt<"bottom"lip><"clear">',
        buttons = list(
          'print',
          list(
            extend =  "csv",
            title = "file.csv"
          ),
          list(
            extend =  "pdf",
            title = "file.pdf"
          )
        )
      ),
      rownames = FALSE,
      selection = "none"
    )
  })
  
  ### Extra features
  # Download plots
  output$download_plots_button <-
    downloadHandler(
      filename = function() {
        isolate(tab <- input$tabs)
        paste0(gsub("\\s", "_", tab), "_plots.zip")
      },
      content = function(file) {
        isolate(tab <- input$tabs)
        if (tab == "celltype_tab") {
          plot_names <- c("dotplot", "distplot")
        } else if (tab == "solubility_tab") {
          plot_names <- c("dotplot", "solubility", "protein_violinplot")
        } else if (tab == "mRNA_tab") {
          plot_names <- c("dotplot", "gene_volcano", "gene_violinplot")
        }else if (tab == "enrichment_tab") {
          plot_names <- c("enrichment_barplot")
        }
        files <- sapply(plot_names, function(x) {
          p <- plots[[x]]
          file_name <- paste0(x, ".pdf")
          if (is.null(p)) {
            return(NA)
          }
          ggsave(file = file_name, plot = p)
          return(file_name)
        })
        files <- files[!is.na(files)]
        zip(file, files)
      }
    )
  
  #deal with selection from marker's table
  observeEvent(input$markers_table_rows_selected, {
    row_selected <- input$markers_table_rows_selected
    
    isolate(cell_type <- values$cell_type)
    dt <-
      markers_table[cluster == cell_type, -c(which(colnames(markers_table) == "cluster")), with =
                      F]
    
    new_gene_name <- dt[row_selected, gene]
    values$gene <- new_gene_name
  })
  
  
  # add clicking ability to gene_volcano
  observeEvent(event_data("plotly_click", source = "gene_volcano"), {
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      event <- event_data("plotly_click", source = "gene_volcano")
      if (is.null(event)) {
        return(NULL)
      }
      selected <- event$curveNumber == 1
      dt <- isolate(copy(plots$gene_volcano$data))
      dt[, selection := (colour != "none")]
      dt <- dt[selection == selected]
      row <- event$pointNumber + 1
      new_gene_name <- dt[row, Gene]
      values$gene <- new_gene_name
    })
  })
})