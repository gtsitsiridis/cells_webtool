shinyServer(function(input, output, session) {
  ### Helper functions
  check_save <- function(plot) {
    # Check if exists
    if (class(plot) == "try-error") {
      plot  <-
        plot(
          0,
          col = 'white',
          xlim = c(-1, 1),
          ylim = c(-1, 1),
          yaxt = 'n',
          xaxt = 'n',
          ann = F
        )
      text(0, 0, 'Not detected')
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
    violinplot = NULL
  )
  
  values <- reactiveValues(gene = genes[1],
                           cell_type = cell_types[1])
  
  observeEvent(values$gene, {
    new_gene_name <- values$gene
    updateSelectInput(session,
                      "gene", "Query gene/protein:", genes, selected = new_gene_name)
  })
  observeEvent(input$gene, {
    new_gene_name <- input$gene
    values$gene <- new_gene_name
  })
  observeEvent(values$cell_type, {
    new_cell_type <- values$cell_type
    updateSelectInput(session,
                      "cell_type", "Query cell type:", cell_types,
                      selected = new_cell_type)
  })
  observeEvent(input$cell_type, {
    new_cell_type <- input$cell_type
    values$cell_type <- new_cell_type
  })
  output$description <- renderUI(HTML("<center><h1>Add description of the tool</h1></center>"))
  
  ### Create plots
  output$distplot <- renderPlot({
    gene_name <- values$gene
    p <- try(
      distplot(gene_name),
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
  
  output$violinplot <- renderPlot({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      gene_name <- values$gene
      cell_type <- values$cell_type
      p <- try(genBoxplot(gene_name, cell_type), silent = T)
      class(p)[3] <- "violinplot"
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
  output$protein_volcano <- renderPlot({
    p <- 1
    class(p) <- "try-error"
    class(p)[3] <- "protein_volcano"
    check_save(p)
  })
  
  ### Extra features
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