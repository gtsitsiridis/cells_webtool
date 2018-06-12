shinyServer(function(input, output, session) {
  output$tab1_distPlot <- renderPlot(width = 500, height = 500, {
    gene <-
      h5read("data/scaledData.h5", name = as.character(input$gene))
    gene.min <- quantile(gene, 0.01)
    gene.max <- quantile(gene, 0.99)
    gene[which(gene > gene.max)] <- gene.max
    gene[which(gene < gene.min)] <- gene.min
    farben <-
      color.scale(gene,
                  extremes = c("grey", "darkblue"),
                  alpha = 0.5)
    try(plot(
      tsne_coord,
      col = farben,
      pch = 19,
      main = input$gene,
      cex = 0.6
    ),
    silent = T)
  })
  
  output$tab2_volcanoPlot <- renderPlotly({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      cell_type <- input$cell_type
      gene_name <- input$gene
      try({
        p <- plot_volcano(de_table, cell_type, gene = gene_name)
        p <-
          ggplotly(p,
                   tooltip = c("Gene", "-log10(pvalue)", "log2FoldChange"))
      }, silent = T)
    })
  })
  
  output$tab2_boxplot <- renderPlot({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      gene_name <- input$gene
      expression <-
        h5read("data/scaledData.h5", name = as.character(gene_name))
      cell_type <- input$cell_type
      
      try(genBoxplot(gene_name, cell_type), silent = T)
    })
  })
  
  output$tab3_dotplot <- renderPlot({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      gene_name <- input$gene
      try(dotPlot(gene_name, cell_info), silent = T)
    })
  })
  output$tab2_dotplot <- renderPlot({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      gene_name <- input$gene
      try(dotPlot(gene_name, cell_info), silent = T)
    })
  })
  output$tab3_solubilityPlot <- renderPlot({
    withProgress(session = session, value = 0.5, {
      setProgress(message = "Calculation in progress")
      protein <- input$gene
      try(genLinePlot(protein), silent = T)
    })
  })
  
})
