

getMarkersTable <- function(cell_type = "Alveolar_macrophage") {
  dt <-
    markers_table[cluster == cell_type,-c(which(colnames(markers_table) == "cluster")), with =
                    F]
  dt <- cbind(gene = dt$gene, dt[, 1:5])
  dt
}


# Volcano plot
plot_volcano <- function(de_table, cell_type, gene_name) {
  # extract cell type info
  pval_keyword <- "p_val"
  fc_keyword <- "avg_logFC"
  pval_column.name <- paste(cell_type, pval_keyword, sep = "__")
  fc_column.name <- paste(cell_type, fc_keyword, sep = "__")
  de.dt <-
    de_table[, c("rn", pval_column.name, fc_column.name), with = F]
  colnames(de.dt) <- c("Gene", "pvalue", "log2FoldChange")
  
  # remove NAs
  de.dt <- de.dt[complete.cases(de.dt)]
  
  # add colour for top 10
  de.dt <- de.dt[, colour := "none"]
  de.dt[Gene == gene_name, colour := "selected"]
  # de.dt <- de.dt[order(pvalue)]
  # de.dt[(1:10), colour := "top"]
  
  p <-
    ggplot(data = de.dt,
           aes(
             x = log2FoldChange,
             y = -log10(pvalue),
             col = colour,
             label = Gene
           )) +
    # , colour = threshold)) +
    geom_point() +
    labs(x = "Fold change (log2)", y = "-log10 p-value") +
    scale_color_manual(values = c(selected = "red", none = "black")) +
    guides(col = F) +
    ggtitle(cell_type)
  p +
    geom_point(
      data = de.dt[de.dt$colour == "selected",],
      aes(x = log2FoldChange, y = -log10(pvalue)),
      colour = "red",
      size = 2
    ) +
    geom_text_repel(data = de.dt[de.dt$colour == "selected",],
                    aes(label = Gene),
                    colour = "red",
                    size = 5)
}

# Protein bulk age DE boxplot
genBoxplot_protein <- function(protein = "Bpifa1") {
  if (length(intersect(protein, rownames(protein_bulk))) == 0)
    return(emptyPlot())
  
  expression <- log(protein_bulk[protein, ])
  
  dt <-
    data.frame(expression, grouping = c(rep("24m", 4), rep("3m", 4)))
  
  ggplot(dt, aes(
    factor(grouping, levels = c("3m", "24m")),
    expression,
    col = grouping,
    fill = grouping
  )) +
    geom_boxplot() + geom_jitter(colour = "black") +
    scale_color_manual(values = c(`3m` = "blue", `24m` = "red")) +
    scale_fill_manual(values = c(`3m` = "white", `24m` = "white")) +
    xlab("") + ylab("MS intensity") + ggtitle(protein)
}

# Dot plot
dotPlot <- function (gene_name = "Scgb1a1") {
  # Defaults
  cols.use = c("lightgrey", "blue")
  plot.legend = FALSE
  do.return = FALSE
  x.lab.rot = FALSE
  
  scale.func <- switch(EXPR = "radius",
                       'size' = scale_size,
                       'radius' = scale_radius,
                       stop("'scale.by' must be either 'size' or 'radius'"))
  
  MinMax <- function (data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    return(data2)
  }
  
  PercentAbove <- function (x, threshold) {
    return(length(x = x[x > threshold]) / length(x = x))
  }
  
  # Load gene expression
  expression <-
    h5read("data/AgingData.h5", name = as.character(gene_name))
  expression <- (expression - mean(expression)) / sd(expression)
  
  # remove cell types from cell info
  
  data.to.plot <- data.frame(expression)
  colnames(x = data.to.plot) <- gene_name
  data.to.plot$id <- cell_info$celltype
  data.to.plot <-
    data.to.plot %>% gather(key = genes.plot, value = expression,-c(id))
  data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>%
    dplyr::summarize(avg.exp = mean(expm1(x = expression)),
                     pct.exp = PercentAbove(x = expression,  threshold = 0))
  data.to.plot <-
    data.to.plot %>% ungroup() %>% group_by(genes.plot) %>%
    mutate(avg.exp.scale = scale(x = avg.exp)) %>% mutate(avg.exp.scale = MinMax(
      data = avg.exp.scale,
      max = 2.5,
      min = -2.5
    ))
  data.to.plot$pct.exp[data.to.plot$pct.exp < 0] <- NA
  data.to.plot <- as.data.frame(data.to.plot)
  colnames(data.to.plot) <-
    c("Cell_type", "Gene", "AvgExpr", "PctExpressed", "AvgRelExpr")
  
  bad <-
    c("red_blood_cells",
      "Gamma-Delta_T_cells",
      "low_quality_cells")
  data.to.plot <-
    data.to.plot[-match(bad, as.character(data.to.plot$Cell_type)),]
  data.to.plot <-
    data.to.plot[-which(is.na(data.to.plot$Cell_type)),]
  
  celltype_order <-
    rev(
      c(
        "Alveolar_macrophage",
        "Mki67+_proliferating_cells",
        "Natural_Killer_cells",
        "Plasma_cells",
        "B_cells",
        "Cd4+_T_cells",
        "CD8+_T_cells",
        "Interstitial_macrophages",
        "non-classical_monocyte_(Ly6c2-)",
        "classical_monocyte_(Ly6c2+)",
        "Cd103+/Cd11b-_dendritic_cells",
        "CD209+/Cd11b+_dendritic_cells",
        "Ccl17+/Cd103-/Cd11b-_dendritic_cells",
        "Megakaryocytes",
        "Neutrophils",
        "Eosinophils",
        "Fn1+_macrophage",
        "lymphatic_endothelial_cells",
        "Vcam1+_endothelial_cells",
        "vascular_endothelial_cells",
        "Capillary_endothelial_cells",
        "Mesothelial_cells",
        "Smooth_muscle_cells",
        "Interstitial_Fibroblast",
        "Lipofibroblast",
        "Type1_pneumocytes",
        "Type_2_pneumocytes",
        "Ciliated_cells",
        "Club_cells",
        "Goblet_cells"
      )
    )
  data.to.plot$Cell_type <-
    factor(data.to.plot$Cell_type, levels = celltype_order)
  
  p <-
    ggplot(data = data.to.plot, mapping = aes(x = Gene, y = Cell_type)) +
    geom_point(mapping = aes(size = PctExpressed, color = AvgRelExpr)) +
    scale.func(range = c(0, 10), limits = c(NA, NA)) +
    theme(
      axis.text.y = element_text(size = 13),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.position = c(0.75, 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  p
}


# mRNA boxplot
genBoxplot <-
  function(gene_name = "Scd1",
           cell_type = "Type_2_pneumocytes") {
    expression <-
      h5read("data/AgingData.h5", name = as.character(gene_name))
    
    dt <-
      cbind(expression = expression, cell_info[, .(grouping, celltype)])
    dt <- dt[celltype == cell_type]
    
    ggplot(dt, aes(
      factor(grouping, levels = c("3m", "24m")),
      expression,
      col = grouping,
      fill = grouping
    )) + geom_violin()  + geom_jitter(colour = "black") +
      scale_color_manual(values = c(`3m` = "blue", `24m` = "red")) +
      scale_fill_manual(values = c(`3m` = "blue", `24m` = "red")) +
      xlab("") + ylab("UMI counts [log2]") + ggtitle(gene_name)
  }

# tSNE plot
# genTSNEplot <- function(gene = "Ear2"){
#   expr <-
#     h5read("data/AgingData.h5", name = as.character(gene))
#
#   if(sum(expr) == 0) return(emptyPlot())
#
#   expr <- (expr - mean(expr)) / sd(expr)
#   expr.min <- quantile(expr, 0.01)
#   expr.max <- quantile(expr, 0.99)
#   if(expr.min == expr.max) expr <- h5read("data/AgingData.h5", name = as.character(gene))
#   #expr[which(expr > expr.max)] <- expr.max
#   #expr[which(expr < expr.min)] <- expr.min
#   farben <- color.scale(expr, extremes = c("grey", "darkblue"), alpha = 0.5)
#   plot(
#     tsne_coord,
#     col = farben,
#     pch = 19,
#     main = gene,
#     cex = 0.6)
# }

# tSNE plot
genTSNEplot <- function(gene_name = 'Frem1') {
  gene <-
    h5read(expression.file, name = gene_name)
  gene.min <- quantile(gene, 0.01)
  gene.max <- quantile(gene, 0.99)
  gene[which(gene > gene.max)] <- gene.max
  gene[which(gene < gene.min)] <- gene.min
  H5close()
  dt <- cbind(tsne_coord, expression = gene)
  
  if (all(gene == 0)) {
    high = "grey"
    
  }
  ggplot(dt) + geom_point(aes(tSNE_1, tSNE_2, col = gene), alpha = .5) +
    guides(col = F) +
    ggtitle(gene_name) +    scale_color_continuous(low = "grey", high = high)
  
}

# Solubility plot
genLinePlot <- function(protein = "Frem1") {
  if (length(intersect(protein, rownames(protein_fractions))) == 0)
    return(emptyPlot())
  
  age <- c(rep("young", 16), rep("old", 16))
  fractions <-
    c(
      rep("FR1", 4),
      rep("FR2", 4),
      rep("FR3", 4),
      rep("ECM", 4),
      rep("FR1", 4),
      rep("FR2", 4),
      rep("FR3", 4),
      rep("ECM", 4)
    )
  fractions <- factor(fractions, c("FR1", "FR2", "FR3", "ECM"))
  
  expr_tmp <- protein_fractions[protein,]
  expr_tmp <- log2(expr_tmp)
  means <-
    c(rep(mean(expr_tmp[1:16], na.rm = T), 16), rep(mean(expr_tmp[17:32], na.rm = T), 16))
  expr_tmp <- expr_tmp - means
  data <- data.frame(expression = expr_tmp, age, fractions)
  
  res <- summary(aov(expression ~ age * fractions, data = data))
  pval <- signif(res[[1]]$`Pr(>F)`[3], 2)
  title <- paste(protein, '(ANOVA interaction P:', pval, ")")
  
  agg <-
    ddply(data, .(age, fractions), function(x)
      c(
        mean = mean(x$expression, na.rm = T),
        se = sd(x$expression, na.rm = T) / sqrt(length(x$expression))
      ))
  agg$lower <- agg$mean + agg$se
  agg$upper <- agg$mean - agg$se
  
  ggplot(agg, aes(y = mean, x = fractions, colour = age)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .3) +
    geom_point() + ggtitle(title) +
    stat_summary(fun.y = mean,
                 geom = "smooth",
                 aes(group = age),
                 lwd = 1) +
    scale_color_manual(values = c(old = "red", young = "blue")) +
    ylab("Normalized MS-Intensity") + xlab("") +
    geom_hline(yintercept = 0, lty = 2) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 2), "cm"))
}

# Some meta functions
RNA_panel <- function(gene, celltype) {
  boxplot_rna <- genBoxplot(gene_name = gene, cell_type = celltype)
  dotplot_rna <- dotPlot(gene_name = gene)
  vulcano_rna <-
    plot_volcano(de_table = de_table,
                 cell_type = celltype,
                 gene_name = gene)
  
  lay <- rbind(c(1, 2, 3))
  grid.arrange(dotplot_rna, vulcano_rna, boxplot_rna, layout_matrix = lay)
}

Solubility_panel <- function(gene) {
  dotplot_rna <- dotPlot(gene_name = gene)
  solubility_prot <- genLinePlot(gene)
  
  lay <- rbind(c(1, 2, 2))
  grid.arrange(dotplot_rna, solubility_prot, layout_matrix = lay)
}

emptyPlot <- function() {
  df <- data.frame(x = 5, y = 5, text = "Not detected")
  p<- ggplot(df, aes(x, y, label = text)) +
    geom_point(col = "white") + xlim(0, 10) + ylim(0, 10) + geom_text()+ theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major= element_blank()
    ) +
    theme(plot.margin = unit(c(2, 2, 2, 2), "cm")) 
  class(p)[4] <- "empty_plot"
  p
}