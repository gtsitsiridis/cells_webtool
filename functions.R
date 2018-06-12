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
  de.dt[Gene == gene_name, colour :="selected"]
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
    guides(col = F)
  p + aes(
    x = log2FoldChange,
    y = -log10(pvalue),
    col = colour,
    label = Gene
  )
}

# genePlot <-
#   function (gene_name,
#             expression,
#             cell_info,
#             cell_type,
#             ident.include = NULL,
#             nCol = NULL,
#             do.sort = FALSE,
#             y.max = NULL,
#             same.y.lims = FALSE,
#             size.x.use = 16,
#             size.y.use = 16,
#             size.title.use = 20,
#             adjust.use = 1,
#             point.size.use = 1,
#             cols.use = NULL,
#             group.by = NULL,
#             y.log = FALSE,
#             x.lab.rot = FALSE,
#             y.lab.rot = FALSE,
#             legend.position = "right",
#             single.legend = TRUE,
#             remove.legend = FALSE,
#             do.return = FALSE,
#             return.plotlist = FALSE)
#   {
#     data.use <- data.frame(expression)
#     colnames(data.use) <- gene_name
#     cells.to.include <- which(cell_info$celltype == cell_type)
#     data.use <- data.use[cells.to.include, , drop = FALSE]
#     ident.use <-
#       factor(cell_info$grouping, levels = c("3m", "24m"))[cells.to.include]
#     gene.names <- gene_name
#     remove.legend <- TRUE
#
#     y.max <- max(data.use)
#
#     cols.use <- c("blue", "red")
#     x <- gene_name
#     Seurat:::SingleVlnPlot(
#       feature = x,
#       data = data.use[,
#                       x, drop = FALSE],
#       cell.ident = ident.use,
#       do.sort = do.sort,
#       y.max = y.max,
#       size.x.use = size.x.use,
#       size.y.use = size.y.use,
#       size.title.use = size.title.use,
#       adjust.use = adjust.use,
#       point.size.use = point.size.use,
#       cols.use = cols.use,
#       gene.names = gene.names,
#       y.log = y.log,
#       x.lab.rot = x.lab.rot,
#       y.lab.rot = y.lab.rot,
#       legend.position = legend.position,
#       remove.legend = remove.legend
#     )
#   }

# Seurat::DotPlot
# dotPlot(gene_name, expression, cell_info)
dotPlot <-
  function (gene_name,
            cell_info)
  {
    # Defaults
    cols.use = c("lightgrey", "blue")
    col.min = -2.5
    col.max = 2.5
    dot.min = 0
    dot.scale = 6
    scale.by = "radius"
    scale.min = NA
    scale.max = NA
    plot.legend = FALSE
    do.return = FALSE
    x.lab.rot = FALSE
    scale.func <- switch(EXPR = scale.by,
                         'size' = scale_size,
                         'radius' = scale_radius,
                         stop("'scale.by' must be either 'size' or 'radius'"))
    # Load gene expression
    expression <-
      h5read("data/scaledData.h5", name = as.character(gene_name))
    
    # remove cell types from cell info
    inds <- which(
      cell_info$celltype %in% c(
        "red_blood_cells",
        "Gamma-Delta_T_cells",
        "low_quality_cells "
      )
    )
    inds <- c(inds, which(is.na(cell_info$celltype)))
    cell_info <- cell_info[-inds]
    expression <- expression[-inds]
    
    data.to.plot <- data.frame(expression)
    colnames(x = data.to.plot) <- gene_name
    data.to.plot$id <- cell_info$celltype
    data.to.plot <- data.to.plot %>% gather(key = genes.plot,
                                            value = expression, -c(id))
    data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>%
      dplyr::summarize(
        avg.exp = mean(expm1(x = expression)),
        pct.exp = Seurat:::PercentAbove(x = expression,
                                        threshold = 0)
      )
    data.to.plot <-
      data.to.plot %>% ungroup() %>% group_by(genes.plot) %>%
      mutate(avg.exp.scale = scale(x = avg.exp)) %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale,
                                                                                   max = col.max, min = col.min))
    data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
    data.to.plot <- as.data.frame(data.to.plot)
    
    p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot,
                                                   y = id)) + geom_point(mapping = aes(size = pct.exp,
                                                                                       color = as.numeric(avg.exp.scale))) + scale.func(range = c(0, dot.scale),
                                                                                                                                        limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(),
                                                                                                                                                                                  axis.title.y = element_blank())
    p
  }

genBoxplot <-
  function(gene_name = "Scd1",
           cell_type = "Type_2_pneumocytes") {
    expression <-
      h5read("data/scaledData.h5", name = as.character(gene_name))
    
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
      xlab("") + ylab("Scaled expression") + ggtitle(gene_name)
  }


genLinePlot <- function(protein) {
  library(readxl)
  library(plyr)
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
  
  expr_tmp <- protein_fractions[protein, ]
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
    geom_hline(yintercept = 0, lty = 2)
}
