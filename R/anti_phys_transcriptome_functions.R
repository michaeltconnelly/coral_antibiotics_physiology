#anti_phys_transcriptome_functions.R
#author: "Mike Connelly"
#date: "05/07/2020"

# Create ggplot2 theme for journal submission
# ISME figure guidelines
# 1-column width: 85 mm, 2-column width: 175 mm
# All text should be sans-serif typeface, preferably Helvetica or Arial.
# Maximum text size is 7pt. Minimum text size is 5pt.
theme_ismej <- function(base_size = 6) {
  (theme_foundation(base_size=base_size)
   + theme(
     # 
     plot.background = element_rect(colour = NA),
     panel.background = element_blank(),
     panel.grid.major = element_blank(), 
     panel.grid.minor = element_blank(), 
     panel.border = element_rect(color = "black", fill = NA),
     # 
     plot.title = element_text(face = "plain", size = rel(1), hjust = 0),
     plot.subtitle = element_text(face = "plain", size = rel(0.8)),
     # 
     axis.title = element_text(face = "plain",size = rel(1)),
     axis.title.x = element_text(vjust = -2),
     axis.title.y = element_text(angle = 90, vjust = 2),
     text = element_text(),
     axis.text = element_text(size = rel(0.8)), 
     axis.text.x = element_text(angle = 0),
     axis.text.y = element_text(), 
     #axis.line = element_line(colour = "black"),
     axis.ticks = element_line(),
     # 
     legend.title = element_text(size = rel(1)),
     legend.text = element_text(size = rel(0.8), margin = margin(t = 2, b = 2, unit = "mm")),
     legend.key = element_rect(color = NA),
     legend.background = element_rect(fill = NA, colour = NA),
     legend.position = "right",
     legend.direction = "vertical",
     legend.spacing.x = unit(2, "mm"),
     legend.spacing.y = unit(0, "mm"),
     legend.key.size = unit(2, "mm"),
     # 
     plot.margin = unit(c(2,2,2,2), "mm"),
     # 
     strip.background = element_rect(color = "black", fill = "grey"),
     strip.text = element_text(face = "plain")
   ))
}
###

# EAPSI AXH Transcriptome Analysis Functions ---------------------------------------------------------------------------------
### PCA plot with custom PC axes----------------------------------------------------------------------------------
plotPCA.custom <-  function(object, intgroup="Treatment", ntop=500, returnData=FALSE, pcs = c(1,2))
{
  stopifnot(length(pcs) == 2)    ### added this to check number of PCs ####
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]]
  }
  # assemble the data for the plot
  ########## Here we just use the pcs object passed by the end user ####
  d <- data.frame(PC1=pca$x[,pcs[1]], PC2=pca$x[,pcs[2]], group=group, intgroup.df, name=colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  # extract loadings
}

### PCA plot formatted with aesthetics ---------------------------------------------------------------------------
ggPCA <- function(vsd, samples, condcolors, ntop = 500,  pclab = c(1,2)) {
  #
  PCAtmtdata <- plotPCA.custom(vsd, intgroup = c("Treatment", "Genotype"), ntop = ntop, returnData = TRUE,  pcs = c(pclab[1],pclab[2]))
  #set factor orders 
  PCAtmtdata$Genotype <- factor(PCAtmtdata$Genotype, levels = genotype_levels, ordered = TRUE)
  PCAtmtdata$Treatment <- factor(PCAtmtdata$Treatment, levels = treatment_levels, ordered = TRUE)
  #
  PCAtmtpercentVar <- round(100 * attr(PCAtmtdata, "percentVar"), 1)
  #
  PCAplot <-  PCAtmtdata %>% ggplot(aes(PC1,PC2)) +
    geom_point(aes(fill = Treatment), shape = 21, size = 2, alpha = 1, stroke = 0.5, color = "black", show.legend = TRUE) +
    stat_ellipse(aes(color = Treatment), type = "norm") +
    xlab(paste0( "PC", pclab[1], " (", PCAtmtpercentVar[pclab[1]], "%)")) + 
    ylab(paste0( "PC", pclab[2], " (", PCAtmtpercentVar[pclab[2]], "%)")) + 
    coord_fixed(1) + 
    scale_fill_manual(values=treatcolors[2:3], name="Treatment") +
    scale_color_manual(values=treatcolors[2:3], name="Treatment") +
    # scale_shape_manual(values=colshapes, name="Colony") +
    theme(legend.position = "right") +
    guides(color = guide_legend(override.aes = list(color = treatcolors[2:3], alpha = 1, stroke = 1)))
  #
PCAplot
}

### PCA plot formatted with modified aesthetics - treatment convex hulls ---------------------------------------------------------------------------
ggPCA_mod1 <- function(vsd, samples, condcolors, ntop = 500,  pclab = c(1,2)) {
  #
  PCAtmtdata <- plotPCA.custom(vsd, intgroup = c("Treatment", "Genotype"), ntop = ntop, returnData = TRUE,  pcs = c(pclab[1],pclab[2]))
  #set factor orders 
  PCAtmtdata$Genotype <- factor(PCAtmtdata$Genotype, levels = genotype_levels, ordered = TRUE)
  PCAtmtdata$Treatment <- factor(PCAtmtdata$Treatment, levels = treatment_levels, ordered = TRUE)
  #
  PCAtmtpercentVar <- round(100 * attr(PCAtmtdata, "percentVar"), 1)
  #
  PCAtmthull <- PCAtmtdata %>%
    group_by(Treatment) %>% 
    dplyr::slice(chull(PC1, PC2))
  #
  PCAtmtcent <- PCAtmtdata %>% 
    dplyr::group_by(Treatment) %>% 
    dplyr::summarise(c1 = mean(`PC1`), c2 = mean(`PC2`)) %>%    
    full_join(PCAtmtdata)
  #
  PCAtmtlab <- PCAtmtdata %>% 
    dplyr::group_by(Treatment) %>% 
    dplyr::summarise(c1 = mean(`PC1`), c2 = mean(`PC2`)) %>% 
    dplyr::mutate(Code = ifelse(`Treatment` == "Control", "C", "A"))
  #
  PCAplot <-  PCAtmtdata %>% ggplot(aes(PC1,PC2)) +
    # spider segments
    # geom_segment(data = PCAtmtcent, mapping = aes(x = `PC1`, y = `PC2`, xend = c1, yend = c2), lwd = 0.25, col = "dark grey") +
    # convex hull
    geom_polygon(data = PCAtmthull,
                 aes(fill = Treatment, color = Treatment),
                 alpha = 0.3,
                 show.legend = FALSE) + 
    # sample points
    geom_point(aes(fill = Treatment), shape = 21, size = 2, alpha = 1, stroke = 0.5, color = "black", show.legend = TRUE) +
    # treatment labels
    geom_label(data = PCAtmtlab, size = 2, aes(x = c1, y = c2, label = `Code`), box.padding = 0.15, alpha = 0.8, segment.alpha = 0) +
    xlab(paste0( "PC", pclab[1], " (", PCAtmtpercentVar[pclab[1]], "%)")) +
    ylab(paste0( "PC", pclab[2], " (", PCAtmtpercentVar[pclab[2]], "%)")) +
    coord_fixed(PCAtmtpercentVar[pclab[2]]/PCAtmtpercentVar[pclab[1]]) + 
    scale_fill_manual(values=treatcolors[2:3], name="Treatment") +
    scale_color_manual(values=treatcolors[2:3], name="Treatment") +
    # scale_shape_manual(values=colshapes, name="Colony") +
    # theme(legend.position = "none") +
    guides(color = guide_legend(override.aes = list(color = "black", fill = treatcolors[2:3], alpha = 1, stroke = 0.5)), fill = NULL)
  #
  PCAplot
}

### PCA plot formatted with modified aesthetics - spider lines ---------------------------------------------------------------------------
ggPCA_mod2 <- function(vsd, samples, condcolors, ntop = 500,  pclab = c(1,2)) {
  #
  PCAtmtdata <- plotPCA.custom(vsd, intgroup = c("Treatment", "Genotype"), ntop = ntop, returnData = TRUE,  pcs = c(pclab[1],pclab[2]))
  #set factor orders 
  PCAtmtdata$Genotype <- factor(PCAtmtdata$Genotype, levels = genotype_levels, ordered = TRUE)
  PCAtmtdata$Treatment <- factor(PCAtmtdata$Treatment, levels = treatment_levels, ordered = TRUE)
  #
  PCAtmtpercentVar <- round(100 * attr(PCAtmtdata, "percentVar"), 1)
  #
  #
  PCAtmthull <- PCAtmtdata %>%
    # dplyr::select(-group, -Treatment, -name) %>% 
    dplyr::group_by(Genotype) %>% 
    dplyr::slice(chull(PC1, PC2))
  #
  PCAtmtcent <- PCAtmtdata %>% 
    dplyr::group_by(Genotype) %>% 
    dplyr::summarise(c1 = mean(`PC1`), c2 = mean(`PC2`)) %>%    
    full_join(PCAtmtdata)
  #
  PCAtmtlab <- PCAtmtdata %>% 
    dplyr::group_by(Genotype) %>% 
    dplyr::summarise(c1 = mean(`PC1`), c2 = mean(`PC2`))
  #
  PCAplot <-  PCAtmtdata %>% ggplot(aes(PC1,PC2)) +
    # convex hull
    # geom_polygon(data = PCAtmthull, aes(Group = `Genotype`),
                 # color = "black", fill = "white",
                 # alpha = 0.3,
                 # show.legend = FALSE) + 
    # spider segments
    geom_segment(data = PCAtmtcent, mapping = aes(x = `PC1`, y = `PC2`, xend = c1, yend = c2), lwd = 0.5, col = "darkgrey") +
    # genotype centroid points
    # geom_point(data = PCAtmtcent, size = 1, aes(x = c1, y = c2), fill = "black", color = "black", show.legend = FALSE) + 
    # sample points
    geom_point(aes(fill = Treatment), shape = 21, size = 2, alpha = 1, stroke = 0.5, color = "black", show.legend = TRUE) +
    # Genotype centroid labels
    geom_label(data = PCAtmtlab, size = 1.5, aes(x = c1, y = c2, label = `Genotype`), alpha = 0.8) +
    # stat_ellipse(aes(Group = `Genotype`), color = "black", type = "norm") +
    xlab(paste0( "PC", pclab[1], " (", PCAtmtpercentVar[pclab[1]], "%)")) + 
    ylab(paste0( "PC", pclab[2], " (", PCAtmtpercentVar[pclab[2]], "%)")) +
    # ylab(NULL) +
    coord_fixed(PCAtmtpercentVar[pclab[2]]/PCAtmtpercentVar[pclab[1]]) + 
    scale_fill_manual(values=treatcolors[2:3], name="Treatment") +
    scale_color_manual(values=treatcolors[2:3], name="Treatment") +
    # scale_shape_manual(values=colshapes, name="Colony") +
    # theme(legend.position = "none") +
    guides(color = guide_legend(override.aes = list(color = "black", fill = treatcolors[2:3], alpha = 1, stroke = 0.5)), fill = NULL)
  #
  PCAplot
  # return(PCAtmthull)
}

### PCoA plot formatted with aesthetics ------------------------------------------------------------------------------
ggPCoA <- function(mds) {
  #
  # Plot with spiders
  PCoAplot <- ggplot(mds, aes(fill = Treatment)) +
    # geom_segment(mapping = aes(x = `1`, y = `2`, xend = c1, yend = c2), lwd = 0.25, col = "dark grey") +
    # treatment centroid points
    # geom_point(size = 3, aes(x = c1, y = c2, fill = Treatment), color = "black", shape = 21, stroke = 1, show.legend = TRUE) + #FALSE
    # sample points
    geom_point(size = 2, aes(x = `1`, y = `2`, fill = Treatment), shape = 21, size = 2, stroke = 0.5, color = "black", show.legend = TRUE) + #FALSE
    stat_ellipse(aes(x = `1`, y = `2`, color = Treatment), type = "norm") +
    scale_color_manual(values = treatcolors[2:3], name = "Treatment") +
    # scale_shape_manual(values = c(21,22), name = "Treatment") +
    scale_fill_manual(values = treatcolors[2:3], name = "Treatment") +
    labs(x = xlab, y = ylab) +
    coord_fixed(1) #+
    # guides(color = guide_legend(override.aes = list(color = condcolors_AxH, alpha = 1, stroke = 1)),
    #        fill = guide_legend(override.aes = list(fill = condcolors_AxH, shape = 21, alpha = 1, stroke = 0.5)),
    #        shape = guide_legend(override.aes = list(shape = colshapes, alpha = 1, stroke = 0.5)))
  PCoAplot
}

### Function to plot color bar -------------------------------------------------------------------------------------
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }	
}

### Volcano plot for differential gene expression -----------------------------------------------------------------
volcanoplot <- function(res) {
  
  ##Highlight genes that have a padj < 0.05
  res$threshold <- ifelse(res$padj < 0.05 & res$log2FoldChange > 0, "Upregulated", ifelse(res$padj < 0.05 & res$log2FoldChange < 0, "Downregulated", "NA"))
  res$log10padj <- -log10(res$padj)
  dat_genes <- data.frame(cbind(res$log2FoldChange, res$log10padj, res$threshold), stringsAsFactors = FALSE)
  colnames(dat_genes) <- c("log2FoldChange", "log10padj", "threshold")
  row.names(dat_genes) <- res$IDGeneInfo
  #dat_genes <- dat_genes[order(dat_genes$log2FoldChange, decreasing = TRUE),]
  dat_genes$log2FoldChange <- as.numeric(dat_genes$log2FoldChange)
  dat_genes$log10padj <- as.numeric(dat_genes$log10padj)
  dat_genes$threshold <- factor(dat_genes$threshold, levels = c("Upregulated", "Downregulated", "NA"), ordered = TRUE)
  #Create volcanoplot
  gVolcano <- dat_genes %>% 
    ggplot(aes(log2FoldChange, log10padj)) + 
    geom_point(aes(color = threshold), alpha=0.7, size=2) +
    scale_color_manual(values = DEGcolors) +
    scale_x_continuous(limits = c(-6,6), breaks = seq(-10,10,2)) + 
    ylim(c(0, 10)) +
    xlab("log2 fold change") +
    ylab("-log10 p-value") + 
    #geom_text_repel(data = dat_genes_LPS_ctrl[1:15, ], aes(label = rownames(dat_genes_LPS_ctrl[1:15, ])), color = "black", size = 2.5, box.padding = unit(0.35, "lines")) +
    theme_bw() +
    theme(legend.position = "none", 
          plot.title = element_text(size = 12, hjust = 0, vjust = 1),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10))
  print(gVolcano)
}

### Gene expression boxplot functions ------------------------------------------------------------------------------------
ggboxplot <- function(gene) {
  plotTitle <- gene_annotation %>% filter(ID == gene) %>% select(Gene_Info)
  subTitle <- gene
  df <- plotCounts(dds, gene = gene, intgroup = "Treatment", returnData = TRUE)
  df$Treatment <- factor(df$Treatment, levels = c("control", "Heat", "Antibiotics", "Antibiotics.Heat"), ordered = TRUE)
  gbplot <- df %>% ggplot(aes(x = Treatment, y = count, color = Treatment, fill = Treatment)) +
    geom_boxplot(show.legend = FALSE) + 
    geom_point(aes(shape = samples$Colony), size = 2, show.legend = FALSE) +
    scale_color_manual(values = condcolors_AxH) +
    scale_fill_manual(values = condfillcolors_AxH) +
    scale_shape_manual(values = colshapes) +
    scale_y_continuous(trans = log2_trans(),
                       breaks = trans_breaks("log2", function(x) 2^x),
                       labels = trans_format("log2", math_format(2^.x))) +
    ggtitle(plotTitle, subtitle = subTitle)
  print(gbplot)
}


genoboxplot <- function(gene) {
  plotTitle <- gene_annotation %>% filter(ID == gene) %>% select(Gene_Info)
  subTitle <- gene
  df <- plotCounts(dds, gene = gene, intgroup = "Treatment", returnData = TRUE)
  df$Treatment <- factor(df$Treatment, levels = c("control", "Heat", "Antibiotics", "Antibiotics.Heat"), ordered = TRUE)
  gbplot <- df %>% ggplot(aes(x = Treatment, y = count, color = Treatment)) +
    geom_boxplot(aes(fill = Treatment), show.legend = FALSE) + 
    geom_point(aes(shape = samples$Colony), size = 2, show.legend = FALSE) +
    facet_grid(.~samples$Colony) +
    ### add interaction line between genotypes
    stat_summary(fun=mean, geom="path", colour="black", size=0.8, aes(group = samples$Colony)) +
    stat_summary(fun=mean, geom="point", colour="black", size=3, aes(shape = samples$Colony, group = samples$Colony), show.legend = FALSE) +
    ###
    scale_color_manual(values = condcolors_AxH) +
    scale_fill_manual(values = condfillcolors_AxH) +
    scale_shape_manual(values = colshapes) +
    scale_y_continuous(trans = log2_trans(),
                       breaks = trans_breaks("log2", function(x) 2^x),
                       labels = trans_format("log2", math_format(2^.x))) +
    ggtitle(plotTitle, subtitle = subTitle)
  print(gbplot)
  ###
  g_bplot <- ggplot_gtable(ggplot_build(gbplot))
  strip_both <- which(grepl('strip-', g_bplot$layout$name))
  k <- 1
  for (i in strip_both) {
    j <- which(grepl('rect', g_bplot$grobs[[i]]$grobs[[1]]$childrenOrder))
    g_bplot$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colcolors[k]
    k <- k+1
  }
  gbp <- grid.draw(g_bplot)
  print(gbp)
}

### Functions to plot KOG delta ranks heatmap ---------------------------------------------------------------------
#from https://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps
draw_main <- function (text, ...) 
{
  res = textGrob(text, gp = gpar(fontface = "plain", ...), hjust = 0)
  return(res)
}
#overwrite default draw_main with left-justified version 
assignInNamespace(x="draw_main", value="draw_main",
                  ns=asNamespace("pheatmap"))

draw_colnames_330 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 1, hjust = 0, rot = 330, gp = gpar(...))
  return(res)}

assignInNamespace(x="draw_colnames", value="draw_colnames_330",
                  ns=asNamespace("pheatmap"))

KOGheatmap <- function(deltaranks, pvals, ...) {
  deltaranks <- as.matrix(deltaranks)
  paletteLength <- 100
  myColor <- colorRampPalette(rev(c("red", "white", "blue")))(paletteLength)
  myBreaks <- c(seq(min(deltaranks), 0, length.out = ceiling(paletteLength/2) + 1), 
                seq(max(deltaranks)/paletteLength, max(deltaranks), 
                    length.out = floor(paletteLength/2)))
  pheatmap(mat = deltaranks, display_numbers = pvals,
           cluster_cols = FALSE, cluster_rows = FALSE,
           treeheight_row = 15, treeheight_col = 15,
           number_color = "black",
           border_color = "black", scale = "none",
           color = myColor, breaks = myBreaks, ...)
}