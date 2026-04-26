rm(list = ls(all.names = TRUE)) # clears global environ.

# Load packages
library(tidyverse)
library(dplyr) 
library(broom)
library(tibble)
library(anndataR)
library(ComplexHeatmap)
library(circlize)
library(GSVA)
library(GSEABase)
library(BiocParallel)
library(rentrez)
library(org.Hs.eg.db)

set_entrez_key("")

setwd("./Output")

Data <- read_h5ad("combined_emb_normalized.h5ad", mode = 'r')

# Isolate metadata
Metadata <- Data$obs
Var <- Data$var
Uns <- Data$uns

###############################################################

if(!file.exists('./Entrez_summaries.csv')){

  # Load the human db
  Human_genes <- org.Hs.eg.db
  
  # Find the entrez IDs for the symbols in the dataset
  Entrez <- select(Human_genes, 
                   keys = as.character(Var$feature_name),
                   columns = c("ENTREZID", "SYMBOL"),
                   keytype = "SYMBOL")
  
  # Set the symbol column to feature_name to match Var
  colnames(Entrez)[1] <- 'feature_name'
  
  # Save the Entrez ID information to manually add Entrez ID for old gene names which do not match current db
  if(!file.exists('./Entrez_ids.csv')){write.csv(file = 'Entrez_ids.csv', Entrez)}
  
  # Load in manually adjusted Entrez IDs
  Entrez <- read.csv('./Entrez_ids.csv')[-1]
  
  gene_summary <- list()
  
  for(i in 1:50){
    
    start = ((i - 1) * 50) + 1
    stop = 50 * i
    
    # Query NCBI to get gene summary data for each Entrez ID
    Query <- entrez_summary(db = "gene", id = Entrez$ENTREZID[start:stop]) %>% sapply(., "[[", "summary")
    
    # Change into data frame with Entrez ID (which is saved as list object name) as row name
    gene_summary[[i]] <- data.frame(Query, row.names = names(Query))
    
  }
  
  # bind results into single data frame
  gene_summary <- do.call(bind_rows, gene_summary)

  # Set row name (Entrez ID) to a column and rename as ENTREZID
  gene_summary <- gene_summary %>% rownames_to_column()
  colnames(gene_summary)[1] <- 'ENTREZID'
  
  # Merge the summary data
  Entrez <- merge(Entrez, gene_summary, by = 'ENTREZID', all.x = TRUE)
  
  write.csv(file = 'Entrez_summaries.csv', Entrez)
  
}

Entrez <- read.csv('./Entrez_summaries.csv')[-1]

Var <- Var %>% rownames_to_column()
Var <- merge(Var, Entrez, by = 'feature_name', all.x = TRUE)
Var <- data.frame(Var[-2], row.names = Var$rowname)

###############################################################

# Reformat group designation
Metadata$ann_finest_level <- sub(pattern = ' \\(', replacement = '_', Metadata$ann_finest_level)
Metadata$ann_finest_level <- sub(pattern = '\\-', replacement = '_', Metadata$ann_finest_level)
Metadata$ann_finest_level <- sub(pattern = '\\)', replacement = '', Metadata$ann_finest_level)
Metadata$ann_finest_level <- sub(pattern = ' ', replacement = '_', Metadata$ann_finest_level)

category_mapping = c('Basal_resting' = 'Basal', 
                     'Suprabasal' = 'Basal', 
                     'Hillock_like' = 'Basal',
                     
                     'Goblet_subsegmental' = 'Non_nasal_secretory',
                     'Club_non_nasal' = 'Non_nasal_secretory',
                     'Goblet_bronchial' = 'Non_nasal_secretory',
                     
                     'Club_nasal' = 'Nasal_secretory', 
                     'Goblet_nasal' = 'Nasal_secretory', 
                     
                     'Multiciliated_nasal' = 'Ciliated', 
                     'Multiciliated_non_nasal' = 'Ciliated',
                     'Deuterosomal' = 'Ciliated')

Metadata$Group <- category_mapping[Metadata$ann_finest_level]

###############################################################

colors <- list(
  Cell_Type = c('Basal_resting' = '#1f77b4', 
                'Suprabasal' = '#ffbb78', 
                'Hillock_like' = '#bcbd22',
                
                'Goblet_subsegmental' = '#e377c2',
                'Club_non_nasal' = '#2ca02c',
                'Goblet_bronchial' = '#9467bd',
                
                'Club_nasal' = '#ff7f0e', 
                'Goblet_nasal' = '#8c564b', 
                
                'Multiciliated_nasal' = '#17becf', 
                'Multiciliated_non_nasal' = '#aec7e8',
                'Deuterosomal' = '#d62728', 
                
                'Basal' = '#1f77b4', 
                'Non_nasal_secretory' = '#d62728', 
                'Nasal_secretory' = '#2ca02c',
                'Ciliated' = '#ff7f0e'), 
  
  Cluster = c('#000000', '#0057E9', '#E11845', '#D3D3D3', '#F2CA19'))

names(colors[['Cluster']]) <- seq(1, length(colors[['Cluster']]), by = 1)

###############################################################

use_top_genes = TRUE
Abs_max = FALSE
Use_summary = TRUE
do_plot = FALSE
x_genes = 100
Subcluster_n = 5

Sample_groups <- list('All', 'Nasal_secretory', 'Non_nasal_secretory', 'Basal', 'Ciliated')

Heatmaps <- lapply(Sample_groups[1], function(Single_sample_group){
  
  # Isolate samples in grouping
  ifelse(Single_sample_group == 'All',
         current_samples <- Metadata, 
         current_samples <- Metadata[Metadata$Group == Single_sample_group, ])
  
  # Import LR model coefficients and rename groups to match metadata format
  Predictive_genes <- read.csv(file = paste0('./LR Models/LR Model Gene Coefficients/0.5_Predictive_genes_', Single_sample_group, '.csv'), row.names = 1)
  rownames(Predictive_genes) <- sub(pattern = ' \\(', replacement = '_', rownames(Predictive_genes))
  rownames(Predictive_genes) <- sub(pattern = '\\-', replacement = '_', rownames(Predictive_genes))
  rownames(Predictive_genes) <- sub(pattern = '\\)', replacement = '', rownames(Predictive_genes))
  rownames(Predictive_genes) <- sub(pattern = ' ', replacement = '_', rownames(Predictive_genes))
  
  # Calculate the absolute sum of coefficients for each gene
  Top_predictive_genes <- apply(Predictive_genes, 2, function(x) sum(abs(x))) %>% data.frame() %>% setNames(c('Total_coef'))
  
  ###############################################################
  
  # Determine which maximum positive coefficient for each gene and its corresponding group
  Max_coef <- apply(Predictive_genes, 2, function(x) x[which.max(x)])
  Max_coef_names <- apply(Predictive_genes, 2, function(x) rownames(Predictive_genes)[which.max(x)])
  
  ###############################################################
  
  # Create a data frame to store the coefficient data and add the gene type, i.e. protein coding vs lncRNA
  Gene_graphing_data <- data.frame(Top_predictive_genes, Max_coef, Max_coef_names)
  Gene_graphing_data$Type <- Var[rownames(Gene_graphing_data), 'feature_type']
  
  # If the LR coefficients only have one column, add the reference group and the top genes of that group (top negative coefficients)
  if(nrow(Predictive_genes) != length(unique(current_samples$ann_finest_level)) & Single_sample_group != 'All'){
    
    Unused_name <- unique(current_samples$ann_finest_level)[!(unique(current_samples$ann_finest_level) %in% rownames(Predictive_genes))]
    Gene_graphing_data$Max_coef_names[Gene_graphing_data$Max_coef < 0] <- Unused_name
    Gene_graphing_data$Max_coef <- abs(Gene_graphing_data$Max_coef)
    
  }
  
  # Determine the percent the max group coefficient makes up of the absolute coefficient sum
  Gene_graphing_data$Percent_max_of_total <- abs(Gene_graphing_data$Max_coef)/Gene_graphing_data$Total_coef
  
  ###############################################################
  
  # Isolate the top x genes for each group
  Top_x_genes <- lapply(unique(Gene_graphing_data$Max_coef_names), function(Single_group){
    
    Subset <- Gene_graphing_data[which(Gene_graphing_data$Max_coef_names == Single_group), ] %>% 
      arrange(desc(abs(Max_coef))) %>%
      dplyr::slice(1:x_genes)
    
  })
  
  # Arrange result by total coefficient
  Gene_graphing_data <- do.call(bind_rows, Top_x_genes) %>% arrange(desc(Total_coef))
  
  # Convert the max coefficient group into a factor
  Gene_graphing_data$Max_coef_names <- factor(Gene_graphing_data$Max_coef_names)
  
  ###############################################################
  
  # Isolate a slice of the count data corresponding to the selected samples and top genes
  Count_data <- data.frame(Data$X[rownames(current_samples), rownames(Gene_graphing_data)])

  Group_order <- list()
  Gene_order <- list()
  Gene_sets <- list()
  clu <- list()
  Dendrograms <- list()
  Max_clusters = 0

  for(i in 1:length(levels(Gene_graphing_data$Max_coef_names))){

    # Isolate the current group
    Current_group_name <- levels(Gene_graphing_data$Max_coef_names)[[i]]
    Group_order[[i]] <- Current_group_name
    
    # Isolate the top genes for the current group
    Gene_subset <- rownames(Gene_graphing_data[which(Gene_graphing_data$Max_coef_names == Current_group_name), ])
    Gene_order[[i]] <- Gene_subset
    
    # Isolate the samples in the current group
    ifelse(Single_sample_group == 'All',
           Sample_subset <- rownames(current_samples[which(current_samples$Group == Current_group_name), ]),
           Sample_subset <- rownames(current_samples[which(current_samples$ann_finest_level == Current_group_name), ]))
    
    ###############################################################
    
    # Isolate the current count data for the current samples and genes
    Count_data_subset <- Count_data[Sample_subset, Gene_subset]
    
    # Calculate the spearman correlation of gene counts, cluster genes by distance, and build a dendrogram of the results
    Count_cor <- cor(Count_data_subset, method = c("spearman")) %>% data.frame()
    Count_cor_hclust <- hclust(as.dist(1 - Count_cor))
    Count_cor_dend <- as.dendrogram(Count_cor_hclust)
    Dendrograms[[i]] <- Count_cor_dend
    
    ###############################################################
    
    # Cut the current dendrogram into k trees
    # Add the current number of clusters to the results so each group is not 1:k trees and instead 1+k:k+k
    clu[[i]] <- dendextend::cutree(Count_cor_dend, k = Subcluster_n) %>% data.frame() %>% setNames(c('Cluster_id'))
    Genes_in_clusters <- clu[[i]] + Max_clusters
    Genes_in_clusters <- Genes_in_clusters %>% rownames_to_column()
    
    # For each cut tree, make a gene set containing the genes in that tree
    for(i in min(Genes_in_clusters$Cluster_id):max(Genes_in_clusters$Cluster_id)){
      
      Genes_for_gene_set <- Genes_in_clusters[which(Genes_in_clusters$Cluster_id == i), 'rowname']
      Gene_sets[[i]] <- GeneSet(Genes_for_gene_set, setName = paste(Current_group_name, i - Max_clusters, sep = '_'))
      
    }
    
    # Add the number of cut trees to the cluster tracker
    Max_clusters = Max_clusters + Subcluster_n

  }
  
  ###############################################################
  
  Merged_dend <- Dendrograms[[1]]

  for(i in 2:length(Dendrograms)){

    Merged_dend <- merge(Merged_dend, Dendrograms[[i]])
    
  }
  
  ###############################################################
  
  Gene_order <- unlist(Gene_order)
  
  Clustered_genes <- do.call(bind_rows, clu)
  colnames(Clustered_genes) <- c('Cluster')

  Gene_graphing_data <- merge(Clustered_genes, Gene_graphing_data, by = 0)
  Gene_graphing_data <- data.frame(Gene_graphing_data[-1], row.names = Gene_graphing_data$Row.names)
  
  Gene_graphing_data$Gene_symbol <- Var[rownames(Gene_graphing_data), 'feature_name']
  Gene_graphing_data$Gene_summary <- Var[rownames(Gene_graphing_data), 'Query']
  Gene_graphing_data <- Gene_graphing_data[match(unlist(Gene_order), rownames(Gene_graphing_data)), ]
  
  # Reorder count data to match gene clustering order
  Unscaled_count_data <- Count_data[, rownames(Gene_graphing_data)]
  
  boxplot_limits <- round(max(Unscaled_count_data), 0)
  
  ###############################################################
  
  Group_colors <- colors[['Cell_Type']][names(colors[['Cell_Type']]) %in% levels(Gene_graphing_data$Max_coef_names)]
  Cluster_colors <- colors[['Cluster']][names(colors[['Cluster']]) %in% unique(Gene_graphing_data$Cluster)]
  
  if(Use_summary == TRUE){
  
    Count_data <- data.frame(scale(Unscaled_count_data))
    
    if(Single_sample_group == 'All'){
      
      Plotting_data <- data.frame(current_samples$Group, Count_data) %>%
        setNames(c('Group', colnames(Count_data))) %>%
        group_by(Group) %>%
        summarize(across(everything(), median))
      
      Plotting_data <- data.frame(Plotting_data[-1], row.names = Plotting_data$Group)
      
    }
    
    ###############################################################
    
    else{
  
      Plotting_data <- data.frame(current_samples$ann_finest_level, Count_data) %>%
        setNames(c('ann_finest_level', colnames(Count_data))) %>%
        group_by(ann_finest_level) %>%
        summarize(across(everything(), median))
      
      Plotting_data <- data.frame(Plotting_data[-1], row.names = Plotting_data$ann_finest_level)
      
    }
    
    col_split <- levels(Gene_graphing_data$Max_coef_names)
    col_order <- levels(Gene_graphing_data$Max_coef_names)
    
    colAnn <- HeatmapAnnotation(Cell_Type = levels(Gene_graphing_data$Max_coef_names),
                                col = list(Cell_Type = Group_colors), 
                                show_annotation_name = FALSE, 
                                show_legend = FALSE)
    
    data_name = 'Median Log Count Z-Score'
    Greatest_value <- max(abs(round(min(Plotting_data), 0)), abs(round(max(Plotting_data), 0)))

  }
  
  ###############################################################
  
  else{
    
    Plotting_data <- data.frame(scale(Unscaled_count_data))
    
    ifelse(Single_sample_group == 'All',
           col_split <- current_samples$Group,
           col_split <- current_samples$ann_finest_level)
    
    col_order <- NULL
    
    colAnn <- HeatmapAnnotation(Cell_Type = col_split,
                                col = list(Cell_Type = Group_colors), 
                                show_annotation_name = FALSE, 
                                show_legend = FALSE)
    
    data_name = 'Log Count Z-Score'
    Greatest_value <- 2

  }
  
  ###############################################################
  
  row_anno_plot <- rowAnnotation("Log Counts" = anno_boxplot(x = Unscaled_count_data,
                                                          which = "row",
                                                          axis = TRUE, 
                                                          size = unit(0, "mm"), 
                                                          ylim = c(0, boxplot_limits)),
                                 
                                 width = unit(3, "cm"), 
                                 show_legend = FALSE)
  
  ###############################################################

  rowAnn <- rowAnnotation(Cluster = Gene_graphing_data$Cluster,
                          Maximum_Coefficient = Gene_graphing_data$Max_coef,
                          Percent_max_of_total = Gene_graphing_data$Percent_max_of_total * 100,
                          Cell_Type = Gene_graphing_data$Max_coef_names,
                          
                          col = list(Cluster = Cluster_colors, 
                                     Maximum_Coefficient = colorRamp2(c(0, 0.5, 1, 1.5, 2), hcl_palette = "spectral"),
                                     Percent_max_of_total = colorRamp2(c(0, 50, 100), hcl_palette = "Blues"),
                                     Cell_Type = Group_colors), 
                          
                          simple_anno_size = unit(0.75, "cm"), 
                          show_annotation_name = FALSE, 
                          show_legend = FALSE)
  
  #####
  
  lgd1 = Legend(at = 1:length(levels(Gene_graphing_data$Max_coef_names)),
                labels = names(Group_colors),
                legend_gp = gpar(fill = Group_colors), 
                title = "Cell Type")
  
  # if(Abs_max == TRUE){lgd2 = Legend(col_fun = colorRamp2(c(-1, 0, 1), hcl_palette = "Blue-Red3"), title = "Max Coefficient")}
  
  lgd2 = Legend(col_fun = colorRamp2(c(0, 0.5, 1, 1.5, 2), hcl_palette = "spectral", rev = FALSE), title = "Max Coefficient") 
  lgd3 = Legend(col_fun = colorRamp2(c(0, 50, 100), hcl_palette = "Blues", rev = TRUE), title = "Percent of Total Coefficient")
  
  ###############################################################
  
  legends = packLegend(list = list(lgd1, lgd2, lgd3))
  
  # myCol <- colorRamp2(c(-Greatest_value, -Greatest_value/2, 0, Greatest_value/2, Greatest_value), hcl_palette = "Blue-Red 2")
  myCol <- colorRamp2(c(-2, -1, 0, 1, 2), hcl_palette = "Blue-Red 2")
  
  ###############################################################

  if(do_plot == TRUE){
  
    heatmap <- Heatmap(as.matrix(t(Plotting_data)),
                       
                       # split the genes / rows
                       column_split = col_split,
                       column_order = col_order,
                       row_split = length(levels(Gene_graphing_data$Max_coef_names)),
                       cluster_column_slices = FALSE,
                       
                       row_gap = unit(2.5, "mm"),
                       column_gap = unit(2.5, 'mm'),
                       border = TRUE,
                       width =  unit(350, 'mm'),
  
                       name = data_name,
                       col = myCol,
                       
                       # parameters for the color-bar that represents gradient of expression
                       heatmap_legend_param = list(
                         color_bar = 'continuous',
                         legend_direction = 'horizontal',
                         legend_width = unit(25, 'cm'),
                         legend_height = unit(25.0, 'cm'),
                         title_position = 'topcenter',
                         title_gp = gpar(fontsize = 24, fontface = 'bold'),
                         labels_gp = gpar(fontsize = 24, fontface = 'bold'), 
                         at = c(-2, -1, 0, 1, 2), 
                         labels = c(-2, -1, 0, 1, 2)),
                       
                       # row (gene) parameters
                       row_title = NULL,
                       row_title_side = 'right',
                       row_title_gp = gpar(fontsize = 24,  fontface = 'bold'),
                       row_title_rot = 0,
                       show_row_names = FALSE,
                       cluster_rows = Merged_dend,
                       show_row_dend = TRUE,
                       row_dend_width = unit(25,'mm'),
  
                       # column (sample) parameters
                       cluster_columns = FALSE,
                       show_column_dend = FALSE,
                       column_title = NULL,
                       show_column_names = FALSE,
                       column_dend_height = unit(25,'mm'),
                       
                       # specify top and bottom annotations
                       top_annotation = colAnn,
                       left_annotation = rowAnn, 
                       right_annotation = row_anno_plot)
    
    png(file = paste0('./Heatmaps/Gene_Heatmaps/', Single_sample_group, " Heatmap.png"), height = 1500, width = 2000)
    draw(heatmap, gap = unit(20, "mm"), heatmap_legend_side = 'bottom')
    dev.off()
    
    ######
    
    png(file = paste0('./Heatmaps/Gene_Heatmaps/', Single_sample_group, " Heatmap Legend.png"), height = 250, width = 350)
    draw(legends)
    dev.off()
    
  }
  
  ###############################################################
  ###############################################################
  ###############################################################
  
  # Create an annotated data frame with the current samples
  phenoData <- new("AnnotatedDataFrame", data = current_samples)
  
  # Create ExpressionSet object containing all count data
  gsva_counts <- ExpressionSet(assayData = as.matrix(t(Data$X[rownames(current_samples), ])), phenoData = phenoData)
  
  # Perform GSVA on gene sets with minsize of 2
  gsva_results <- gsva(gsvaParam(gsva_counts, GeneSetCollection(Gene_sets), minSize = 2),
                       verbose = FALSE,
                       BPPARAM = SerialParam(progressbar = TRUE))
  
  # Isolate GSVA results
  GSVA_metadata <- gsva_results@phenoData@data
  GSVA_counts <- data.frame(t(gsva_results@assayData[["exprs"]]))
  
  ifelse(Single_sample_group == 'All',
         col_split <- GSVA_metadata$Group,
         col_split <- GSVA_metadata$ann_finest_level)
  
  ###############################################################
  
  Group <- list()
  Dendrograms <- list()
  
  for(i in 1:length(levels(Gene_graphing_data$Max_coef_names))){
    
    # Isolate the current group
    Current_group_name <- levels(Gene_graphing_data$Max_coef_names)[[i]]
    Group[[i]] <- Current_group_name
    
    # Isolate the samples in the current group
    ifelse(Single_sample_group == 'All',
           Sample_subset <- rownames(GSVA_metadata[which(GSVA_metadata$Group == Current_group_name), ]),
           Sample_subset <- rownames(GSVA_metadata[which(GSVA_metadata$ann_finest_level == Current_group_name), ]))
    
    #################################
    
    # Isolate the current count data for the current samples and genes
    gsva_results_subset <- GSVA_counts[Sample_subset, grepl(Current_group_name, colnames(GSVA_counts), perl = TRUE)]
    
    # Calculate the spearman correlation of gene counts, cluster genes by distance, and build a dendrogram of the results
    Count_cor <- cor(gsva_results_subset, method = c("spearman")) %>% data.frame()
    Count_cor_hclust <- hclust(as.dist(1 - Count_cor))
    Dendrograms[[i]] <- as.dendrogram(Count_cor_hclust)
    
  }
  
  Merged_dend <- Dendrograms[[1]]
  
  for(i in 2:length(Dendrograms)){
    
    Merged_dend <- merge(Merged_dend, Dendrograms[[i]])
    
  }
  
  ###############################################################
  
  GSVA_counts$Group_var <- col_split
  
  GSVA_25 <- GSVA_counts %>% group_by(Group_var) %>% summarize(across(everything(), ~quantile(., 0.25)))
  GSVA_75 <- GSVA_counts %>% group_by(Group_var) %>% summarize(across(everything(), ~quantile(., 0.75)))
  
  GSVA_counts <- GSVA_counts %>% group_by(Group_var) %>% summarize(across(everything(), median))
  GSVA_counts <- data.frame(GSVA_counts[-1], row.names = GSVA_counts$Group_var)
  
  GSVA_counts <- t(GSVA_counts) %>% data.frame()
  plot_cat <- 'Median GSVA Score'
  
  # Create a group designation for GSVA groupings
  count_groups <- sapply(unique(as.character(Gene_graphing_data$Max_coef_names)), function(x) sum(str_count(rownames(GSVA_counts), x)))
  Groups <- rep(levels(Gene_graphing_data$Max_coef_names), count_groups)
  cluster_number <- as.numeric(unlist(str_extract_all(rownames(GSVA_counts), "\\d+")))
  
  Group_colors <- colors[['Cell_Type']][names(colors[['Cell_Type']]) %in% levels(Gene_graphing_data$Max_coef_names)]
  Cluster_colors <- colors[['Cluster']][names(colors[['Cluster']]) %in% unique(cluster_number)]
  
  colAnn <- HeatmapAnnotation(Cell_Type = unlist(Group),
                              col = list(Cell_Type = Group_colors), 
                              show_annotation_name = FALSE, 
                              show_legend = FALSE)
  
  ###############################################################
  
  rowAnn <- rowAnnotation(Cluster = cluster_number,
                          Cell_Type = Groups,
                          
                          col = list(Cluster = Cluster_colors,
                                     Cell_Type = Group_colors), 
                          
                          simple_anno_size = unit(0.75, "cm"), 
                          show_annotation_name = FALSE, 
                          show_legend = FALSE)
  
  lgd1 = Legend(at = 1:length(levels(Gene_graphing_data$Max_coef_names)),
                labels = names(Group_colors),
                legend_gp = gpar(fill = Group_colors), 
                title = "Cell Type")
  
  legends = packLegend(list = list(lgd1))
  
  ###############################################################
  
  myCol <- colorRamp2(c(-1, -0.5, 0, 0.5, 1), hcl_palette = "Blue-Red 2")
  
  #####
  
  if(do_plot == TRUE){
  
    heatmap <- Heatmap(as.matrix(GSVA_counts),
                       
                       # split the genes / rows
                       column_split = unlist(Group),
                       column_order = unlist(Group),
                       row_split = length(Dendrograms),
                       cluster_column_slices = FALSE,
                       cluster_row_slices = FALSE,
                       
                       row_gap = unit(2.5, "mm"),
                       column_gap = unit(2.5, 'mm'),
                       border = TRUE,
                       width =  unit(350, 'mm'),
                       
                       name = plot_cat,
                       col = myCol,
                       
                       # Add median GSVA values
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(paste0(sprintf("%.4f", GSVA_counts[i, j]),
                                          '\n',
                                          '(',
                                          sprintf("%.4f", GSVA_25[-1][j, i]),
                                          '   -  ',
                                          sprintf("%.4f", GSVA_75[-1][j, i]),
                                          ')'), 
                                   
                                   x, y, gp = gpar(fontsize = 20, fontface = 'bold'))},
                       
                       # parameters for the color-bar that represents gradient of expression
                       heatmap_legend_param = list(
                         color_bar = 'continuous',
                         legend_direction = 'horizontal',
                         legend_width = unit(25, 'cm'),
                         legend_height = unit(25.0, 'cm'),
                         title_position = 'topcenter',
                         title_gp = gpar(fontsize = 24, fontface = 'bold'),
                         labels_gp = gpar(fontsize = 24, fontface = 'bold'), 
                         at = c(-1, -0.5, 0, 0.5, 1), 
                         labels = c(-1, -0.5, 0, 0.5, 1)),
                       
                       # row (gene) parameters
                       row_title = NULL,
                       row_title_rot = 0,
                       show_row_names = FALSE,
                       row_names_side = 'left',
                       row_names_gp = gpar(fontsize = 24, fontface = 'bold'),
                       cluster_rows = Merged_dend,
                       show_row_dend = TRUE,
                       row_dend_width = unit(25,'mm'),
                       
                       # column (sample) parameters
                       cluster_columns = FALSE,
                       show_column_dend = FALSE,
                       column_title = NULL,
                       show_column_names = FALSE,
                       column_dend_height = unit(25,'mm'),
                       
                       # specify top and bottom annotations
                       top_annotation = colAnn,
                       left_annotation = rowAnn)
    
    png(file = paste0('./Heatmaps/GSVA_Heatmaps/', Single_sample_group, " GSVA Heatmap.png"), height = 1500, width = 2000)
    
    draw(heatmap, gap = unit(20, "mm"), heatmap_legend_side = 'bottom')
    draw(legends, x = unit(58, "cm"), y = unit(25, "cm"))
    
    dev.off()
    
  }
  
  var_subset <- Var[rownames(Gene_graphing_data), c("feature_type", "ENTREZID")]
  Gene_graphing_data <- cbind(Gene_graphing_data, var_subset)
  
  write.csv(file = paste0('./Heatmaps/Gene_Data/', Single_sample_group,'_Cluster_Results.csv'), Gene_graphing_data)
  
})

names(Heatmaps) <- Sample_groups[1]

save(Heatmaps, Sample_groups, colors, Var, Metadata, file = 'AUC_Data.RData')
