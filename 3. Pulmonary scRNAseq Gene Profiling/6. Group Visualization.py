import math
import alphashape
import numpy as np
import scipy.sparse
import scanpy as sp
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from shapely.geometry import Polygon, MultiPolygon

def determine_distance(coordinate_1, coordinate_2):
    '''
    Given the coordinates of two points, determine the distance between them.

    Parameters:
    -----------
    coordinate_1, coordinate_2: Pair of coordinates for calculation of distances.
    '''    
    x0, y0 = coordinate_1
    x1, y1 = coordinate_2
    distance = math.sqrt((x1 - x0)**2 + (y1 - y0)**2)

    return distance

###############

def find_perimeter(anndata, category, alpha = 0.95, dim_red = 'X_umap', plot = False):
    '''
    Find the perimeter of a series of related points using alphashape hulls.

    Parameters:
    -----------
    anndata: Input AnnData object containing counts for hull calculations.
    category: Name of column in anndata.obs for grouping for hull calculations.
    alpha: Alpha value for alphashape calculation. Defaul to 0.95.
    dim_red: Name of data in anndata.obsm containing data from dimensionality reduction. Default to X_umap
    plot: Whether to plot the results. Default to False.
    '''

    if isinstance(anndata, sp.AnnData) != True:
        raise TypeError('Input must be an AnnData object')
    
    if category not in anndata.obs.columns:
        raise TypeError('category must be the name of a column in anndata.obs')
    
    ##############

    hull_group = []
    hull_points = []

    # Isolate data for each category and calculate alphashape hulls
    for cluster in np.unique(anndata.obs[category]):
        hull_group.append(cluster)
        index = np.where(anndata.obs[category] == cluster)
        data_subset = anndata.obsm[dim_red][index]
        center_point = [data_subset[:, 0].mean(), data_subset[:, 1].mean()]

        distance = []
        for point in data_subset:
            distance.append(determine_distance(center_point, point))

        sort_distance = np.argsort(distance)
        ordered_data_subset = data_subset[sort_distance][0:math.ceil(len(distance)*0.95)]
        
        hull = alphashape.alphashape(ordered_data_subset, alpha)

        if isinstance(hull, MultiPolygon):
            largest_shape = max(hull.geoms, key = lambda p: p.area)

        elif isinstance(hull, Polygon):
            largest_shape = hull

        hull_pts = largest_shape.exterior.coords.xy
        coords = np.column_stack((hull_pts[0].tolist(), hull_pts[1].tolist()))
        hull_points.append(coords)

    if plot == True:

        # Get min and max values for x and y axis limits
        min_x = anndata.obsm[dim_red][:, 0].min()
        max_x = anndata.obsm[dim_red][:, 0].max()
        min_y = anndata.obsm[dim_red][:, 1].min()
        max_y = anndata.obsm[dim_red][:, 1].max()

        # Plot data points and related hulls calculated above
        fig, ax = plt.subplots()
        ax.scatter(data_subset[:, 0], data_subset[:, 1])
        ax.add_patch(patches.Polygon(coords, facecolor = 'none', edgecolor = 'red', linewidth = 1))
        plt.axis([min_x, max_x, min_y, max_y])
        plt.title(cluster)

    return hull_group, hull_points 

###############

def group_expression_metrics(anndata, count_layer_name, grouping_name):
    '''
    Function to calculate expression metrics for specified sample groups.

    Parameters:
    -----------   
    anndata: AnnData object containing scRNAseq count data and metadata.
    count_layer_name: Name of the AnnData layer containing count data.
    grouping_name: Name of the column in the AnnData.obs layer containing grouping information.
    '''
    if isinstance(anndata, sp.AnnData) != True:
        raise TypeError('Input must be an AnnData object')
    
    ##############
    
    if count_layer_name == 'X':
        counts = anndata.X

    elif count_layer_name in anndata.layers:
        counts = anndata.layers[count_layer_name]

    else:
        raise TypeError('count_layer_name must be X or the name of a layer in anndata.')
    
    if isinstance(counts, scipy.sparse.spmatrix):
        counts = counts.toarray()
    
    ##############
    
    if grouping_name in anndata.obs.columns:
        metadata = anndata.obs[grouping_name]

    else:
        raise TypeError('grouping_name must be the name of a column in anndata.obs')
    
    ##############

    percent_counts = [] # Percent of cells in group which have any counts for a given gene
    group_means = [] # Mean counts per group for a given gene
    group_medians = [] # Median counts per group for a given gene

    for group in np.unique(metadata):
        group_counts = counts[np.where(metadata == group)]
        percent_counts.append((np.sum(group_counts != 0, axis = 0)/group_counts.shape[0])*100)
        group_means.append(np.mean(group_counts, axis = 0))
        group_medians.append(np.median(group_counts, axis = 0))

    percent_counts = pd.DataFrame(percent_counts, index = np.unique(metadata))
    means = pd.DataFrame(group_means, index = np.unique(metadata))
    medians = pd.DataFrame(group_medians, index = np.unique(metadata))

    return percent_counts, means, medians

###################################################################

directory = './Output'

combined_emb = sp.read_h5ad(f'{directory}/combined_emb_normalized.h5ad')

# Create mask for data with 0 counts
mask = (combined_emb.layers['counts'] != 0).toarray()
mask = np.logical_not(mask)

##################################################################

# Isolate groups and add an 'All' category
Groups = np.unique(combined_emb.obs.Group).tolist()
Groups.insert(0, 'All')

for sample_group in Groups:

    # Isolate data for current group
    if sample_group == 'All':
        grouping_cat = 'Group'
        Sample_subset = combined_emb
        Sample_index = list(range(len(combined_emb.obs)))
        
    else:
        grouping_cat = 'ann_finest_level'
        Sample_subset = combined_emb[combined_emb.obs.Group == sample_group, :]
        Sample_index = np.where(combined_emb.obs.Group == sample_group)[0].tolist()

    # Calculate percent counts and subset genes with less than 25% of samples demonstrating counts
    percent_counts, _, _ = group_expression_metrics(Sample_subset, 'counts', grouping_cat)
    Sample_subset = Sample_subset[:, percent_counts.max(axis = 0) > 25].copy()

    # Rank genes for group
    sp.tl.rank_genes_groups(Sample_subset, grouping_cat, method = "wilcoxon")

    # Get gene ranks with p-value and log fold change cutoffs
    Gene_ranks = sp.get.rank_genes_groups_df(Sample_subset, group = None, gene_symbols = 'feature_name', pval_cutoff = 0.05)
    Gene_ranks = Gene_ranks.loc[np.abs(Gene_ranks.logfoldchanges) > 1, :]

    ####################################

    linked_dict_symbols = {}
    top_genes = []
    top_gene_symbols = []
    
    # Isolate top 10 genes for each subgroup
    for group in np.unique(Gene_ranks.group):
        Gene_ranks_subset = Gene_ranks.loc[Gene_ranks.group == group, :]
        linked_dict_symbols[group] = Gene_ranks_subset.feature_name[0:10].tolist()
        top_genes.append(Gene_ranks_subset.names[0:10].tolist())
        top_gene_symbols.append(Gene_ranks_subset.feature_name[0:10].tolist())

    top_genes = [items for sublist in top_genes for items in sublist]
    top_gene_symbols = [items for sublist in top_gene_symbols for items in sublist]

    ####################################
    
    # Isolate unique genes
    seen_gene = set()
    duplicated_gene = {}

    for gene in top_genes:
        if gene in seen_gene:
            gene_symbol = combined_emb.var.loc[gene, 'feature_name']
            duplicated_gene[gene_symbol] = gene
        else:
            seen_gene.add(gene)

    ####################################

    Gene_index = np.where(combined_emb.var.index.isin(np.unique(top_genes)))[0].tolist() # Save the index for each gene to plot
    mask_subset = mask[np.ix_(Sample_index, Gene_index)] # Subset the mask to match the genes to plot
    Top_gene_subset = combined_emb[Sample_index, Gene_index].copy()
    Top_gene_subset.X[mask_subset] = 0 # Set masked values to 0

    ####################################
    
    Top_gene_subset.var.index = Top_gene_subset.var.feature_name

    # Create ordered dot plot for top gene exptession
    sp.pl.dotplot(Top_gene_subset, 
                  linked_dict_symbols, 
                  gene_symbols = 'feature_name',
                  groupby = grouping_cat, 
                  show = False, 
                  dendrogram = False, 
                  categories_order = linked_dict_symbols.keys())

    # If sample_group is 'All', set below grouping and ordering
    if sample_group == 'All':

        sp.pl.dotplot(Top_gene_subset, 
                      linked_dict_symbols, 
                      gene_symbols = 'feature_name',
                      groupby = 'ann_finest_level', 
                      show = False, 
                      dendrogram = False, 
                      categories_order = ['Hillock-like', 'Basal resting', 'Suprabasal',
                                          'Deuterosomal', 'Multiciliated (nasal)', 'Multiciliated (non-nasal)',
                                          'Club (nasal)', 'Goblet (nasal)', 
                                          'Goblet (subsegmental)', 'Club (non-nasal)', 'Goblet (bronchial)'])

    ####################################

    Sample_subset.var.index = Sample_subset.var.feature_name

    sp.pp.pca(Sample_subset) # Calculate PCA
    sp.pp.neighbors(Sample_subset) # Calculate neighbors
    sp.tl.umap(Sample_subset) # Calculate UMAP

    # Calculate alphashape hull for current samples and plot
    hull_group, hull_points = find_perimeter(Sample_subset, grouping_cat, alpha = 0.95)

    loc = ['on data', 'right margin', 'none']

    fig, axes = plt.subplots(1, 3, figsize = (15, 5))

    for i, plot_var in enumerate([grouping_cat, 'anatomical_region_ccf_score', 'ref_or_query']):

        sp.pl.umap(Sample_subset, color = plot_var, ax = axes[i], show = False, legend_loc = loc[i])

        edge_color = Sample_subset.uns[grouping_cat + '_colors']

        for j, hull in enumerate(hull_points):
                            
            axes[i].add_patch(patches.Polygon(hull,
                                              alpha = 0.75,
                                              facecolor = 'none',
                                              edgecolor = edge_color[j],
                                              linewidth = 2))

    ####################################

    # For each top gene, plot hulls with coloring based on gene expression
    for key in linked_dict_symbols.keys():

        genes = linked_dict_symbols[key]

        fig, axes = plt.subplots(1, 10, figsize = (50, 5))

        for i, gene in enumerate(genes):

            sp.pl.umap(Sample_subset, color = gene, ax = axes[i], show = False, legend_loc = 'right margin')
            
            edge_color = Sample_subset.uns[grouping_cat + '_colors']

            for j, hull in enumerate(hull_points):
                    
                axes[i].add_patch(patches.Polygon(hull,
                                            alpha = 0.75,
                                            facecolor = 'none',
                                            edgecolor = edge_color[j],
                                            linewidth = 2))

####################################

sp.pp.pca(combined_emb) # Calculate PCA
sp.pp.neighbors(combined_emb) # Calculate neighbors
sp.tl.umap(combined_emb) # Calculate UMAP

combined_emb.var.index = combined_emb.var.feature_name
hull_group, hull_points = find_perimeter(combined_emb, 'Group', alpha = 0.95)

search_group = '' # String corresponding to gene name/family to find. For example, 'CYP' to find Cytochrome P450 genes
genes = combined_emb.var.loc[combined_emb.var['feature_name'].str.contains(search_group, case = False), 'feature_name'].tolist()

# Groups of related genes for plotting
# For example, group1 could be all CYP genes, group2 all mucins, group3 all tubulin genes, etc.
group1 = []
group2 = []
group3 = []

Genes = [group1, group2, group3]

# For each group of genes, plot the UMAP results, group hulls, and color by gene
# Create panels based on number of input genes
for sublist_genes in Genes:
        
        sorted_sublist = sorted(sublist_genes, key = lambda x: int("".join([i for i in x if i.isdigit()])))

        for loop_n in range(math.ceil(len(sorted_sublist) / 5)):

            start = loop_n * 5

            if start + 5 <= len(sorted_sublist):
                stop = start + 5
                fig, axes = plt.subplots(1, 5, figsize = (25, 5))

            else:
                stop = len(sorted_sublist)
                n_panels = len(sorted_sublist) % 5
                fig, axes = plt.subplots(1, n_panels, figsize = (n_panels * 5, 5))

            for i, gene in enumerate(sorted_sublist[start:stop]):

                sp.pl.umap(combined_emb, color = gene, ax = axes[i], show = False, legend_loc = 'right margin')
                
                face_colors = combined_emb.uns['Group_colors']
                edge_color = combined_emb.uns['Group_colors']

                for j, hull in enumerate(hull_points):
                        
                    axes[i].add_patch(patches.Polygon(hull,
                                                alpha = 0.75,
                                                facecolor = 'none',
                                                edgecolor = edge_color[j],
                                                linewidth = 2))
