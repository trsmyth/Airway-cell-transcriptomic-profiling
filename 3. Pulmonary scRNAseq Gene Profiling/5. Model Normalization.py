import scvi
import os.path
import numpy as np
import scipy.sparse
import scanpy as sp
import pandas as pd

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

###############

directory = './Output'

surgery_model = scvi.model.SCANVI.load(f"{directory}/scvi models/surgery_model.pt")
combined_emb = sp.read_h5ad(f'{directory}/combined_emb_raw.h5ad')

types_of_interest = ['Basal resting', 'Suprabasal', 'Hillock-like',
                     'Club (nasal)', 'Club (non-nasal)',
                     'Goblet (bronchial)', 'Goblet (nasal)', 'Goblet (subsegmental)',
                     'Multiciliated (nasal)', 'Multiciliated (non-nasal)', 'Deuterosomal']

# Isolate groups of interest
combined_emb = combined_emb[combined_emb.obs.ann_finest_level.isin(types_of_interest)].copy()

##################################################################

# Create 'Group' category based on known relation between cell types. Some of this was 
# determined in seperate runs and introduced here for ease of use moving forward.
category_mapping = {'Basal resting' : 'Basal', 
                    'Suprabasal' : 'Basal',
                    'Hillock-like' : 'Basal', 

                    'Goblet (subsegmental)' : 'Non_nasal_secretory',
                    'Club (non-nasal)' : 'Non_nasal_secretory',
                    'Goblet (bronchial)' : 'Non_nasal_secretory',

                    'Club (nasal)' : 'Nasal_secretory', 
                    'Goblet (nasal)' : 'Nasal_secretory', 

                    'Multiciliated (nasal)' : 'Ciliated', 
                    'Multiciliated (non-nasal)' : 'Ciliated',
                    'Deuterosomal' : 'Ciliated'}

combined_emb.obs['Group'] = combined_emb.obs['ann_finest_level'].map(category_mapping)
combined_emb.layers['counts'] = combined_emb.X.copy()
del combined_emb.X # Delete the X layer to free memory

###############

unnormalized_percent_counts, unnormalized_means, unnormalized_medians = group_expression_metrics(combined_emb, 'counts', 'ann_finest_level')

# Create a mask for locations that are zero
mask = (combined_emb.layers['counts'] != 0).toarray()
mask = np.logical_not(mask)

##################################################################
'''
At this point, combined_emb.X are raw counts.
As the scvi/scANVI/scARCHS model uses the 'counts' layer, the raw
count stored in X are copied and saved in this layer for normalization.

get_normalized_expression() takes the raw count data stored in 'counts' and
performes library scaling (1e4 library size), de-noising, and batch correction
by projecting count data onto the batch of choice (or average of multiple batches
if provided as a list of batches). This provides count data which is scaled to 
the specified library size and batch corrected **BUT NOT LOG TRANSFORMED**.
'''
##################################################################

# Determine the number of cells of each group per experiment and save as a csv to determine which
# experiments serve as good transform_batch below
dataset_counts = []

for dataset in np.unique(combined_emb.obs.dataset):
    dataset_counts.append(combined_emb.obs[combined_emb.obs.dataset == dataset].groupby('ann_finest_level', observed = True).count().iloc[:, 0])

dataset_counts = pd.DataFrame(dataset_counts, index = np.unique(combined_emb.obs.dataset))
dataset_counts['row_total'] = dataset_counts.sum(axis = 1)
dataset_counts = dataset_counts.sort_values(by = ['row_total'], ascending = False)

if os.path.isfile(f'{directory}/Dataset_cell_counts.csv') == False:
    dataset_counts.to_csv(f'{directory}/Dataset_cell_counts.csv')

##################################################################

# Get normalized count data using mean of three listed batches .These three batches together provide 
# strong coverage of each cell type while maintaining overall trends. Use batch size of 10k for 
# computational ease and scale to library size of 1e4.
Normalized_counts = surgery_model.get_normalized_expression(combined_emb, 
                                                            library_size = 1e4, 
                                                            transform_batch = ['Barbry_Leroy_2020', 'Nawijn_2021', 'Jain_Misharin_2021_10Xv2'], 
                                                            batch_size = 10000)

combined_emb.X = Normalized_counts # Set X to normalized counts
combined_emb = combined_emb[:, unnormalized_percent_counts.max(axis = 0) > 25].copy() # Isolate genes with at least 25% expression in a single group
mask = mask[:, unnormalized_percent_counts.max(axis = 0) > 25] # Subset the mask to match the isolated genes
combined_emb.layers['norm_counts'] = combined_emb.X # Set X layer to normalized counts
del Normalized_counts # Delete normalized counts to free memory

##################################################################

sp.pp.log1p(combined_emb)

# Analyze top expressed genes per group over all groups (ann_finest_level) or
# more focused 'Group' category which collapses related cell types together
for sample_group in ["ann_finest_level", "Group"]:

    # Rank genes by current group to find marker genes
    sp.tl.rank_genes_groups(combined_emb, sample_group, method = "wilcoxon")
    sp.pl.rank_genes_groups_tracksplot(combined_emb, n_genes = 10, gene_symbols = 'feature_name')
    combined_emb.uns[sample_group + '_rank_genes_groups'] = combined_emb.uns['rank_genes_groups']

    ####################################

    # Isolate top 10 genes per group
    linked_dict = {field: combined_emb.uns['rank_genes_groups']['names'][field][:10].tolist() for field in combined_emb.uns['rank_genes_groups']['names'].dtype.names}
    top_genes = [items for sublist in linked_dict.values() for items in sublist]

    ####################################
    
    # Find uniquely highly expressed genes
    seen_gene = set()
    duplicated_gene = {}

    for gene in top_genes:
        if gene in seen_gene:
            gene_symbol = combined_emb.var.loc[gene, 'feature_name']
            duplicated_gene[gene_symbol] = [key for key, value in linked_dict.items() if gene in value]
        else:
            seen_gene.add(gene)

    ####################################

    Top_gene_subset = combined_emb[:, np.unique(top_genes)].copy()
    Gene_index = np.where(combined_emb.var.index.isin(np.unique(top_genes)))[0].tolist() # Save the index for each gene to plot
    mask_subset = mask[:, Gene_index] # Subset the mask to match the genes to plot

    ####################################

    linked_dict_symbols = {}

    for group, gene_names in zip(linked_dict.keys(), linked_dict.values()):
        if group not in linked_dict_symbols: # If the current group is not already in the dict
            linked_dict_symbols[group] = [] # Create a storage list for the values to be linked to that group
        symbols = Top_gene_subset.var.loc[gene_names, 'feature_name'].tolist()

        for gene in symbols:
            linked_dict_symbols[group].append(gene) # Add the gene names to the storage list

    Gene_results = pd.DataFrame(linked_dict_symbols)
    Gene_results.to_csv(f'{sample_group}_top_genes.csv')
    
    group_order = ['Deuterosomal', 'Multiciliated (nasal)', 'Multiciliated (non-nasal)',
                   'Hillock-like', 'Basal resting', 'Suprabasal',
                   'Club (nasal)', 'Goblet (nasal)', 
                   'Goblet (subsegmental)', 'Club (non-nasal)', 'Goblet (bronchial)']
    
    if sample_group == 'Group':
        row_order = ['Ciliated', 'Basal', 'Nasal_secretory', 'Non_nasal_secretory']
        linked_dict_symbols = {key: linked_dict_symbols[key] for key in row_order}

    else:
        row_order = group_order
        linked_dict_symbols = {key: linked_dict_symbols[key] for key in group_order}

    Top_gene_subset.var.index = Top_gene_subset.var.feature_name

    ####################################

    # Create dot plots with and without masks for expression
    for mask_true in [False, True]:

        if mask_true == True:

            Top_gene_subset.X[mask_subset] = 0 # Set masked values to 0

        sp.pl.dotplot(Top_gene_subset, 
                      linked_dict_symbols, 
                      gene_symbols = 'feature_name',
                      groupby = sample_group, 
                      show = False, 
                      dendrogram = False, 
                      categories_order = row_order)
        
        if sample_group == 'Group':
        
            sp.pl.dotplot(Top_gene_subset, 
                        linked_dict_symbols, 
                        gene_symbols = 'feature_name',
                        groupby = 'ann_finest_level', 
                        show = False, 
                        dendrogram = False, 
                        categories_order = group_order)

##################################################################

# Save combined embedings and normalized counts
combined_emb.write_h5ad(f'{directory}/combined_emb_normalized.h5ad')