import sys
import scvi
import torch
import numpy as np
import pandas as pd
import scanpy as sp
import scarches as sca
import decoupler as dc
from scipy import stats

# https://docs.scarches.org/en/latest/hlca_map_classify.html

torch.cuda.is_available()
torch.set_float32_matmul_precision("high")

columns_to_retain = ['dataset', 
                     'donor_id', 
                     'development_stage', 
                     'sex', 
                     'tissue',
                     'ann_finest_level',
                     'anatomical_region_ccf_score',
                     '_scvi_batch', 
                     '_scvi_labels']

###############################################################

directory = './Output'

scANVI_model = scvi.model.SCANVI.load(f"{directory}/total_scanvi_model.pt")
Airway_data = sp.read_h5ad(f"{directory}/total_scanvi_model.pt/adata.h5ad")
del Airway_data.uns, Airway_data.varm, Airway_data.obsp, Airway_data.obsm

Airway_data.obsm["X_scANVI"] = scANVI_model.get_latent_representation(Airway_data)

labels = Airway_data.obs[['ann_finest_level', 'ann_level_1', 'ann_level_2', 
                         'ann_level_3', 'ann_level_4', 'ann_level_5']]

Airway_data.obs = Airway_data.obs.loc[:, columns_to_retain]
Airway_data.obs["ref_or_query"] = "ref"
Airway_data.X = Airway_data.layers['counts'] # Set X to raw counts

###############################################################

surgery_model = scvi.model.SCANVI.load(f"{directory}/surgery_model.pt")
adata_query = sp.read_h5ad(f"{directory}/surgery_model.pt/adata.h5ad")
del adata_query.uns, adata_query.obsm

'''
I noticed some samples which received ciliated cell labels after label transfer were marked
as digestive system tissue samples at a later stage of analysis and wanted to confirm these 
samples were not improperly included. This demonstrates that these ciliated samples are more 
specifically nasopharynx samples and should be retained.
'''

print(np.unique(adata_query.obs['tissue_general']))
Tissue_check = adata_query.obs[adata_query.obs['tissue_general'] == 'digestive system']
Tissue_check = Tissue_check.loc[:, ['tissue', 'tissue_general']]

adata_query.obsm["X_scANVI"] = surgery_model.get_latent_representation(adata_query)

adata_query.obs["tissue"] = adata_query.obs["tissue_general"]
adata_query.obs["anatomical_region_ccf_score"] = None
adata_query.obs = adata_query.obs.loc[:, columns_to_retain]
adata_query.obs["ref_or_query"] = "query"

###############################################################

def label_transfer_probs(query_adata, query_adata_emb, label_keys, knn_model, ref_adata_obs, region_keys = None):
    
    '''
    Modified implementation of scArches label transfer.
    This function calculates the class probabilities for each class rather than 
    returing only the majority class and corresponding uncertanties. Additionally,
    if a region_key is included, the median region score, which corresponds to where
    the reference sample was isolated, is calculated.

    Input:
    query_adata: Query data anndata object
    query_adata_emb: Query data embeddings. Expects the name of the anndata layer containing embeddings.
    label_keys: Name of column containing class labels for transfer.
    knn_model: KNN model to use for label transfer
    ref_adata_obs : Reference data obs anndata object, i.e. anndata.obs
    region_keys: Name of column containing region labels for transfer, if desired.

    Output: Pandas data frame with resulting class probabilities for each specified label_keys.

    '''

    # Calculate distances and indices of neighbors
    query_emb = query_adata.obsm[query_adata_emb]
    distances, indices = knn_model.kneighbors(X = query_emb)

    # Calculate Gaussian weights as per scArches logic
    stds = np.std(distances, axis = 1)
    stds = (2.0 / stds) ** 2
    stds = stds.reshape(-1, 1)

    # Exponential kernel for weighted distance
    # I.e. weigh data points based on distance from center with exponential decrease based on distance
    weights = np.exp(-np.true_divide(distances, stds))
    weights = weights / np.sum(weights, axis = 1, keepdims = True)

    # Aggregate weights by class
    y_train_labels = ref_adata_obs[label_keys].values
    unique_labels = np.unique(y_train_labels).tolist()

    probs = np.empty((len(query_adata.obs), len(unique_labels)))

    # Get region ccf scores
    if region_keys != None:
        y_train_ccf = ref_adata_obs[region_keys].values
        ccf_median = []

    for i in range(len(weights)):

        # Get labels of neighbors for cell
        neighbor_labels = y_train_labels[indices[i]]

        # Calculate mean and median of neighbor ccf scores
        if region_keys != None:
            neighbor_ccf = y_train_ccf[indices[i]]
            ccf_median.append(np.median(neighbor_ccf))

        for j, label in enumerate(neighbor_labels):
            probs[i, unique_labels.index(label)] += weights[i, j]

    # Calculate uncertanty and determine label of highest prediction
    probs_df = pd.DataFrame(probs, columns = unique_labels)
    probs_df.insert(0, column = 'Label', value = probs_df.idxmax(axis = 1))
    probs_df.insert(1, column = 'Uncertanty', value = np.round(1.0 - probs.max(axis = 1), 6))
    probs_df.loc[probs_df['Uncertanty'] <= 0, 'Uncertanty'] = 0

    if region_keys != None:
        probs_df.insert(2, column = 'median_ccf_score', value = ccf_median)

    return probs_df

###############################################################

knn_transformer = sca.utils.knn.weighted_knn_trainer(train_adata = Airway_data,
                                                     train_adata_emb = "X_scANVI",  # location of joint embedding
                                                     n_neighbors = 50)

probs = label_transfer_probs(query_adata = adata_query,
                             query_adata_emb = "X_scANVI",  # location of embedding
                             label_keys = "ann_finest_level",  # obs column name for transfer labels
                             knn_model = knn_transformer,
                             ref_adata_obs = Airway_data.obs, 
                             region_keys = 'anatomical_region_ccf_score')

probs.to_csv(f'{directory}/Label_transfer_results.csv')

###############################################################

labels = probs['Label']
uncert = probs['Uncertanty']

# Use 20% uncertanty cutoff. I.e. > 80% certanty or above for cell type classification
labels = labels.mask(uncert > 0.2, "Unknown")

adata_query.obs['ann_finest_level'] = labels.tolist()
adata_query.obs['anatomical_region_ccf_score'] = probs['median_ccf_score'].tolist()

adata_query = adata_query[adata_query.obs['ann_finest_level'] != 'Unknown']

combined_emb = sp.concat((Airway_data, adata_query), index_unique = None, join = "outer")
combined_emb.var = Airway_data.var
del combined_emb.layers['counts']

combined_emb.write_h5ad(f'{directory}/combined_emb_raw.h5ad')

###############################################################

combined_emb = sp.read_h5ad(f"{directory}/combined_emb_raw.h5ad")

import seaborn as sns

sp.pp.neighbors(combined_emb, use_rep = "X_scANVI")
sp.tl.umap(combined_emb)

# Calculate leiden clustering with designated resolutions and plot umap
for resolution in [0.1, 0.25]:

    sp.tl.leiden(combined_emb, flavor = "igraph", n_iterations = 2, resolution = resolution, key_added = f'{resolution}_leiden')
    sp.pl.umap(combined_emb, color = [f'{resolution}_leiden'], frameon = False, ncols = 1)

sp.pl.umap(combined_emb, color = ["dataset"], frameon = False, legend_loc = None)

sp.pl.umap(combined_emb, color = ["ref_or_query"], frameon = False)

############################

# Plot umap for each finest level annotation
sp.pl.umap(combined_emb, color = ["ann_finest_level"], 
           frameon = False, legend_loc = "on data", 
           legend_fontsize = 4)

# Create a color palette for better visualization of annotations
combined_palette = sns.color_palette("tab10") + sns.color_palette("Accent") + sns.color_palette("Dark2") + sns.color_palette("Set2") + sns.color_palette("Set1")

# Set finest level annotation colors to specified colors
combined_emb.uns['ann_finest_level_colors'] = combined_palette[0:len(combined_emb.uns['ann_finest_level_colors'])]

# Replot with the new color scheme
sp.pl.umap(combined_emb, color = ["ann_finest_level"], 
           frameon = False, legend_loc = "on data", 
           legend_fontsize = 4)

############################

# Plot umap with anatomical region ccf score as color
sp.pl.umap(combined_emb, 
           color = ["anatomical_region_ccf_score"], 
           frameon = True)