import torch
import scvi
import scanpy as sp
import seaborn as sns
import scarches as sca

torch.cuda.is_available()
torch.set_float32_matmul_precision("high")

input_directory = './Input'
output_directory = './Output'

scanvi_model = scvi.model.SCANVI.load(f"{output_directory}/total_scanvi_model.pt")
Airway_data = sp.read_h5ad(f"{output_directory}/total_scanvi_model.pt/adata.h5ad")

# Return mean and variance of the latent distribution
Airway_data.obsm["X_scANVI"] = scanvi_model.get_latent_representation(Airway_data)

plot = True

if plot == True:

    # Calculate neighbors, umap, and leiden clustering
    sp.pp.neighbors(Airway_data, use_rep = "X_scANVI")
    sp.tl.umap(Airway_data)

    for resolution in [0.1, 0.25]:

        sp.tl.leiden(Airway_data, flavor = "igraph", n_iterations = 2, resolution = resolution, key_added = f'{resolution}_leiden')
        sp.pl.umap(Airway_data, color = [f'{resolution}_leiden'], frameon = False, ncols = 1)

    sp.pl.umap(Airway_data,
            color = ["dataset"],
            frameon = True,
            legend_loc = 'none',
            ncols = 1)
    
    sp.pl.umap(Airway_data,
            color = ["ann_finest_level"],
            frameon = True,
            legend_loc = 'none',
            ncols = 1)
    
    combined_palette = sns.color_palette("tab10") + sns.color_palette("Accent") + sns.color_palette("Dark2") + sns.color_palette("Set2") + sns.color_palette("Set1")

    Airway_data.uns['ann_finest_level_colors'] = combined_palette[0:len(Airway_data.uns['ann_finest_level_colors'])]

    sp.pl.umap(Airway_data, color = ["ann_finest_level"], 
            frameon = True, legend_loc = "on data", 
            legend_fontsize = 4)

######################

query_data = sp.read_h5ad(f"{input_directory}/Over_18_airway_samples.h5ad")
query_data.var.index = query_data.var.feature_id

adata_query = sca.models.SCANVI.prepare_query_anndata(adata = query_data, 
                                                      reference_model = f"{output_directory}/total_scanvi_model.pt", 
                                                      inplace = False)

adata_query.layers["counts"] = adata_query.X.copy()
adata_query.obs['dataset'] = adata_query.obs['dataset_id']

surgery_model = sca.models.SCANVI.load_query_data(adata_query,
                                                  f"{output_directory}/total_scanvi_model.pt",
                                                  freeze_dropout = True)

surgery_model.registry_["setup_args"]

surgery_epochs = 500

early_stopping_kwargs_surgery = {
    "early_stopping_monitor": "elbo_train",
    "early_stopping_patience": 10,
    "early_stopping_min_delta": 0.001,
    "plan_kwargs": {"weight_decay": 0.0},
}

surgery_model.train(max_epochs = surgery_epochs, batch_size = 1024, **early_stopping_kwargs_surgery)
surgery_model.save(f"{output_directory}/surgery_model.pt", save_anndata = True)

history = surgery_model.history_
elbo = history["elbo_train"]

import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize = (10, 6))

plt.plot(np.arange(len(elbo)), elbo, label = "ELBO")
plt.xlabel("Epoch")
plt.ylabel("ELBO")
plt.title("ELBO during scARCHS Training")
plt.legend()
plt.show()