import torch
import scvi
import numpy as np
import scanpy as sp

torch.cuda.is_available()
torch.set_float32_matmul_precision("high")

input_directory = './Input'
output_directory = './Output'
h5_file = 'Airway_data.h5ad'

Airway_data = sp.read_h5ad(f'{input_directory}/{h5_file}')

# Select highly variable genes
sp.experimental.pp.highly_variable_genes(Airway_data,
                                         n_top_genes = 5000,
                                         layer = 'counts',
                                         flavor = 'pearson_residuals', 
                                         batch_key = 'dataset')

Airway_data = Airway_data[:, Airway_data.var['highly_variable']].copy()

########################################################

# Set up scvi anndata. This uses raw count data with batch and covariate information.
scvi.model.SCVI.setup_anndata(Airway_data, 
                              layer = "counts", 
                              batch_key = "dataset", 
                              categorical_covariate_keys = ['sex'])

# Create a scvi model and train it on the airway data with gpu acceleration
model = scvi.model.SCVI(Airway_data, 
                        n_layers = 2, 
                        n_latent = 30, 
                        encode_covariates = True,
                        deeply_inject_covariates = False,
                        use_layer_norm = "both",
                        use_batch_norm = "none",
                        gene_likelihood = "nb")

model.view_anndata_setup(Airway_data)

early_stopping_kwargs = {
    "early_stopping": True,
    "early_stopping_monitor": "elbo_train",
    "early_stopping_patience": 5,
    "check_val_every_n_epoch": 1,
}

model.train(max_epochs = 500, batch_size = 1024, **early_stopping_kwargs)
model.save(f"{output_directory}/total_scvi_model.pt", save_anndata = False)

####

history = model.history_
elbo = history["elbo_train"]

import matplotlib.pyplot as plt

plt.figure(figsize = (10, 6))

plt.plot(np.arange(len(elbo)), elbo, label = "ELBO")
plt.xlabel("Epoch")
plt.ylabel("ELBO")
plt.title("ELBO during scVI Training")
plt.legend()
plt.show()

#######################

early_stopping_kwargs = {
    "early_stopping": True,
    "early_stopping_monitor": "elbo_train",
    "early_stopping_patience": 10,
    "check_val_every_n_epoch": 1,
}

# Create the scanvi model from the trained scvi model
scanvi_model = scvi.model.SCANVI.from_scvi_model(model, 
                                                 adata = Airway_data,
                                                 labels_key = "ann_finest_level",
                                                 unlabeled_category = "Unknown")

scanvi_model.view_anndata_setup(Airway_data)

# Train the scanvi model
scanvi_model.train(max_epochs = 500, n_samples_per_label = 50, **early_stopping_kwargs)
scanvi_model.save(f"{output_directory}/total_scanvi_model.pt", save_anndata = True)

####

history = scanvi_model.history_
elbo = history["elbo_train"]

plt.figure(figsize = (10, 6))

plt.plot(np.arange(len(elbo)), elbo, label = "ELBO")
plt.xlabel("Epoch")
plt.ylabel("ELBO")
plt.title("ELBO during scANVI Training")
plt.legend()
plt.show()