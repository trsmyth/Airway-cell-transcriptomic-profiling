import torch
import scanpy as sp
import anndata as ad

torch.cuda.is_available()
torch.set_float32_matmul_precision("high")

input_directory = './Input'
h5_file = 'Human Lung Cell Atlas Core.h5ad'

Data = sp.read_h5ad(f'{input_directory}/{h5_file}')

del Data.uns, Data.obsm, Data.obsp

# Set count data to raw count data
Data.X = Data.raw.X 

# List gene and cell types to retain
separator = "|" 
gene_types_of_interest = ['protein_coding', 'lncRNA']

Cells_of_interest = ['Basal', 'Club', 'Goblet', 'Multiciliated', 
                     'SMG', 'Suprabasal', 'Deuterosomal', 'Hillock', 
                     'AT0', 'AT1', 'AT2', 'Ionocyte', 'pre-TB secretory']

# Add separator to cell and gene types
gene_types_of_interest = separator.join(gene_types_of_interest)
Cells_of_interest = separator.join(Cells_of_interest)

# Isolate samples and genes according to above lists
Airway_data = Data[Data.obs['ann_finest_level'].str.contains(Cells_of_interest), 
                   Data.var['feature_type'].str.contains(gene_types_of_interest)].copy()

with open('Gene_names_in_HLCA.txt', "w") as file_handler:
    for item in Airway_data.var.index:
        file_handler.write(f"{item}")

Airway_data = ad.AnnData(X = Airway_data.X, 
                         obs = Airway_data.obs, 
                         var = Airway_data.var)

Airway_data.layers["counts"] = Airway_data.X.copy()

Airway_data.write_h5ad(f'{input_directory}/Airway_data.h5ad')