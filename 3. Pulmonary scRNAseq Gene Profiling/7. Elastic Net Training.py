import time
import pickle
import cupy as cp
import numpy as np
import scipy.sparse
import pandas as pd
import anndata as ad
from sklearn import metrics
import matplotlib.pyplot as plt
from cuml.preprocessing import StandardScaler
from cuml.linear_model import LogisticRegression
from cuml.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, balanced_accuracy_score

%matplotlib

np.random.seed(12345)
cp.random.seed(12345)

############################

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

############################

directory = './Output'

h5_file = 'combined_emb_normalized.h5ad'

h5_data = ad.read_h5ad(f'{directory}/{h5_file}')
var = h5_data.var

# Set group/ann_finest_level to numeric cat codes
metadata = pd.DataFrame(h5_data.obs)
metadata['group_cat_code'] = metadata['Group'].astype('category').cat.codes
metadata['ann_cat_code'] = metadata['ann_finest_level'].astype('category').cat.codes

#####################

# Get names of genes which have expression in > 25% of samples in given group
Gene_index = []

for grouping in ['Basal', 'Nasal_secretory', 'Non_nasal_secretory', 'Ciliated']:

    h5_data_subset = h5_data[h5_data.obs.Group == grouping]
    percent_counts = group_expression_metrics(h5_data_subset, 'counts', 'ann_finest_level')
    Gene_index.append(percent_counts.loc[:, percent_counts.max(axis = 0) > 25].columns.tolist())

#####################

# Split data 70/30 as train and test datasets. 
X_train, X_test, y_train, y_test = train_test_split(h5_data.X, 
                                                    metadata['ann_finest_level'].astype('category').cat.codes, 
                                                    test_size = 0.3, 
                                                    random_state = 12345)

scaler = StandardScaler()
fit = scaler.fit(X_train) # Fit scaler
X_train_scaled = fit.transform(X_train) # Transform train data
X_test_scaled = fit.transform(X_test) # Transform test data

del X_train, X_test, h5_data # delete unneeded data to free memory

#####################

for i, grouping in enumerate(['All', 'Basal', 'Nasal_secretory', 'Non_nasal_secretory', 'Ciliated']):

    train_subset = metadata.loc[y_train.index.to_numpy()]
    test_subset = metadata.loc[y_test.index.to_numpy()]

    if grouping == 'All':
        group_category = 'Group'
        metadata_subset = metadata.copy()        
        category_map = dict(zip(metadata[group_category], metadata['group_cat_code']))

        X_train_subset = X_train_scaled
        y_train_subset = train_subset['group_cat_code']

        X_test_subset = X_test_scaled
        y_test_subset = test_subset['group_cat_code']

    # Isolate train and test count data and labels for samples related to current group
    else:
        group_category = 'ann_finest_level'
        metadata_subset = metadata[metadata['Group'] == grouping].copy()
        category_map = dict(zip(train_subset[group_category], train_subset['ann_cat_code']))

        train_subset_index = np.where(train_subset['Group'] == grouping)[0].tolist()
        X_train_subset = X_train_scaled[np.ix_(train_subset_index, Gene_index[i-1])]
        y_train_subset = train_subset.iloc[train_subset_index]['ann_cat_code']

        train_subset_index = np.where(test_subset['Group'] == grouping)[0].tolist()
        X_test_subset = X_test_scaled[np.ix_(train_subset_index, Gene_index[i-1])]
        y_test_subset = test_subset.iloc[train_subset_index]['ann_cat_code']

    ####################################################################################

    multiclass_model = LogisticRegression(C = 1, 
                                          l1_ratio = 0.5,
                                          max_iter = 5000,
                                          class_weight = 'balanced', 
                                          penalty = 'elasticnet',
                                          solver = 'qn')

    start = time.perf_counter()

    # Train the multiclass model
    multiclass_model.fit(X_train_subset, y_train_subset)

    end = time.perf_counter()
    print(f'Fitting model took {((end - start)/60):.2f} minutes')

    #########################################################################

    Coefficient_results = pd.DataFrame(multiclass_model.coef_)

    if grouping == 'All':
        Coefficient_results.columns = var.index

    else:
        Coefficient_results.columns = var.index[Gene_index[i-1]]

    # If more than two classes are present, coefficients for each group are presented
    # If two classes, coefficients are for the positive class (class 1)
    if len(np.unique(metadata_subset[group_category])) > 2:
        Coefficient_results.index = np.unique(metadata_subset[group_category])

    else:
        Coefficient_results.index = [k for k, v in category_map.items() if v == multiclass_model.classes_[1]]

    Predictive_cols = (Coefficient_results == 0).all(axis = 0)
    Predictive_cols = Predictive_cols.index[Predictive_cols == False]
    Predictive_genes = Coefficient_results[Predictive_cols]
    
    if grouping == 'All':
        print(f'{grouping} grouping had {Predictive_genes.shape[1]} genes with non-zero coefficients.')

    else:
        print(f'{grouping} grouping had {Predictive_genes.shape[1]} out of {len(Gene_index[i-1])} genes with non-zero coefficients.')

    Predictive_genes.to_csv(f'{directory}/LR Model Gene Coefficients/{multiclass_model.l1_ratio}_Predictive_genes_{grouping}.csv')

    with open(f'{directory}/LR Model Saves/{multiclass_model.l1_ratio}_logistic_regression_model_{grouping}.sav', 'wb') as file:
        pickle.dump(multiclass_model, file)

    #########################################################################
    
    # Get predictions for test data
    predictions = multiclass_model.predict(X_test_subset)

    # Assemble actual and prediction results for test data
    Results = pd.DataFrame({'Actual' : cp.asarray(y_test_subset).get(), 'Prediction' : cp.asarray(predictions).get()})
    
    # Invert category map to map results
    inverted_dict = {value: key for key, value in category_map.items()}
    Results["Actual"] = Results["Actual"].map(inverted_dict)
    Results["Prediction"] = Results["Prediction"].map(inverted_dict)
    
    # Create and plot confusion matrix
    confusion_matrix = metrics.confusion_matrix(Results["Actual"], Results["Prediction"])
    fig, ax = plt.subplots(figsize = (10, 8))
    matrix_plot = metrics.ConfusionMatrixDisplay(confusion_matrix = confusion_matrix, display_labels = np.unique(metadata_subset[group_category]))
    matrix_plot.plot(ax = ax, xticks_rotation = 90)
    plt.show()

    #########################################################################

    # Calculate model accuracy
    accuracy_scores = []

    for group in np.unique(metadata_subset[group_category]):

        current_group = Results[Results["Actual"] == group]
        accuracy_scores.append(accuracy_score(current_group["Actual"], current_group["Prediction"]))

    final_accuracy = pd.DataFrame(np.round(accuracy_scores, 4), 
                                  index = np.unique(metadata_subset[group_category]), 
                                  columns = ['Accuracy']).to_csv(f'{directory}/LR Model Metrics/{grouping}_model_accuracy.csv')

    #####

    # Calculate model metrics
    precision = precision_score(y_true = Results["Actual"], y_pred = Results["Prediction"], average = 'macro')
    recall = recall_score(y_true = Results["Actual"], y_pred = Results["Prediction"], average = 'macro')
    f1 = f1_score(y_true = Results["Actual"], y_pred = Results["Prediction"], average = 'macro')
    balanced = balanced_accuracy_score(y_true = Results["Actual"], y_pred = Results["Prediction"])

    #####

    final_metrics = pd.DataFrame({'Precision' : np.round(precision, 4), 
                                  'Recall' : np.round(recall, 4), 
                                  'F1 Score' : np.round(f1, 4), 
                                  'Balanced Accuracy' : np.round(balanced, 4)}, 
                                  index = ['']).to_csv(f'{directory}/LR Model Metrics/{grouping}_model_metrics.csv')
