from transformers import set_seed, pipeline, AutoTokenizer, AutoModelForSequenceClassification
import pandas as pd
import numpy as np
import GEOparse

set_seed(12345)

input_directory = "./Input"
output_directory = "./Output"
model_directory = "./top_model"

#####################################################################

Prediction = pd.read_csv(f'{input_directory}/Prediction.csv')

cat_codes = dict(zip(['LABEL_0', 'LABEL_1'], ['Healthy_Control', 'Not_Healthy_Control']))

#####################################################################

# Check num_labels if CUDA runtime error raised
fine_tuned_model = AutoModelForSequenceClassification.from_pretrained(model_directory, 
                                                                      num_labels = len(cat_codes), 
                                                                      local_files_only = True).to('cuda')

tokenizer = AutoTokenizer.from_pretrained(model_directory)

def tokenize(examples):

    return tokenizer(examples, add_special_tokens = False)

#####################################################################

import time

start_time = time.time()

batch_size = 150

clf = pipeline("text-classification", 
               model = fine_tuned_model, 
               tokenizer = tokenizer, 
               top_k = 1, 
               device = 0, 
               batch_size = batch_size)

Pipeline_sample_prediction = clf(Prediction['Combined_Methods'].values.tolist())

end_time = time.time()
print(f'Prediction took {end_time - start_time}')

Pipeline_prediction = []

for pipeline_dict in Pipeline_sample_prediction:

    label = [d.get('label') for d in pipeline_dict if 'label' in d] # Isolate the predicted label
    Pipeline_prediction.append(cat_codes[label[0]]) # Convert the label to the sample type matching that label

######################

data = {'sample' : Prediction["sample"], 
        'series_id' : Prediction["series_id"], 
        'Cell_Type' : Prediction["Prediction"],
        "Prediction": Pipeline_prediction,
        "Combined_Methods" : Prediction["Combined_Methods"].values.tolist()}

Prediction = pd.DataFrame(data)

Prediction.to_csv(f'{output_directory}/Prediction.csv')

#####################################################################

columns = ['sample', 'series_id', 'Cell_Type', 'Prediction', 'Combined_Methods']

Prediction = Prediction[Prediction["Prediction"] != 'Not_Healthy_Control']
Prediction = Prediction[columns]
Prediction = Prediction.reset_index()

# Isolate unique experimental series IDs
unique_series_ids = np.unique(Prediction["series_id"])

series_id = []

# Some samples have two series id designations. Split multi-series designations so each can be searched.
for series in unique_series_ids:
    split_series = series.split(',')
    for split in split_series:
        series_id.append(split)

gse_output = []

for series in series_id:
    
    gse = GEOparse.get_GEO(geo = series, destdir = "./GEO") # Get GEO data for current series
    number_samples = len(gse.metadata["sample_id"]) # Get number of samples linked to the series
    
    GSE_Metadata = []
    GSE_categories = ['pubmed_id', 'contact_name', 'contact_email', 'geo_accession', 'title', 'summary', 'overall_design']
 
    for category in GSE_categories: # For each desired GEO data category
        if category in gse.metadata.keys(): # If the category exists for the current series
            if len(gse.metadata[category]) > 1:
                GSE_Metadata.append([" ".join(gse.metadata[category])] * number_samples) # Join data and set category for each sample
            else:
                GSE_Metadata.append(gse.metadata[category] * number_samples) # Set category for each sample
        else:
            GSE_Metadata.append([" "] * number_samples) # Set category to blank for each sample

    GSM_categories = ['title', 'geo_accession', 
                      'source_name_ch1', 'characteristics_ch1', 
                      'treatment_protocol_ch1', 'growth_protocol_ch1', 
                      'extract_protocol_ch1', 'description']
    
    for category in GSM_categories: # For each specified category
        sample = []
        for gsm_name, gsm in gse.gsms.items(): # For each sample name and its meatdata
            if category in gsm.metadata: # If the category is in the metadata
                if len(gsm.metadata[category]) > 1:
                    sample.append(" ".join(gsm.metadata[category])) # Join the metadata
                else:
                    sample.append(gsm.metadata[category]) # Isolate the metadata
            else:
                sample.append(" ") # Set category metadata to blank

        GSE_Metadata.append(sample) # Add sample metadata to gse metadata
    gse_output.append(GSE_Metadata) # Add gse metadata to output list

# Groups/transposes metadata from across all series together and merges grouped metadata into a flat list for each series.
gse_list = [[item for subsublist in sublist for item in subsublist] for sublist in [[item for item in sublist] for sublist in zip(*gse_output)]]

df_colnames = ['PMID', 'contact_name', 'contact_email', 'series_id',
               'experimental_series_title', 'summary', 'overall_design',
               'title', 'geo_accession', 'source_name_ch1', 
               'characteristics_ch1', 'treatment_protocol_ch1', 
               'growth_protocol_ch1', 'extract_protocol_ch1', 'description']

# Create a dict with column names as keys and gse_list metadata as values
data = dict(zip(df_colnames, gse_list))

GSE_data = pd.DataFrame(data)

GSE_data.to_csv(f'{output_directory}/GSE_data.csv')

#####################################################################

# Set geo_accession metadata to a categorical value to sort, drop na/duplicates, and reset index
GSE_data['geo_accession'] = pd.Categorical(GSE_data['geo_accession'], categories = Prediction['sample'].values, ordered = True)
GSE_data = GSE_data.sort_values('geo_accession')
GSE_data = GSE_data.dropna(subset = ['geo_accession'])
GSE_data = GSE_data.drop_duplicates('geo_accession', keep = 'first')
GSE_data = GSE_data.reset_index()

# Add cell type and prediction data to metadata
GSE_data['Cell_Type'] = Prediction['Cell_Type']
GSE_data['Prediction'] = Prediction['Prediction']

columns = ['PMID', 'series_id', 'geo_accession', 'Cell_Type', 'Prediction', 'contact_name', 
           'contact_email', 'experimental_series_title', 'summary', 'overall_design', 
           'title', 'source_name_ch1', 'characteristics_ch1', 'treatment_protocol_ch1', 
           'growth_protocol_ch1', 'extract_protocol_ch1', 'description']

GSE_data = GSE_data[columns]

GSE_data.to_csv(f'{output_directory}/Sample_Classification_Results.csv')