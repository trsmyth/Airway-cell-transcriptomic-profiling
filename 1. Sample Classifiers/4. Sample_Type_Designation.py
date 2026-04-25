from transformers import set_seed, pipeline, AutoTokenizer, AutoModelForSequenceClassification
from Prediction_Functions import sliding_chunked_prediciton
import pandas as pd
import numpy as np
import pyreadr

set_seed(12345)

cat_code_directory = "./Output"
model_directory = "./Output/Cell_classifier_model"

input_directory = "./Input"
output_directory = "./Output"

#####################################################################

# Import the text data and sample ID information
dataset = pyreadr.read_r(f"{cat_code_directory}/Training_Data.RData")

# Define dataset
dataset = dataset["Training_Data"]

# Create category codes for sample IDs
dataset['Code'] = dataset['ID'].astype('category').cat.codes
dataset['Code'] = 'LABEL_' + dataset['Code'].astype(str)

# Map unique cat codes to sample IDs
cat_codes = np.unique(dict(zip(dataset['Code'], dataset['ID'])))
cat_codes = [code for code in cat_codes]

################################################################################################################################
################################################################################################################################
################################################################################################################################

# Import the total sample metadata dataset
total_dataset = pyreadr.read_r(f"{input_directory}/ARCHS4v2.6_Samples.RData")
total_dataset = total_dataset["Samples"]

# Check num_labels if CUDA runtime error raised
fine_tuned_model = AutoModelForSequenceClassification.from_pretrained(model_directory, 
                                                                      num_labels = len(cat_codes[0]), 
                                                                      local_files_only = True).to('cuda')

tokenizer = AutoTokenizer.from_pretrained(model_directory)

def tokenize(examples):

    return tokenizer(examples, add_special_tokens = False)

#####################################################################

'''
A small number of samples are over 512 tokens in length, meaning they would either 
have to be truncated, likely losing the context needed to correctly classify the 
sample, or the entire dataset would have to be run through the sliding chunk prediction 
which can only predict one sample at a time rather than batched for parallel processing. 
Below, metadata is tokenized and each example is checked for tokenized length. 
If tokenized length is within the BERT token limit (plus a little buffer), the untokenized 
metadata is added to a list for batched classification using the pipeline. Samples with 
metadata above the 512 limit are added to a list for a sliding chunk prediction which 
cannot batch inputs. This reduces prediction time from ~1h for ~900,000 samples to ~10m
with a chunk size of 150.

Token length sorting: 34s
Pipeline: 531s
Sliding chunked prediciton: 31s
'''

# Tokenize input data and isolate input_ids for sorting by length
tokenized_data = tokenize(total_dataset["Combined_Methods"].values.tolist())['input_ids']

pipeline_input = []
pipeline_sample = []

sliding_chunk_input = []
sliding_chunk_sample = []

# If tokenized data is above 500 tokens, add text data to a list for 
# sliding chunked prediction. 500 is chosen to provide some buffer 
# for special/start/end token additions in pipeline tokenizer
for num, tokens in enumerate(tokenized_data):
    if len(tokens) <= 500:
        pipeline_input.append(total_dataset["Combined_Methods"][num])
        pipeline_sample.append(total_dataset["sample"][num])

    else:
        sliding_chunk_input.append(total_dataset["Combined_Methods"][num])
        sliding_chunk_sample.append(total_dataset["sample"][num])

#####################################################################

import time

start_time = time.time()

batch_size = 150

# Define a text classification pipeline
clf = pipeline("text-classification", 
               model = fine_tuned_model, 
               tokenizer = tokenizer, 
               top_k = 1, 
               device = 0, 
               batch_size = batch_size)

# Run the prediction pipeline
Pipeline_sample_prediction = clf(pipeline_input)

end_time = time.time()
print(f'Prediction took {end_time - start_time}')

Pipeline_prediction = []

for pipeline_dict in Pipeline_sample_prediction:

    label = [d.get('label') for d in pipeline_dict if 'label' in d] # Isolate the predicted label
    Pipeline_prediction.append(cat_codes[0][label[0]]) # Convert the label to the sample type matching that label

######################

max_length = 400
jump_length = 100

# Get top predictions and instances of multiple predictions
sample_prediction, multiple_predictions = sliding_chunked_prediciton(tokenizer = tokenizer, 
                                                                     max_length = max_length, 
                                                                     jump_length = jump_length, 
                                                                     dataset = sliding_chunk_input,
                                                                     prediciton_pipeline = clf,
                                                                     prediction_keys = cat_codes, 
                                                                     default_class = 'Not Pulmonary Epithelium')

######################

Total_prediction = Pipeline_prediction + sample_prediction
Total_multiple_predictions = (['Single Class Predicted'] * len(Pipeline_prediction)) + multiple_predictions
Samples = pipeline_sample + sliding_chunk_sample

# Sort total dataset to match order of predictions
total_dataset['sample'] = total_dataset['sample'].astype(pd.CategoricalDtype(categories = Samples, ordered = True))
total_dataset_sorted = total_dataset.sort_values(by = 'sample')

######################

data = {'sample' : total_dataset_sorted["sample"], 
        'series_id' : total_dataset_sorted["series_id"], 
        "Prediction": Total_prediction,
        "Multiple_Predictions" : Total_multiple_predictions, 
        "Combined_Methods" : total_dataset_sorted["Combined_Methods"].values.tolist()}

Prediction = pd.DataFrame(data)

remove_categories = ['Not Pulmonary Epithelium', 'Cancerous Pulmonary', 'Cell Line', 'Organoid']
Prediction = Prediction[~Prediction["Prediction"].isin(remove_categories)]
Prediction = Prediction.reset_index()

columns = ['sample', 'series_id', 'Prediction', 'Multiple_Predictions', 'Combined_Methods']
Prediction = Prediction[columns]

Prediction.to_csv(f'{output_directory}/Prediction.csv')