import pyreadr
import evaluate
import numpy as np
from datasets import Dataset
from datasets.dataset_dict import DatasetDict
from sklearn.model_selection import train_test_split
from Prediction_Functions import find_token_length_distribution
from transformers import set_seed, AutoTokenizer, AutoModelForSequenceClassification, TrainingArguments, Trainer

########################################################

def tokenize(examples):

    return tokenizer(examples["text"], padding = "max_length", truncation = True)

def compute_metrics(eval_pred):

    logits, labels = eval_pred # Isolate logits and labels for predictions
    predictions = np.argmax(logits, axis = -1) # convert the logits to their predicted class

    return metric.compute(predictions = predictions, references = labels)

########################################################

set_seed(12345)

input_directory = "./Input"
save_directory = "./Output/Healthy_control_classifier_model"

tokenizer = AutoTokenizer.from_pretrained("google-bert/bert-base-uncased")

########################################################

# Import the text data and sample ID information
dataset = pyreadr.read_r(f"{input_directory}/Healthy_control_designation.RData")

# Define dataset
dataset = dataset["Healthy_control_designation"]

# print(find_token_length_distribution(dataset["Combined_Methods"].values.tolist(), tokenizer))

########################################################

# Create category codes for sample IDs
dataset['Code'] = dataset['ID'].astype('category').cat.codes

# Map unique cat codes to sample IDs
cat_codes = np.unique(dict(zip(dataset['Code'], dataset['ID'])))

# Check num_labels if CUDA runtime error raised
model = AutoModelForSequenceClassification.from_pretrained("google-bert/bert-base-uncased", 
                                                           num_labels = len(cat_codes[0])).to('cuda')

# Create train/test split
X_train, X_test, y_train, y_test = train_test_split(dataset["Combined_Methods"],
                                                    dataset["Code"], 
                                                    test_size = 0.25, 
                                                    shuffle = True)

# Create train/test dict
dataset = {'train': Dataset.from_dict({'label' : y_train, 'text' : X_train}),
           'test': Dataset.from_dict({'label' : y_test, 'text' : X_test})}

# Convert to dataset dict object
dataset = DatasetDict(dataset)

metric = evaluate.load("accuracy")

# Tokenize data
dataset = dataset.map(tokenize, batched = True)

# Define training arguments
training_args = TrainingArguments(num_train_epochs = 5,
                                  output_dir = "Models - Classifier",
                                  eval_strategy = "epoch",
                                  push_to_hub = False, 
                                  save_total_limit = 2, # Only keep best and latest model
                                  load_best_model_at_end = True, # Load best model at end of each epoch
                                  save_strategy = "epoch") # Save after every epoch

# Define trainer
trainer = Trainer(model = model,
                  args = training_args,
                  train_dataset = dataset["train"],
                  eval_dataset = dataset["test"],
                  compute_metrics = compute_metrics)

trainer.train() # Train the model
model.save_pretrained(save_directory) # Save the model
tokenizer.save_pretrained(save_directory) # Save the tokenizer