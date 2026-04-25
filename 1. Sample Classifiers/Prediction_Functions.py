from collections import Counter
from tqdm import tqdm

##################################################################################################################################

def find_token_length_distribution(dataset, tokenizer):

    # BERT has a 512 token maximum
    # Longformer has a 4096 token maximum
    # ModernBERT has a 8192 token maximum
    token_length = []
    Within_BERT_Limit = 0
    over_512 = 0
    over_1000 = 0
    over_2000 = 0
    over_4096 = 0
    over_8192 = 0

    for example in dataset:

        number = len(tokenizer(example)['input_ids'])

        token_length.append(number)

        if number <= 512:
            Within_BERT_Limit += 1
        if number > 512:
            over_512 += 1
        if number > 1000:
            over_1000 += 1
        if number > 2000:
            over_2000 += 1
        if number > 4096:
            over_4096 += 1
        if number > 8192:
            over_8192 += 1

    Output = str(f"Total dataset: {len(dataset)}\n"
                f"Within BERT Limit: {Within_BERT_Limit}, {(Within_BERT_Limit/len(dataset)):.2%}\n"
                f"Over 512: {over_512}, {(over_512/len(dataset)):.2%}\n"
                f"Over 1000: {over_1000}, {(over_1000/len(dataset)):.2%}\n"
                f"Over 2000: {over_2000}, {(over_2000/len(dataset)):.2%}\n"
                f"Over 4096: {over_4096}, {(over_4096/len(dataset)):.2%}\n"
                f"Over 8192: {over_8192}, {(over_8192/len(dataset)):.2%}")

    return Output

##################################################################################################################################

def sliding_chunked_prediciton(tokenizer, 
                               dataset, 
                               prediciton_pipeline, 
                               prediction_keys, 
                               default_class,
                               max_length = 400, 
                               jump_length = 100):

    # Convert text data to list of strings
    tokens = tokenizer(dataset, add_special_tokens = False)["input_ids"]

    progress_bar = tqdm(total = len(tokens))

    top_prediction = []
    multiple_prediction_index = []

    # For each tokenized string
    for num, single_str in enumerate(tokens):

        # If the length of the tokenized string is greater than the max length
        if len(single_str) > max_length:

            # Define the chunking window based on the length of the string andf the jump length
            chunked_windows = [single_str[i : i + max_length] for i in range(0, len(single_str), jump_length)]
            
            prediction = [] # List to store predictions
            score = [] # List to store prediction scores

            # For each window in the chunked windows
            for window in chunked_windows:
                inputs = tokenizer.decode(window, skip_special_tokens = True) # Decode the tokenized string
                outputs = prediciton_pipeline(inputs) # Feed decoded string into the prediction pipeline
                prediction.append(outputs[0][0]['label']) # Save the top prediction
                score.append(outputs[0][0]['score']) # Save the top prediction score

            Label_count = Counter(prediction) # Count the number of each prediction
            Count_keys = Label_count.keys() # Isolate the prediction keys

            # If only one label is predicted, save the prediction and mark it as having a single class
            if len(Label_count) == 1:    

                top_prediction.append(prediction_keys[0][prediction[score.index(max(score))]])
                multiple_prediction_index.append('Single Class Predicted')

            # Else, count the number of each prediction
            else:
                
                Count_occurance = []
                Average_prediction = []
                Position_key = []

                # For each found prediction
                for key in Count_keys:
                    Position_key.append(key) # Save the key at the position
                    Count_occurance.append(Label_count[key]) # Save the number of counts at the position
                    
                    # Identify the index of predictions which match the current key
                    value_index = [index for index, value in enumerate(prediction) if value == key] 
                    # Identify the scores of predictions which match the above indices
                    score_value = [value for index, value in enumerate(score) if index in value_index]

                    # If the score value for the current key is greater than 0, calculate the average score
                    if len(score_value) > 0:
                        Average_prediction.append(sum(score_value) / len(score_value))

                    # Else, the average prediction for the current key is 0
                    else:
                        Average_prediction.append(0)

                sorted_list = sorted(zip(Count_occurance, Average_prediction, Position_key), key = lambda x: (-x[0], -x[1]))

                if prediction_keys[0][sorted_list[0][2]] != default_class:
                    top_prediction.append(prediction_keys[0][sorted_list[0][2]])

                else:
                    top_prediction.append(prediction_keys[0][sorted_list[1][2]])

                multiple_prediction_index.append([num, sorted_list])
        
        # If the tokenized length is less than the token limit, run the string through the pipeline normally
        else:
            multiple_prediction_index.append('Single Class Predicted')
            inputs = tokenizer.decode(single_str, skip_special_tokens = True)
            outputs = prediciton_pipeline(inputs)
            top_prediction.append(prediction_keys[0][outputs[0][0]['label']])

        # Update progress bar by 1
        progress_bar.update(1)

    # Return the top prediction and the details surrounding multiple predictions
    return top_prediction, multiple_prediction_index