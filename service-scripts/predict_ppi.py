#!/usr/bin/env python
# Carla M. Mann
# 2025-06-23

import transformers
import datasets
import pandas as pd
import numpy as np
import torch
import sys
import json
import os
import Bio.SeqIO
import Bio.Seq
import Bio.SeqRecord
import ast
import datetime
import os
import warnings

# Suppress warnings when divide by zero while calculating metrics
#warnings.filterwarnings("ignore", category = RuntimeWarning)
# Suppress deprecation warnings
#warnings.filterwarnings("ignore", category = DeprecationWarning)
# Suppress the seqeval "0"/"1" are not recognized as labels error
#warnings.filterwarnings("ignore", category = UserWarning)
# Suppress warnings about modifying dataframes
pd.set_option("mode.chained_assignment", None)
# Suppress warnings about future deprecations
#warnings.filterwarnings("ignore", category = FutureWarning)

# Disable the progress bar so that it isn't displayed
transformers.utils.logging.disable_progress_bar()

# Handle reading in fasta files
def read_fasta(file):
    reject_list = []
    
    try:
        seqs = list(Bio.SeqIO.parse(file, "fasta"))
        
        # Strip any * from the end of protein sequences
        # And check if they are protein
        for record in seqs:
            clean_seq = str(record.seq).replace("*", "")
            record.seq = Bio.Seq.Seq(clean_seq)
            if not is_protein(record):
                reject_list.append(record.id)
            
    except Exception as e:
        err = str(e)
        sys.stderr.write(err)

    if len(reject_list) > 0:
        print(reject_list)
        sys.stderr.write("Warning: one or more input sequences appear to be DNA:\n")
        for item in reject_list:
            print(item, file=sys.stderr)

    return(seqs)

def is_protein(seq_record: Bio.SeqRecord) -> bool:
    # Define sets
    dna_letters = set("ACGTUN")
    protein_letters = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
    
    seq_letters = set(str(seq_record.seq).upper())
    
    if seq_letters <= dna_letters:
        # Only DNA letters → not protein
        return False
    elif seq_letters <= protein_letters:
        # Only protein letters (but not pure DNA letters) → protein
        return True
    else:
        # Mixed/ambiguous/invalid → conservative False
        return False
        
def get_seq_combinations(seqs1, seqs2):
    id_list = []
    s1_list = []
    s2_list = []
    l1_list = []
    l2_list = []

    for record1 in seqs1:
        for record2 in seqs2:
            id_list.append(record1.name + "_" + record2.name)
            s1_list.append(str(record1.seq).upper())
            s2_list.append(str(record2.seq).upper())
            l1 = "0" * len(str(record1.seq))
            l2 = "0" * len(str(record2.seq))
            l1_list.append(l1)
            l2_list.append(l2)

    rdf = pd.DataFrame(zip(id_list, s1_list, s2_list, l1_list, l2_list), columns = ["ids", "seq1", "seq2", "labels1", "labels2"])
    return(rdf)
    
def split_into_tokens(string):
    ret_list = list(string)
    ret_list = [s for s in string if s != " "]
    return(ret_list)

# Splits a string of numerical characters into individual int digits
def split_numeric_string_into_tokens(string):
    ret_list = list(string)
    ret_list = [int(s) for s in string if s!= " "]
    return(ret_list)
    
def preprocess_fn(examples, tokenizer):
    seq1 = [split_into_tokens(seq) for seq in examples["seq1"]]
    seq2 = [split_into_tokens(seq) for seq in examples["seq2"]]
    
    lab1 = [split_numeric_string_into_tokens(labs) for labs in examples["labels1"]]
    lab2 = [split_numeric_string_into_tokens(labs) for labs in examples["labels2"]]

    labels = [l1 + [-100] + l2 for l1, l2 in zip(lab1, lab2)]
    inputs = tokenizer(seq1, seq2, is_split_into_words = True)
    inputs["labels"] = labels
    return(inputs)
    
def filter_long_examples(example):
    '''Filter function applied by Huggingface dataset's "filter" function to kick out items from the dataset that are too long (>1024 token>
    '''
    return(len(example["input_ids"]) <= 1024)
    
def process_output(outputs, batch, tokenizer):
    '''Proccesses output from the model for predicting interfacial residues from two protein sequences.
    Reports the logits, so that performance can be examined using different thresholds after the fact.

    Parameters
    ----------
    outputs : torch tensor 
        The output tensor from a Huggingface model
    batch : 
        The output batch from a Huggingface model
    tokenizer: Huggingface tokenizer

    Returns
    -------
    ret_frame : pandas dataframe
        A dataframe with ten columns for protein sequence 1, protein sequence 2, the predicted 
        labels for protein 1, the predicted labels for protein 2, the actual ground truth labels 
        for protein 1, ground truth labels for protein 2, and the number of true positive, true 
        negative, false positive, and false negative predictions: ["Prot1", "Prot2", "Labels1", 
        "Labels2", "Actual1", "Actual2", "TP", "TN", "FP", "FN"]
    '''
    tp, tn, fp, fn = [], [], [], []

    # Get the predictions and compare to the batch labels
    logits = outputs.logits
    pos_probs = torch.sigmoid(logits[:, :, 1]) # Convert logits to probabilities
    pred_labels = (pos_probs * 10**3).round() / (10**3) # Round the logit to 3 decimal points for space
    
    # Move the result tensor to CPU.
    #
    # .tolist() rather than list(np.array(...)): under NumPy 2 the latter calls
    # Tensor.__array__(dtype=None, copy=False), which torch does not accept, so
    # every batch emitted a DeprecationWarning. It also yields np.int64 scalars,
    # and these values are appended straight into the returned frame (act1/act2
    # below are not .tolist()'d the way lab1/lab2 are) -- so they reached the
    # output as numpy scalars. .tolist() drops numpy from the path entirely and
    # gives plain ints; slicing and the -100 comparison behave identically.
    #
    batch_labels = [batch["labels"][i].detach().cpu().tolist() for i in range(len(batch["input_ids"]))]

    # Get the original input sequence from the tokenized sequence in the bath
    decoded = [tokenizer.decode(b) for b in batch["input_ids"]]
    prot1 = [d.split("<eos>")[0].replace("<cls>", "").replace(" ", "") for d in decoded]
    prot2 = [d.split("<eos>")[1].replace(" ", "") for d in decoded]

    # Replace unknown tokens with "X" amino acid
    prot1 = [p.replace("<unk>", "X") for p in prot1]
    prot2 = [p.replace("<unk>", "X") for p in prot2]

    
    pdf = pd.DataFrame(zip(pred_labels, batch_labels, prot1, prot2))
    pdf.columns = ["preds", "actual", "prot1", "prot2"]

    prot1_l, prot2_l = [], []
    lab1_l,  lab2_l  = [], []
    act1_l,  act2_l  = [], []

    for i in range(len(pdf)):
        # Make a temporary dataframe to hold the contents of the other dataframe
        tdf = pd.DataFrame(zip(pdf.preds[i], pdf.actual[i]), columns = ["preds", "actual"])

        # Get the actual and predicted labels
        preds = pdf.preds[i]
        actual = pdf.actual[i]

        # Get the label indices that are part of each sequence
        l1 = len(pdf.prot1[i])
        l2 = len(pdf.prot2[i])
        
        # Split the ground truth and predicted labels to correspond to the original sequences
        lab1 = preds[1 : 1 + l1]
        lab2 = preds[1 + l1 + 1 : 1 + l1 + 1 + l2]
        act1 = actual[1 : 1 + l1]
        act2 = actual[1 + l1 + 1 : 1 + l1 + 1 + l2]

        # Find the places that are a -100 (indicating separator/end of sequence),
        # and kick them out
        tdf_n100 = tdf[tdf.actual != -100]
        
        tp.append(len(tdf_n100[(tdf_n100.preds == 1) & (tdf_n100.actual == 1)]))
        tn.append(len(tdf_n100[(tdf_n100.preds == 0) & (tdf_n100.actual == 0)]))
        fp.append(len(tdf_n100[(tdf_n100.preds == 1) & (tdf_n100.actual == 0)]))
        fn.append(len(tdf_n100[(tdf_n100.preds == 0) & (tdf_n100.actual == 1)]))

        prot1_l.append(pdf.prot1[i])
        prot2_l.append(pdf.prot2[i])
        lab1_l.append(lab1.tolist())
        lab2_l.append(lab2.tolist())
        act1_l.append(act1)
        act2_l.append(act2)

    # Make dataframe
    ret_frame = pd.DataFrame(zip(prot1_l, prot2_l, lab1_l, lab2_l, act1_l, act2_l, tp, tn, fp, fn), 
                             columns = [  "Prot1",  "Prot2", 
                                        "Labels1", "Labels2", 
                                        "Actual1", "Actual2", 
                                        "TP", "TN", "FP", "FN"])
    return(ret_frame)

def process_output_s3(outputs, batch, tokenizer):
    '''Proccesses output from the model for predicting interfacial residues from two protein sequences.
    Reports the logits, so that performance can be examined using different thresholds after the fact.

    Parameters
    ----------
    outputs : torch tensor 
        The output tensor from a Huggingface model
    batch : 
        The output batch from a Huggingface model
    tokenizer : Huggingface tokenizer

    Returns
    -------
    ret_frame : pandas dataframe
        A dataframe with thirteen columns for protein sequence 1, protein sequence 2, protein sequnce 3, 
        predicted labels for protein 1, predicted labels for protein 2, predicted labels for protein 3, 
        actual ground truth labels for protein 1, ground truth labels for protein 2, ground truth labels 
        for protein 3, and the number of true positive, true negative, false positive, and false negative 
        predictions: ["Prot1", "Prot2", "Prot3", "Labels1", "Labels2", "Labels3", "Actual1", "Actual2", 
        "Actual3", "TP", "TN", "FP", "FN"]
    '''
    tp, tn, fp, fn = [], [], [], []

    # Get the predictions and compare to the batch labels
    logits       = outputs.logits
    pos_probs    = torch.sigmoid(logits[:, :, 1]) # Convert logits to probabilities
    pred_labels  = (pos_probs * 10**3).round() / (10**3)
    batch_labels = [batch["labels"][i].detach().cpu().tolist() for i in range(len(batch["input_ids"]))]

    # Get the input sequences
    decoded = [tokenizer.decode(b) for b in batch["input_ids"]]
    prot1 = [d.split("<eos>")[0].replace("<cls>", "").replace(" ", "") for d in decoded]
    prot2 = [d.split("<eos>")[1].replace(" ", "") for d in decoded]
    prot3 = [d.split("<eos>")[2].replace(" ", "") for d in decoded]

    prot1 = [p.replace("<unk>", "X") for p in prot1]
    prot2 = [p.replace("<unk>", "X") for p in prot2]
    prot3 = [p.replace("<unk>", "X") for p in prot3]

    # Find the places that are a -100, and kick them out of the preds and the batch labels
    pdf = pd.DataFrame(zip(pred_labels, batch_labels, prot1, prot2, prot3))
    pdf.columns = ["preds", "actual", "prot1", "prot2", "prot3"]

    prot1_l, prot2_l, prot3_l = [], [], []
    lab1_l,   lab2_l,  lab3_l = [], [], []
    act1_l,   act2_l,  act3_l = [], [], []

    for i in range(len(pdf)):
        # Make a temporary dataframe to hold the contents of the other dataframe
        tdf = pd.DataFrame(zip(pdf.preds[i], pdf.actual[i]), columns = ["preds", "actual"])

        # Get the predicted labels
        preds = pdf.preds[i]
        actual = pdf.actual[i]

        # Get the label indices that are part of each sequence
        l1 = len(pdf.prot1[i])
        l2 = len(pdf.prot2[i])
        l3 = len(pdf.prot3[i])

        # Split the ground truth and predicted labels to correspond to the original sequences
        # The additional 1s account for separator tokens
        lab1 = preds[1 : 1 + l1]
        lab2 = preds[1 + l1 + 1 : 1 + l1 + 1 + l2]
        lab3 = preds[1 + l1 + 1 + l2 + 1 : 1 + l1 + 1 + l2 + 1 + l3]
        act1 = actual[1 : 1 + l1]
        act2 = actual[1 + l1 + 1 : 1 + l1 + 1 + l2]
        act3 = actual[1 + l1 + 1 + l2 + 1 : 1 + l1 + 1 + l2 + 1 + l3]
        # Find the places that are a -100 (indicating separator/end of sequence),
        # and kick them out
        tdf_n100 = tdf[tdf.actual != -100]
        tp.append(len(tdf_n100[(tdf_n100.preds == 1) & (tdf_n100.actual == 1)]))
        tn.append(len(tdf_n100[(tdf_n100.preds == 0) & (tdf_n100.actual == 0)]))
        fp.append(len(tdf_n100[(tdf_n100.preds == 1) & (tdf_n100.actual == 0)]))
        fn.append(len(tdf_n100[(tdf_n100.preds == 0) & (tdf_n100.actual == 1)]))

        prot1_l.append(pdf.prot1[i])
        prot2_l.append(pdf.prot2[i])
        prot3_l.append(pdf.prot3[i])
        lab1_l.append(lab1.tolist())
        lab2_l.append(lab2.tolist())
        lab3_l.append(lab3.tolist())
        act1_l.append(act1)
        act2_l.append(act2)
        act3_l.append(act3)

    # Make dataframe
    ret_frame = pd.DataFrame(zip(prot1_l, prot2_l, prot3_l, lab1_l, lab2_l, lab3_l, 
                                 act1_l, act2_l, act3_l, tp, tn, fp, fn), 
                                 columns = [  "Prot1",   "Prot2",   "Prot3", 
                                            "Labels1", "Labels2", "Labels3", 
                                            "Actual1", "Actual2", "Actual3", 
                                            "TP", "TN", "FP", "FN"])
    return(ret_frame)

def get_test_results(dataset, model, data_collator):
    '''Sets up and then runs the dataset through the model.

    Parameters
    ----------
    dataset : Huggingface datasetdict
        Dataset containing the examples to run through the model
    model : Huggingface model
        Fine-tuned model for running the predictions
    data_collator : Huggingface datacollator
    
    Returns
    -------
    (output_list, batch_list) : tuple

    output_list : list
        A list of the outputs (predictions) from the model on the dataset
    batch_list : list
        A list of the batched inputs to the model

    '''
    dataloader_params = {
        "batch_size": 8, # Use 8 because higher tends to cause Out of Memory issues
        "collate_fn": data_collator,
        "num_workers": 4,
        "pin_memory": True,
    }

    # Get the data loaded and move to GPU
    dataloader = torch.utils.data.DataLoader(dataset, **dataloader_params)
    device = "cuda"
    model.to(device)
    model.eval()
    
    output_list, batch_list = [], []
    
    # Get the predictions
    with torch.no_grad():
        for batch in dataloader:
            outputs = model(**batch.to(device))
            output_list.append(outputs)
            batch_list.append(batch)
    return((output_list, batch_list))

def evaluate_test_results(output_batch, tokenizer):
    '''Takes the output of the model for two protein sequences and formats it nicely 
    for human (and machine-readable) output.
    
    Parameters
    ----------
    output_batch : tuple
        A tuple of the output and batch returned by the Huggingface model
    tokenizer : Huggingface tokenizer

    Returns
    -------
    pandas dataframe
        A pandas dataframe, where each row consists of processed output for each input set of sequences.
    '''
    outputs = output_batch[0]
    batch = output_batch[1]
    r_frames = []
    for i in range(len(outputs)):
        r_frames.append(process_output(outputs[i], batch[i], tokenizer))
    return(pd.concat(r_frames))
    
def evaluate_test_results_s3(output_batch, tokenizer):
    '''Takes the output of the model for three protein sequences and formats it nicely 
    for human (and machine-readable) output.
    
    Parameters
    ----------
    output_batch : tuple
        A tuple of the output and batch returned by the Huggingface model
    tokenizer : Huggingface tokenizer

    Returns
    -------
    pandas dataframe
        A pandas dataframe, where each row consists of processed output for each input set of sequences.
    ''' 
    outputs = output_batch[0]
    batch = output_batch[1]
    r_frames = []
    for i in range(len(outputs)):
        r_frames.append(process_output_s3(outputs[i], batch[i], tokenizer))

    return(pd.concat(r_frames))
    
def get_inference_results(dset, model, data_collator, tokenizer):
    '''Function to call the downstream functions for running the model in inference on
    a dataset (two interacting proteins), and then evaluating the results.

    Parameters
    ----------
    dset : Huggingface dataset object
    model : Huggingface TokenClassification model
        Fine-tuned Huggingface model
    data_collator : Huggingface data_collator object
        The data_collator for the dataset
    tokenizer : Huggingface tokenizer
        Tokenizer for tokenizing the input sequences

    Returns
    -------
    evalres : pandas dataframe
        A dataframe with columns reporting the sequences of the input proteins, the probability the model assigns to each amino acid of bel>
    '''
    res = get_test_results(dset, model, data_collator)
    evalres = evaluate_test_results(res, tokenizer)
    return(evalres)
    
def get_inference_results_s3(dset, model, data_collator, tokenizer):
    '''Function to call the downstream functions for running the model in inference on
    a dataset (three interacting proteins - heavy/light chain of antibody + antigen), 
    and then evaluating the results.

    Parameters
    ----------
    dset : Huggingface dataset object
    model : Huggingface TokenClassification model
        Fine-tuned Huggingface model
    data_collator : Huggingface data_collator object
        The data_collator for the dataset
    tokenizer : Huggingface tokenizer
        Tokenizer for tokenizing the input sequences

    Returns
    -------
    evalres : pandas dataframe
        A dataframe with columns reporting the sequences of the input proteins, the probability the model assigns to each amino acid of bel>
    '''
    res = get_test_results(dset, model, data_collator)
    evalres = evaluate_test_results_s3(res, tokenizer)
    return(evalres)

def convert_labels_to_class(lab_list, t):
    label_list = []
    for item in lab_list:
        if item < t:
            label_list.append(0)
        else:
            label_list.append(1)

    return(label_list)

def assign_output_labels(mdf, threshold):
    mdf["Classes1"] = mdf["Labels1"].apply(lambda x: convert_labels_to_class(x, threshold))
    mdf["Classes2"] = mdf["Labels2"].apply(lambda x: convert_labels_to_class(x, threshold))
    return(mdf)
    
def produce_class_output(rdf_file, idf_file, threshold):
    rdf = pd.read_csv(rdf_file, delimiter = "\t")
    ids = pd.read_csv(idf_file, delimiter = "\t")

    rdf = rdf[["Prot1", "Prot2", "Labels1", "Labels2"]]
    ids.columns = ["ids", "Prot1", "Prot2"]
    
    mdf = ids.merge(rdf, how = "inner", on = ["Prot1", "Prot2"])
    
    list_columns = ["Labels1", "Labels2"]
    
    for col in list_columns:
        mdf[col] = mdf[col].apply(ast.literal_eval)

    adf = assign_output_labels(mdf, threshold)

    return(adf)

def main(config):
    with open(config, "r") as f:
        config = json.load(f)
    
    output_dir = config["output_data_dir"]
    query_file = config["query"]
    target_file = config["target"]
    seq_type = config["params"]["seq_type"]
    model_path = config["model_path"]
    t_seed = config["params"]["t_seed"]
    pt_model = config["params"]["pt_model"]
    t = config["params"]["threshold"]

    now = datetime.datetime.now()
    ts = now.strftime("%Y_%m%d_%H%M")

    ###### Run in inference #######
    transformers.set_seed(t_seed)

    # Create the tokenizer and collator from ESM2
    # tokenizer = transformers.AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D", do_lower_case = False, is_split_into_tokens = Tru>
    tokenizer = transformers.AutoTokenizer.from_pretrained(pt_model, do_lower_case = False, is_split_into_tokens = True)
    data_collator = transformers.DataCollatorForTokenClassification(tokenizer = tokenizer)

    # Load and tokenize the dataset, filtering out examples > 1024 tokens long
    #tokenized_dset = datasets.Dataset.load_from_disk(config["datafile"])
    #tokenized_dset = tokenized_dset.filter(filter_long_examples)

    # Load fasta files
    seqs1 = read_fasta(query_file)
    seqs2 = read_fasta(target_file)

    # Create the pairings
    sdf = get_seq_combinations(seqs1, seqs2)

    dset = datasets.Dataset.from_pandas(sdf)

    tokenized_dset = dset.map(preprocess_fn, batched = True, fn_kwargs={"tokenizer": tokenizer})

    # Setup the model, using the default config from ESM2, but loading a fine-tuned model
    # mconfig = transformers.models.esm.configuration_esm.EsmConfig.from_pretrained("facebook/esm2_t6_8M_UR50D")
    mconfig = transformers.models.esm.configuration_esm.EsmConfig.from_pretrained(pt_model)
    model = transformers.EsmForTokenClassification.from_pretrained(model_path, config = mconfig)

    model.eval()

    # Ouptut file names
    ifn = os.path.join(output_dir, ts + "_inference_results.tsv")
    ids_fn = os.path.join(output_dir, ts + "_sequence_ids.txt")
    lfn = os.path.join(output_dir, ts + "_class_labels.txt")

    # Differentiate between pseudo-tripartite - antibody heavy and light chain vs antigen (3) - and bi-partite (protein vs protein)
    #
    # seq_type lives under config["params"] (App-PPI.pl nests the whole params
    # hash there); there is no top-level copy, so reading config["seq_type"]
    # here raised KeyError on every job. Use the local read at the top of main().
    #
    # str() because the JSON type is caller-dependent and nothing coerces it:
    # AppScript::preprocess_parameters only special-cases bool and enum, so an
    # int-typed param arrives as whatever the submitter sent. Observed '2' (a
    # string) from a real job, while the app spec's own default is the number 2
    # -- and a bare 3 != "3" is True in Python, which would have taken the
    # bipartite branch silently instead of crashing.
    if str(seq_type) != "3":
        inference_results = get_inference_results(   tokenized_dset.select_columns(["input_ids", "attention_mask", "labels"]), model, data_collator, tokenizer)
    else:
        inference_results = get_inference_results_s3(tokenized_dset.select_columns(["input_ids", "attention_mask", "labels"]), model, data_collator, tokenizer)

    inference_results.to_csv(ifn, index = None, sep = "\t")

    # Convert to pandas DataFrame
    df = tokenized_dset.to_pandas()

    # Select columns
    df_subset = df[["ids", "seq1", "seq2"]]
    df_subset.columns = ["ids", "Prot1", "Prot2"]

    # Write to tab-separated file
    df_subset.to_csv(ids_fn, sep='\t', index = False)

    lab_out = produce_class_output(ifn, ids_fn, t)

    lab_out.to_csv(lfn, sep = "\t", index = False)

    return(0)

if __name__ == "__main__":

    #print("hello world")
    main(sys.argv[1])
