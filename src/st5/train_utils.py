import pandas as pd 
import torch
from transformers import T5Tokenizer, T5Model
from Bio import SeqIO
import re 
import random

def load_model(path: str = '/scratch/gpfs/vv7118/projects/vcmsa/prot_t5_xl_uniref50'):
    tokenizer = T5Tokenizer.from_pretrained(path)
    model = T5Model.from_pretrained("Rostlab/prot_t5_xl_uniref50", torch_dtype=torch.float32)
    return model, tokenizer

def load_metadata(file_path: str = '/scratch/gpfs/vv7118/projects/vcmsa/datasets/cath-domain-list.txt'):
    column_names = [
        'Domain', 'Class', 'Architecture', 'Topology', 'Homologous_superfamily',
        'S35', 'S60', 'S95', 'S100_cluster', 'S00_count', 'Domain_length', 'Resolution'
    ]
    metadata = pd.read_csv(file_path, delim_whitespace=True, comment='#', header=None, names=column_names)
    return metadata

def load_fasta_like_data(file_path = '/scratch/gpfs/vv7118/projects/vcmsa/data/ProtTucker/train66k.fasta', n: int = 'all'):
    """
    Load and parse data from a FASTA-like file using Biopython.
    
    Args:
        file_path (str): Path to the input file.
        
    Returns:
        dict: A dictionary where keys are sequence headers, and values are the sequences.
    """

    data = {}
    # Parse the file using SeqIO
    for record in SeqIO.parse(file_path, "fasta"):
        # Store the header and sequence in the dictionary
        data[record.id] = str(record.seq)
    if n != 'all':
        # select a random sample of n 
        keys = random.sample(list(data.keys()), n)
        data = {k: data[k] for k in keys}
    return data

def process_sequence(sequence):
    sequence = " ".join(list(re.sub(r"[UZOB]", "X", sequence)))
    return sequence 

def embed_sequences(sequence_dict, tokenizer, model, device='cuda:0'):
    """
    Embed all sequences in the input dictionary using a transformer model.
    
    Args:
        sequence_dict (dict): A dictionary where keys are sequence IDs and values are sequences.
        tokenizer: The tokenizer compatible with the transformer model.
        model: The transformer model used to generate embeddings.
        device (str): The device to use ('cpu' or 'cuda').
        
    Returns:
        dict: A dictionary where keys are sequence IDs and values are the embeddings.
    """
    embeddings = {}
    # Move model to the specified device
    model = model.to(device)
    
    for seq_id, sequence in sequence_dict.items():
        # Tokenize the sequence
        sequence = " ".join(list(re.sub(r"[UZOB]", "X", sequence)))
        ids = tokenizer.batch_encode_plus(
            [sequence],  # Wrap the sequence in a list for batch encoding
            add_special_tokens=True,
            padding="longest",
            return_tensors="pt"  # Return PyTorch tensors directly
        )
        
        input_ids = ids['input_ids'].to(device)
        attention_mask = ids['attention_mask'].to(device)

        # Generate embeddings
        with torch.no_grad():
            embedding_repr = model.encoder(input_ids=input_ids, attention_mask=attention_mask)
            # Assuming we're taking the last hidden state or pooling the embeddings
            embeddings[seq_id] = embedding_repr.last_hidden_state.mean(dim=1).cpu().numpy()  # Example: mean pooling

    return embeddings