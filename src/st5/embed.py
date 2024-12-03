from datasets import SequenceDataset
import sys
import pickle 
from train_utils import load_model, load_fasta_like_data, process_sequence

model, tokenizer = load_model()
model = model.to('cuda:0')
file_names = ['CBS.vie.20seqs', "Glyco Hydro 18 D1 20 Seqs", "AAA vie 20 seqs", 'RuBisCO_large_N.vie.20seqs']

import torch
print('loading model')
state_dict = torch.load('/scratch/gpfs/vv7118/projects/vcmsa/models/2024-12-02_16-40-27/best_model.pt', map_location='cuda:0', weights_only=True)
print('loading state dict')
model.load_state_dict(state_dict) 
print('model loaded')
model.eval()

for file in file_names:
    data = load_fasta_like_data(f'/scratch/gpfs/vv7118/projects/vcmsa/data/{file}.fasta')
    n_last_layers = [1, 16]
    for n_last in n_last_layers:
        embs = []
        for name, val in data.items():
            print(name)
            print(n_last)
            seq1 = process_sequence(val)
            inputs = tokenizer(seq1, return_tensors='pt', padding=True, truncation=True)
            inputs = {key: val.to('cuda') for key, val in inputs.items()}
            outputs = model.encoder(input_ids=inputs['input_ids'], 
                                    attention_mask=inputs['attention_mask'], 
                                    output_hidden_states=True)
            # Extract the last 16 layers' hidden states
            hidden_states = outputs.hidden_states[-n_last:]
            # Concatenate the embeddings along the last dimension
            concatenated_embs = torch.cat(hidden_states, dim=-1)
            # Take the mean across the sequence length dimension
            mean_embs = concatenated_embs.mean(dim=1)
            print(mean_embs.shape)
            embs.append({name: mean_embs.cpu().detach().numpy()})


        with open(f'data/{file}_{n_last}.pkl', 'wb') as f:
            pickle.dump(embs, f)