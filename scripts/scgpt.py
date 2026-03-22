import torch
import numpy as np
import scanpy as sc
from pathlib import Path
import logging
from tqdm import tqdm
logging.basicConfig(level=logging.INFO)

from integration_utils import add_metadata, remove_slots, clean_categorical_column
from utils.io import read_anndata, write_zarr_linked

from tdc import tdc_hf_interface
from tdc.model_server.tokenizers.scgpt import scGPTTokenizer

input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
logging.info(f'Device available: {device}')
adata = read_anndata(
    input_file,
    X='layers/counts', 
    obs='obs',
    uns='uns',
    backed=False 
)

clean_categorical_column(adata, batch_key)

logging.info('Loading scGPT model and vocabulary...')
scgpt = tdc_hf_interface("scGPT")
model = scgpt.load()
model = model.to(device)
tokenizer = scGPTTokenizer()

logging.info('Matching genes with scGPT vocabulary...')
data_gene_names = adata.var_names.astype(str).str.upper().to_numpy()
model_vocab = tokenizer.vocab

valid_indices = []
clean_gene_ids = []

for i, gene in enumerate(data_gene_names):
    if gene in model_vocab:
        valid_indices.append(i)
        clean_gene_ids.append(gene)

valid_indices = np.array(valid_indices)
clean_gene_ids = np.array(clean_gene_ids)

if len(valid_indices) == 0:
    raise ValueError("Error: No gene overlap found between the dataset and scGPT vocabulary! Please check if your gene names are standard symbols.")

logging.info(f'Gene matching result: {len(adata.var_names)} -> {len(clean_gene_ids)} valid genes found.')

if hasattr(adata.X, "toarray"):
    X_data = adata.X.toarray()
else:
    X_data = adata.X

X_data_clean = X_data[:, valid_indices]

logging.info('Tokenizing data...')
tokenized_data = tokenizer.tokenize_cell_vectors(X_data_clean, clean_gene_ids)

logging.info('Generating Embeddings (Batch Processing)...')
model.eval()

all_embeddings = []
batch_size = 32
num_cells = len(tokenized_data)

with torch.no_grad():
    for i in tqdm(range(0, num_cells, batch_size), desc="Running Inference"):
        batch_slice = tokenized_data[i : i + batch_size]
        
        for cell_data in batch_slice:
            # cell_data[0]: gene_ids, cell_data[1]: values
            input_ids = torch.tensor(cell_data[0]).to(device).unsqueeze(0)
            input_vals = torch.tensor(cell_data[1]).to(device).unsqueeze(0)
            mask = (input_vals != 0).to(torch.bool)
            
            # Forward Pass
            try:
                output = model(input_ids, input_vals, attention_mask=mask)
                if isinstance(output, dict) and 'cell_emb' in output:
                    emb = output['cell_emb'].cpu().numpy()
                elif isinstance(output, dict) and 'last_hidden_state' in output:
                    emb = output['last_hidden_state'][:, 0, :].cpu().numpy()
                else:
                    emb = output[0].cpu().numpy() if isinstance(output, tuple) else output.cpu().numpy()

                    
                if len(emb.shape) == 1:
                    emb = emb.reshape(1, -1)
                    
                all_embeddings.append(emb)
                
            except Exception as e:
                logging.error(f"Error processing cell index {i}: {e}")
                all_embeddings.append(np.zeros((1, 512)))


final_embeddings = np.concatenate(all_embeddings, axis=0)

logging.info(f'Saving results... Shape: {final_embeddings.shape}')

adata.obsm["X_emb"] = final_embeddings
adata = remove_slots(adata=adata, output_type=params.get('output_type', 'embed'), keep_X=False)

add_metadata(adata, wildcards, params)

write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns'],
)