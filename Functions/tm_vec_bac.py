# title: "TM vec"
# author: "peyman"
# date: "2024-04-09"
# ----------------------------------------
import time
import torch
from transformers import T5EncoderModel, T5Tokenizer
import re
import gc
import os

import numpy as np
import pandas as pd

import torch
from torch.utils.data import Dataset

from tm_vec.embed_structure_model import trans_basic_block, trans_basic_block_Config
from tm_vec.tm_vec_utils import featurize_prottrans, embed_tm_vec, cosine_similarity_tm


import matplotlib.pyplot as plt
import seaborn as sns 
from tmtools.testing import get_pdb_path
# --------------------------------------
pathApp=input('Enter the installation path of the `tm-vec` (ex: /path/to/prot_t5_xl_uniref50): ')
#pathApp = '/home/peymanc/apps/Rostlab/prot_t5_xl_uniref50'
tokenizer = T5Tokenizer.from_pretrained(pathApp, do_lower_case=False )
model = T5EncoderModel.from_pretrained(pathApp)
gc.collect()
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
model = model.to(device)
model = model.eval()
# ------------------------------
# ---------TM-Vec model paths
tm_vec_model_cpnt = pathApp+"/tm_vec_cath_model.ckpt"
tm_vec_model_config = pathApp+"/tm_vec_cath_model_params.json"
# ---------Load the TM-Vec model
tm_vec_model_config = trans_basic_block_Config.from_json(tm_vec_model_config)
model_deep = trans_basic_block.load_from_checkpoint(tm_vec_model_cpnt, config=tm_vec_model_config)
model_deep = model_deep.to(device)
model_deep = model_deep.eval()
# --------------------------------------
path=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
path = path+'/Data/csv/'
# ----------------------------------------------
df = pd.read_csv(path+'Bacteriocin.csv')
# ----------------------------------------------
start = time.time()
n=df.shape[0]
emb=pd.DataFrame(np.zeros((n, 512)))
for i in range(0,n,1):
  print(i)
  seq = df.loc[i, 'Seq']
  seq = np.expand_dims(seq, axis=0)
  protrans_seq = featurize_prottrans(seq, model, tokenizer, device).detach()
  embedded_seq = embed_tm_vec(protrans_seq, model_deep, device)
  emb.loc[i,:] = embedded_seq
#emb.to_csv(path+'TM_vec_emb_bac.csv',index=False)

# Calculate cosine similarity for each pair of rows
tensor_df = torch.tensor(emb.values)
cosine_similarity = torch.zeros(len(emb), len(emb))
for i in range(len(emb)):
    print(i)
    output_seq1_tensor = tensor_df[i].unsqueeze(0).repeat(len(emb), 1)
    output_seq2_tensor = tensor_df
    cosine_similarity[i] = cosine_similarity_tm(output_seq1_tensor, output_seq2_tensor)

# Convert result to DataFrame
cosine_similarity_df = pd.DataFrame(cosine_similarity.numpy(), index=emb.index, columns=emb.index)

cosine_similarity_df.to_csv(path+'disTM_vec_bac2.csv',index=False)
elapsed = (time.time() - start)

print(elapsed)


