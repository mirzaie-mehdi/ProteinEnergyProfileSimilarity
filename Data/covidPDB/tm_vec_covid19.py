# title: "TM vec"
# author: "peyman"
# date: "2024-04-09"
# ----------------------------------------
import time
import torch
from transformers import T5EncoderModel, T5Tokenizer
import re
import gc

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
tokenizer = T5Tokenizer.from_pretrained("/home/peymanc/apps/Rostlab/prot_t5_xl_uniref50", do_lower_case=False )
model = T5EncoderModel.from_pretrained("/home/peymanc/apps/Rostlab/prot_t5_xl_uniref50")
gc.collect()
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
model = model.to(device)
model = model.eval()
# ------------------------------
# ---------TM-Vec model paths
tm_vec_model_cpnt = "/home/peymanc/apps/Rostlab/prot_t5_xl_uniref50/tm_vec_cath_model.ckpt"
tm_vec_model_config = "/home/peymanc/apps/Rostlab/prot_t5_xl_uniref50/tm_vec_cath_model_params.json"
# ---------Load the TM-Vec model
tm_vec_model_config = trans_basic_block_Config.from_json(tm_vec_model_config)
model_deep = trans_basic_block.load_from_checkpoint(tm_vec_model_cpnt, config=tm_vec_model_config)
model_deep = model_deep.to(device)
model_deep = model_deep.eval()
# --------------------------------------
start = time.time()
path = '/home/peymanc/Desktop/proj/ProteinEnergyProfileSimilarity/Data/covidPDB/'
df = pd.read_csv(path+'spike_close.csv')

n=df.shape[0]
disTM = np.zeros((n,n))
for i in range(0,n,1):
  seq_1 = df.loc[i, 'seq']
  seq_1 = np.expand_dims(seq_1, axis=0)
  protrans_seq_1 = featurize_prottrans(seq_1, model, tokenizer, device).detach()
  embedded_seq_1 = embed_tm_vec(protrans_seq_1, model_deep, device)
  for j in range(i,n,1):
    print(i,j)
    seq_2 = df.loc[j, 'seq']
    seq_2 = np.expand_dims(seq_2, axis=0)
    protrans_seq_2 = featurize_prottrans(seq_2, model, tokenizer, device).detach()
    embedded_seq_2 = embed_tm_vec(protrans_seq_2, model_deep, device)
    predicted_tm_score = cosine_similarity_tm(torch.tensor(embedded_seq_1), torch.tensor(embedded_seq_2))
    disTM[i,j] = predicted_tm_score.numpy()[0]
    disTM[j,i] = predicted_tm_score.numpy()[0]
disTM = pd.DataFrame(disTM)
disTM.to_csv(path+'TM_vec.csv')
elapsed = (time.time() - start)
print(elapsed)

