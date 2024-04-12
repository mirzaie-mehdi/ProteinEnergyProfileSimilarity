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
pathApp = '/home/peymanc/apps/Rostlab/prot_t5_xl_uniref50'
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
start = time.time()
path = '/home/peymanc/Desktop/proj/ProteinEnergyProfileSimilarity/Data/csv/'
df = pd.read_csv(path+'Ferritin_Like_seq.csv')
#df['tm_vec'] = ''

n=df.shape[0]
emb=pd.DataFrame(np.zeros((n, 512)))
for i in range(0,n,1):
  print(i)
  seq = df.loc[i, 'seq']
  seq = np.expand_dims(seq, axis=0)
  protrans_seq = featurize_prottrans(seq, model, tokenizer, device).detach()
  embedded_seq = embed_tm_vec(protrans_seq, model_deep, device)
  emb.loc[i,:] = embedded_seq
emb.to_csv(path+'TM_vec_emb_Ferritin_Like_seq.csv',index=False)
elapsed1 = (time.time() - start)

start = time.time()
disTM = np.zeros((n,n))
for i in range(0,n,1):
  emb1 = emb.loc[i,:]
  emb1 = emb1.to_numpy()
  emb1 = emb1.reshape(1, -1)
  for j in range(i+1,n,1):
    print(i,j)
    emb2 = emb.loc[j,:]
    emb2 = emb2.to_numpy()
    emb2 = emb2.reshape(1, -1)
    predicted_tm_score = cosine_similarity_tm(torch.tensor(emb1), torch.tensor(emb2))
    disTM[i,j] = predicted_tm_score.numpy()[0]
    disTM[j,i] = predicted_tm_score.numpy()[0]
disTM = pd.DataFrame(disTM)
disTM.to_csv(path+'disTM_vec_Ferritin_Like_seq.csv',index=False)
elapsed2 = (time.time() - start)
print(elapsed1, elapsed2)

