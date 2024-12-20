# title: "TM score"
# author: "peyman"
# date: "2023-12-10"
# ----------------------------------------
# pip install tmtools
# pip install biopython
# ---------------------------------
import time
import os
import pandas as pd
import numpy as np
from tmtools import tm_align
from tmtools.io import get_structure, get_residue_data
from tmtools.testing import get_pdb_path
# --------------------------------------
path=os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(path+"/spike_seq.csv")
# ---------------------
start = time.time()
n=df.shape[0]
disTM = np.zeros((n,n))
for i in range(0,n,1):
  print(i)
  p1 = df.iloc[i,0]
  p1 =p1[0:4]
  s1 = get_structure(get_pdb_path(path+'/'+p1))
  chain1 = next(s1.get_chains())
  coords1, seq1 = get_residue_data(chain1)
  for j in range(i+1,n,1):
    p2 = df.iloc[j,0]
    p2 = p2[0:4]
    s2 = get_structure(get_pdb_path(path+'/'+p2))
    chain2 = next(s2.get_chains())
    coords2, seq2 = get_residue_data(chain2)
    res = tm_align(coords1, coords2, seq1, seq2)
    disTM[i,j] = res.tm_norm_chain1
    disTM[j,i] = res.tm_norm_chain2
disTM = pd.DataFrame(disTM)
disTM.to_csv(path+'/disTM_score_spike.csv',index=False)
elapsed = (time.time() - start)
print(elapsed)

