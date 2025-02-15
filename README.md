# ProteinEnergyProfileSimilarity
This repository contains all the necessary code and resource data to calculate the energy profile from protein structures or sequences. It includes the following functions and data:

1. A non-redundant dataset of protein chains was employed to derive the distance-dependent knowledge-based potential. These protein structures served as the basis for training and calculating the potential function. The PDB IDs used are listed in the `Train_Energy/list_with_ChainID.txt` file and the resulting distance-dependent potential between pairs of atoms is located in the `Train_Energy/energy_dell_dunbrack.mat` and `Data/csv/energy.csv` files.

2. The predictor matrix to estimate the profile of energy based on protein sequence is located at `csv/Pij.csv`

3. All functions needed to convert PDB files/Protein sequence into energy profiles are found in `Functions/Functions.R`. The `Profile_Energy.Rmd` file contains the code chunks used to calculate energy profiles for each dataset discussed in the paper. Each chunk is briefly described before its corresponding section, and the output from each chunk is already saved as an RDS file in the `Data/rds` folder.

    #### Note:
   To obtain the results of TM-Vec, you need to install the TM-Vec software as instructed at https://github.com/tymor22/tm-vec. The `Functions` folder contains four Python scripts that generate the TM-Vec results. For instance, you can reproduce the results for bacteriocins by running the script `tm_vec_bac.py`. The results of `tm_vec_bac.py`, `tm_vec_CT_Ho_cathID_filtered.py`, `tm_vec_Ferritin_Like_seq.py` scripts are already saved  as `disTM_vec_bac.csv`, `disTM_vec_CT_Ho_cathID_filtered.csv`, `disTM_vec_Ferritin_Like_seq.csv` in the `Data/csv/` folder. The result of `tm_vec_fiveSF.py` is shared via the google drive link:
   https://drive.google.com/drive/folders/1osbntpxhGUFtQsG-beiJDQvAfWeCBMzl?usp=drive_link
   

   To obtain the TM-Vec and TM-score results on the Covid dataset, you need to run the scripts `tm_vec_spike.py` and `tm_score_spike.py` located in the `Data/covidPDB/` folder. The results of these scripts are already saved in this folder as `disTM_vec_spike.csv` and `disTM_score_spike.csv`. The RMSD results for this dataset can be obtained from the script `Functions/RMSD_spike.R`. To run this script, you need to have the `muscle3.8.1551` software installed. The results of `Functions/RMSD_spike.R` script has been already saved at `Data/covidPDB/RMSD_spike.rds`.

5. The file `Results_Figures_Tables.Rmd` contains all chunks used to generate results, figures and tables for every dataset referenced in the manuscript. All necessary files for running the `Results_Figures_Tables.Rmd` file have already been created and saved using the `Profile_Energy.Rmd` file.

 ## Datasets Used in This Project
Below is the list of datasets utilized in this project:

- Training Set for Knowledge-Based Potential: [PDBIDs](Train_Energy/list_with_chainID_rm_Olaps.txt).
- Protein Domains in [ASTRAL40](Data/csv/astral-scopedom-seqres-gd-sel-gs-bib-40-2.08.fa) and [ASTRAL95](Data/csv/astral-scopedom-seqres-gd-sel-gs-bib-95-2.08.fa).
- The list of [Bacteriocin Proteins family](Data/csv/Bacteriocin.csv) available in the BAGEL database that includes 690 proteins.
- PDBIDs of [Ferritin Superfamily](Data/csv/Ferritin_Like_seq.csv) (SCOP ID: a.25.1) .
- The list of [C_terminal and Homing endonucleases](Data/csv/CT_Ho_cathID.csv) (CATH Code: 1.10.8.10 and 3.10.28.10).
- The list of protein domains of [five superfamily](Data/csv/fiveSF.csv) winged helix(SCOP ID: a.4.5), PH domain-like(SCOP ID: a.55.1), NTF-like(SCOP ID: d.17.4), Ubiquitin-like(SCOP ID: d.15.1), and Immunoglobulins(SCOP ID: b.1.1).
- Covid19 spike proteins data set [covidPDBIDs](Data/covidPDB/).
- The list of [Drug-Targets](Data/csv/41467_2019_9186_MOESM4_ESM.xlsx) including 65 antihypertensive drugs and their protein targets IDs.
- The list of 21 mammalian hemoglobins proteins in [Globin](Data/Globin/Globin.csv) family.
- Large-Scale SARS-CoV-2 Proteome Analysis across 28 families. The protein modesl can be find in the [Large_Scale_SARS2](Data/Large_Scale_SARS2) folder according to the sars_proteom column [SARS_Proteom](Data/Large_Scale_SARS2/sars_proteom.csv).

 ## Example:
The `Example/SPE.R` and `Example/CPE.R` scripts are developed to compute SPE and CPE for Alpha-globin and Beta-globin proteins across different species. Additionally, they demonstrate the UMAP of SPE and CPE representaions, as outlined below.

<img src="https://github.com/user-attachments/assets/ac22e365-eff8-4bbb-9456-a3dcd63fac65" width="400" height="300">
<img src="https://github.com/user-attachments/assets/019e55bd-d082-4a32-86d4-ee01c701899e" width="400" height="300">

Code and data sets can be found on Zenodo: https://zenodo.org/records/14765519
