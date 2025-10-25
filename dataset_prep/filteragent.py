from pathlib import Path
from tqdm.notebook import tqdm
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina
from rdkit import DataStructs
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D

def _compute_props(df):
        mol = []
        rb = []
        hbd = []
        hba = []
        tpsa = []
        fsp3 = []
        heavy = []
        ch = []
        for smi in tqdm(df['smiles'], desc="Step 1. Read and pre-filtering dataset based on physical properties"):
            m = Chem.MolFromSmiles(smi)
            mol.append(m)
            rb.append(rdMolDescriptors.CalcNumRotatableBonds(m))
            hbd.append(rdMolDescriptors.CalcNumHBD(m))
            hba.append(rdMolDescriptors.CalcNumHBA(m))
            tpsa.append(rdMolDescriptors.CalcTPSA(m))
            fsp3.append(rdMolDescriptors.CalcFractionCSP3(m))
            heavy.append(m.GetNumHeavyAtoms())

        df['mol'] = mol
        df['rb'] = rb
        df['hbd'] = hbd
        df['hba'] = hba
        df['tpsa'] = tpsa
        df['fsp3'] = fsp3
        df['heavy'] = heavy
        return df

def _get_cluster_representative(df, cutoff):
   fps = df['morgan'] 
   n_fps = len(fps)
   dists = []
   for i in tqdm(range(1, n_fps), desc="Calculating similarity"):
       sims = DataStructs.BulkTanimotoSimilarity(fps[i], list(fps[:i]))
       dists.extend([1 - x for x in sims])

   print("Performing clustering")
   clusters = Butina.ClusterData(data=dists, nPts=n_fps, distThresh=cutoff, isDistData=True)

   # choose one representative per cluster: pick molecule with max sum of similarities to others in cluster (central)
   selected_idx = []
   for cl in clusters:
       if len(cl) == 1:
           selected_idx.append(cl[0])
           continue
       best = None
       best_score = -1
       for idx in cl:
           sims = [DataStructs.TanimotoSimilarity(fps[idx], fps[j]) for j in cl if j != idx]
           score = sum(sims)
           if score > best_score:
               best_score = score
               best = idx
       selected_idx.append(best)

   return df.loc[selected_idx]


class FilterAgent:
    def __init__(self, smiles_dir):
        self.smi_files = list(Path(smiles_dir).rglob("*.smi"))
        self.morgan = AllChem.GetMorganGenerator(radius=2, fpSize=512)
    
    def step1(self, reduce_to=10000):
        li = []

        for filename in self.smi_files:
            df = pd.read_csv(filename, sep='\s+')
            li.append(df)

        df = pd.concat(li, axis=0, ignore_index=True)
        df = _compute_props(df)

        print("Filtering..")
        cond = (
                (df['hbd'] <= 5) & (df['hba'] <= 10) &
                (df['tpsa'] <= 140) &
                (df['fsp3'] >= 0.12) &
                (df['heavy'] <= 50)
                )
        df_filtered = df[cond].reset_index(drop=True)

        li = []
        pool_size = reduce_to//5
        seed=100
        for i in range(0, 10, 2):
            li.append(df_filtered[(df_filtered['rb'] > i) & (df_filtered['rb'] <= i+2)].sample(n=pool_size, random_state=seed).reset_index(drop=True))
    
        self.df = pd.concat(li) # reduce to 10000 with diverse rotatable bonds
    
    def step2(self, ref_mol, not_analog_cutoff=0.7, analog_cutoff=0.5):
        print("Step 2. Remove similar ligands")

        tqdm.pandas(desc="Calculating fingerprints")
        self.df['morgan'] = self.df['mol'].progress_apply(self.morgan.GetFingerprint)

        ref_fp = self.morgan.GetFingerprint(ref_mol)
        
        print("Separate into analogs and non-analogs of co-ligand")
        self.df['sim_to_ref'] = DataStructs.BulkTanimotoSimilarity(ref_fp, list(self.df['morgan']))

        n_analogs = round(len(self.df)*0.1)
        df_analogs = self.df.sort_values('sim_to_ref', ascending=False).head(n_analogs).reset_index(drop=True)
        df_not_analog = self.df[~self.df['zinc_id'].isin(df_analogs['zinc_id'])].reset_index(drop=True)
        print(f"Number of analogs: {len(df_analogs)}, the rest: {len(df_not_analog)}")

        print("Clustering non-analogs")
        df_not_analog_cluster = _get_cluster_representative(df_not_analog, not_analog_cutoff) # molecules with simlarity >= 1-cutoff will be grouped together
        print("Clustering analogs")
        df_analog_cluster = _get_cluster_representative(df_analogs, analog_cutoff) # molecules with simlarity >= 1-cutoff will be grouped together

        self.df = pd.concat([df_analog_cluster, df_not_analog_cluster])
    
    def step3(self, ref_mol, n_ligs=250):

        pharm2d_score = []
        none_3dmols = 0
        fp_ref = Generate.Gen2DFingerprint(ref_mol, Gobbi_Pharm2D.factory)
        for mol in tqdm(self.df['mol'], desc="Step 3. Pharmacopohore screening"):
            fp_mol = Generate.Gen2DFingerprint(mol, Gobbi_Pharm2D.factory)
            sim = DataStructs.TanimotoSimilarity(fp_ref, fp_mol)
            pharm2d_score.append(sim)

        self.df['pharm2d_sim'] = pharm2d_score

        self.df = self.df.sort_values(by='pharm2d_sim', ascending=False).head(n_ligs)
        
    

# run scrub.py filtered.smi -o ligands.sdf --ph_low 7.4 --ph_high 7.4 --ff mmff94
# to get SDF files



