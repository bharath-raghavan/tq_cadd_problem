import re
from rdkit import Chem

def get_first_mol_from_pdbqt(file, return_type='mol'):
    """Helper function to read PDBQT files from Autodock Vina

    This function simply reads and return data on the first molecule in
    a PDBQT file of Autodocl Vina. This is usually the most scored mol
    from a docking simulation.

    Parameters
    ----------
    file : string
        The file name of the PDBQT file.
    return_type : string
        What to return. Can take the following values:
            mol -> returns a RDKit from coordinates + score
            smiles -> returns smiles string + score
            score_only -> return only the docking score

    Returns
    -------
    (RDKitMol, int)
    or,
    (string, int)
    or,
    int
        Return depends on return_type

    Examples
    --------
    >>> get_first_mol_from_pdbqt('ligand.pdbqt', return_type="smiles")
    CC(=O)C, -6
    >>> get_first_mol_from_pdbqt('ligand.pdbqt', return_type="score_only")
    -6
    """
    
    mol_block = ""
    score = None
    smiles = ""
    with open(file, 'r') as f:
        for i, line in enumerate(f):
            if i == 1:
                finds = re.findall(r"REMARK VINA RESULT:\s+([0-9\-.]+)", line)
                if finds != []:
                    score = float(finds[0])

            if line.startswith("REMARK SMILES") and 'IDX' not in line:
                line_ = line.strip().split()
                smiles = line_[2]
                
            if line.startswith("ATOM"):
                line_ = line.strip()
                mol_block += line_[:len(line_)-2]+'\n'

            if line.strip() == 'ENDMDL':
                break
    
    if return_type == 'mol':
        return Chem.MolFromPDBBlock(mol_block, removeHs=True), score
    elif return_type == 'smiles':
        return smiles, score
    elif return_type == 'score_only':
        return score
        
        