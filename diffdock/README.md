# Docking with DiffDock

Docking with DiffDock vs Autodock Vina was investigate for randomly selected compounds from the Zinc20 database. This performed in a 3 step procedure described in:

1. `1_protein_prep.ipynb`: Quick cleanup of PDB file for DiffDock
2. `2_redock/`: The redocking of FPRCK into the binding site. The file `2_redock/run.ipynb` details the procedure, and `2_redock/rmsd_check.ipynb` shows the code used to calculated RMSD difference.
3. `3_docking/`: The docking procedure on random molecules from Zinc20. The details of the runs are in `3_docking/1_run.md`, and the comparison of the DiffDock and Vina scores in `3_docking/2_rescore.ipynb`.

