# Docking with DiffDock

First a directory for DiffDock to write the poses is created:

```
mkdir poses
```

DiffDock can be run with:

```
python3 inference.py --protein_ligand_csv protein_ligand.csv --config=default_inference_args.yaml \
 --no_final_step_noise --loglevel=INFO --samples_per_complex=10 --out_dir=./poses
```

The most highly ranked poses for each molecule is written to `poses/<complex_name>/rank1.sdf`. We can use these to redock them in Vina. The molecules are combined into a single SDF file:

```
cat poses/*/rank1.sdf > in.sdf
```

DiffDock removes all hydrogens from the input file, as they are implicity treated in the model. We will need to add them back for docking with Vina. We use `obabel` so as to not generate additional conformers:

```
obabel in.sdf -Oligands.sdf -ph 7.4
```

And prepare them with Meeko:

```
mk_prepare_ligand.py -i ligands.sdf --multimol_outdir ligands/
```

Finally, a directory to store the Vina docking data is created, and docking is run:

```
mkdir vina_docking

vina --batch ligands/*.pdbqt --receptor ../../vina/1_receptor/1ppb_receptor.pdbqt \
    --config ../../vina/1_receptor/1ppb_receptor.box.txt --exhaustiveness=8 --dir vina_docking/ \
    --num_modes 1 | tee out.log
```

The results are analyzed in `2_rescore.ipynb`.
