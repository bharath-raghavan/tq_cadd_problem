# Redocking

DiffDock requires the input to be specified in a csv file. For the case of docking between 1PPB and the coligand, I use `protein_ligand.csv`. With that, DiffDock can be run:

```
python3 inference.py --protein_ligand_csv protein_ligand.csv --config=default_inference_args.yaml \
 --no_final_step_noise --loglevel=INFO --samples_per_complex=10 --out_dir=./
```

Note that `default_inference_args.yaml` is provided in the DiffDock distribution.
