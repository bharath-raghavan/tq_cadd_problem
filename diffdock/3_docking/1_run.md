mkdir poses

python3 inference.py --protein_ligand_csv protein_ligand.csv --config=default_inference_args.yaml \
 --no_final_step_noise --loglevel=INFO --samples_per_complex=10 --out_dir=./poses
 
cat poses/*/rank1.sdf > in.sdf
obabel in.sdf -Oligands.sdf -ph 7.4
mk_prepare_ligand.py -i ligands.sdf --multimol_outdir ligands/

mkdir vina_docking

vina --batch ligands/*.pdbqt --receptor ../../vina/1_receptor/1ppb_receptor.pdbqt \
    --config ../../vina/1_receptor/1ppb_receptor.box.txt --exhaustiveness=8 --dir vina_docking/ \
    --num_modes 1 | tee out.log
