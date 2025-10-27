# Virtual Screening in Autodock Vina

The procedure used to perform virtual screening of the 250 compounds selected in the dataset prepartion step is detailed here. A 4 steps procedure is performed. All are detailed in dedicated Jupyter notebooks, except for virtual screening step. Instead a markdown file is provided with commands to execute.

The molecules obtained from virtual screening are ranked based on docking score, clustered for strucural similarity and the top 10 are selected. This is detailed in `4_top10.ipynb`, and the results available in `top10.csv`.