# scRNA-Integ

The goal of this potential workflow is to assess performance of tools such as GraphST, alongside more standardized single cell workflows for spatial analysis. In addition, we hope to include a streamline method of cell annotation utilyzing SingleR or other alternative methods in the future. The workflow will be ran as a python script calling from a module of functions that access the tools involved. Portions of the workflow will be done using R, wrapped within python. 

## Significant Updates

| 4.23.23 | Adapted clust.py and clust.sh for utilizing GraphST clustering algorithm. 
- Refined workflow with documentation, OOP , and argparse.

| 5.08.23 | Added preprocessing functionality if input file is a seurat object. (.Rds) 
- Troubleshooting file signature errors with .Rds in preprocess_seurat() 

## Goals 

- Reformat clust.py OOP as a loadable module for later master script. 
- Add more basic spatialomic analysis functionality. Differential expression, cell type clustering. 
- Add annotation functionality to master script. 
- Create refined .yaml file for necessary imports and dependencies. 
