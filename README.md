# scRNA-Integ
Creating an scRNA-seq spatial workflow which integrates cell annotation, taking advantage of cuda. 

The goal of this potential workflow is to assess performance of tools such as GraphST, alongside more standardized single cell workflows for spatial analysis. In addition, we hope to include a streamline method of cell annotation utilyzing SingleR. The workflow will be ran as a python script calling from a module of functions that access the tools involved. Portions of the workflow will be done using R, wrapped within python. 
