# Import base python packages
import os
import torch
import pandas as pd 
from sklearn import metrics
import multiprocessing as mp 
from GraphST import GraphST
import scanpy as sc
import matplotlib.pyplot as plt
import argparse
import anndata

# Import r packages here
#from rpy2.robjects.packages import importr
#from rpy2.robjects import zellkonverter

# Argparse for easy plug-and-play 
parser = argparse.ArgumentParser(description = "Designate R and File Location, and perform GraphST Clustering")
parser.add_argument('-r','--r_dir', required = True, help = "R Directory, cannot be a docker container.")
parser.add_argument('-d','--data', required = True, help = "Data File Directory")
parser.add_argument('-m','--mode', required = True, help = "Processing mode, either cuda or cpu. Cuda recommended.")
parser.add_argument('-t','--type', required= True, help= "Data type")
parser.add_argument('-n','--clust', required= True, help= "Number of clusters")
parser.add_argument('-p','--prep', required = True, help = "Normalize/Spatialize seurat object? (y/n)")
args = parser.parse_args()

# Convert -d file from seurat to annotated h5ad. This workflow assumes your input file is a seurat object. 
os.environ['R_HOME'] = args.r_dir


# Runs through GraphST stereoseq vignette given input settings/files. 
class Clusterplot:
    'Produces spatial cluster plot via GraphST' 
    def __init__(self, rdir, n, type, data, mode):
        'Load in argparse parameters'
        self.rdir = args.r_dir
        self.n = args.clust
        self.type = args.type
        self.data = args.data
        self.mode = args.mode

    # Populate list of parameters
    def setparams(self): 
        'Set parameters and relevant directories for GraphST'
        parameters = {}
        device = torch.device(self.mode)
        os.environ['R_HOME'] = self.rdir
        dataset = f'{self.data}'
        n_clusters = self.n
        file_path = self.data
        datatype = self.type

        parameters['device'] = device
        parameters['R'] = os.environ['R_HOME']
        parameters['data'] = dataset
        parameters['clusters'] = n_clusters
        parameters['data_path'] = file_path
        parameters['datatype'] = datatype

        return parameters
    # Cluster via mclust 
    def run(parameters): 
        'Runs GraphST given parameters'
        adata = sc.read_h5ad(parameters['data_path'])
        adata.var_names_make_unique() 

        model = GraphST.GraphST(adata, parameters['datatype'], parameters['device'])
        adata = model.train()

        from GraphST.utils import clustering
        tool = 'mclust'
        clustering(adata, parameters['clusters'], tool)

        return adata
    # Produce spatial cluster plot. Adjust plot_color[] based on # of clusters. 
    def visualize(adata, parameters):
        plot_color = ["#F56867","#556B2F","#C798EE","#59BE86","#006400","#8470FF","#CD69C9",
                      "#EE7621","#B22222","#FFD700","#CD5555","#DB4C6C","#8B658B","#1E90FF",
                      "#AF5F3C","#CAFF70","#F9BD3F","#DAB370","#877F6C","#268785", '#82EF2D', '#B4EEB4']
        adata.obsm['spatial'][:, 1] = -1*adata.obsm['spatial'][:, 1]
        plt.rcParams["figure.figsize"] = (3,4)

        ax = sc.pl.embedding(adata, basis = "spatial",
                             color = "domain",
                             s = 30,
                             show = False,
                             palette = plot_color,
                             title = 'GraphST')
        ax.axis('off')
        ax.set_title(parameters['data_path'])

        # Figure will be saved as same name as file with .png extension. 
        plt.savefig(f'{parameters["data_path"]}.png', dpi = 300, bbox_inches = 'tight')
        return adata
    
class Preprocessing:
    def preprocess_seurat(parameters):
        "If working with a seurat object, this will perform batch correction, log transformation, and convert the file to h5ad"
        adata = anndata.AnnData()
        seurat = sc.read_h5ad(parameters['data'])
        
        # Extract the expression matrix, meta data, and features. 
        adata.X = seurat.X
        adata.obs = seurat.obs
        adata.var = seurat.var

        # Batch correction via combat 
        adata = sc.pp.batch_correction(adata, batch_key = "batch")

        # Log transform 
        adata = sc.pp.log1p(adata)

        # Write to h5ad file
        h5ad = adata.write_h5ad(f"{parameters['data_path']}.h5ad")

        return h5ad
    
############################################################################
############################ Running the Script ############################       
############################################################################

# Load the argparse input as Parameters vector.
Input = Clusterplot(rdir = args.r_dir, data = args.data, mode = args.mode, type = args.type, n = args.clust)
Parameters = Clusterplot.setparams(Input)

# If working with a seurat input file, set prep argument to "y" to preprocess and convert to h5ad for GraphST input.
if args.prep == "y":
    Preprocessing.preprocess_seurat(parameters=Parameters)
    Parameters['data'] = f"{Parameters['data_path']}.h5ad"
    Model = Clusterplot.run(parameters = Parameters)
    Clusterplot.visualize(adata = Model, parameters= Parameters)

# If working with a .h5ad input already, then continue to GraphST as normal. 
else:
    Model = Clusterplot.run(parameters = Parameters)
    Clusterplot.visualize(adata = Model, parameters= Parameters)

        


