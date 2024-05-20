import numpy as np
import pandas as pd
import h5py
import os
import argparse
import scipy.sparse as sparse
#import pygwalker as pyg
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager as fm, rcParams
import csv
import scanpy as sc
import scanpy.external as sce
import anndata
import bbknn
#import scanorama
#import scrublet as scr
import scvelo as scv
#import cellrank as cr
import loompy as lp
from scipy.stats import pearsonr
#import palantir
#import doubletdetection
from igraph import *
#from MulticoreTSNE import MulticoreTSNE as TSNE
from anndata import read_h5ad
from anndata import read_csv
from matplotlib import rcParams
sc.settings.verbosity = 3
#cr.settings.verbosity = 2
sc.logging.print_versions()
%matplotlib inline
sc.settings.set_figure_params()#300
###################################
umap_adata = sc.read('ipfcopd_umapadata.h5ad')
#############################
#manual_cell_type_annotation#
#############################
sc.pl.umap(umap_adata, color=['EPCAM','KRT8',#Epithelial cells
                              'PECAM1','CDH5',#Endothelial cells
                              'HHIP','ASPN',#Mesenchymal cells
                              'LUM','DCN',#Fibroblast
                              'POSTN','CTHRC1',#Myofibroblast
                              'TAGLN','MYH11',#SMC
                              'PDGFRB','COX4I2',#Pericyte
                              'PTPRC',#Immune cells
                              'S100A8','C1QB',#Myeloid cells
                              'CD3D','CD3E','CD4','CD8A','CD8B',# T cells
                              'CD79A','CD79B'#B cells
                             ], cmap = 'Reds', vmax = 10)#
umap_adata.obs['major_cluster'] = umap_adata.obs['leiden'].map(
    {
        '0': 'Myeloid cells',
        '1': 'Myeloid cells',
        '2': 'Epithelial cells',
        '3': 'T & NK cells',
        ...
    }
)
###############
#marker_visual#
###############
ax = sc.pl.heatmap(umap_adata, celltype_marker, groupby='celltype', cmap='cividis',
                   dendrogram=False, vmax = 1, figsize=(2,5),
                   save = '_celltype_marker.pdf'
                  )#
###############
#density_adata#
###############
sc.tl.embedding_density(umap_adata, basis='umap', groupby='diagnosis')
#'ds_01_don', 'ds_02_sscild', 'ds_03_ipf'
sc.pl.embedding_density(umap_adata, basis='umap',
                        key='umap_density_diagnosis',
                        group='ds_01_don',colorbar_loc=None,
                        save = 'density_01_donor.pdf')
#############
#dotplot_ecm#
#############
umap_adata.layers['scaled'] = sc.pp.scale(umap_adata, copy=True).X
collagen_dict = {'Collagen': ['COL6A6','COL13A1','COL6A1','COL6A2','COL1A2','COL6A3','COL3A1','COL1A1',
                              'COL5A2','COL8A1','COL5A1','COL14A1','COL16A1','COL15A1','COL10A1','COL24A1',
                              'COL7A1','COL27A1','COL4A1','COL4A2','COL18A1','COL12A1','COL5A3','COL21A1'
                             ]
                }
x = sc.pl.dotplot(umap_adata, collagen_dict, groupby='celltype',layer='scaled',
                  var_group_rotation=0, dendrogram=False, dot_max=1,
                  dot_min=0.01, smallest_dot=10, standard_scale='var', cmap='RdGy_r',
                  save = 'collagene_adata.pdf'
                  )
######
#Save#
######
umap_adata.write('ipfcopd_umapadata.h5ad')