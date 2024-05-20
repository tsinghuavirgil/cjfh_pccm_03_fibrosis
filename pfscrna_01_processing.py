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
umap_adata = sc.read('ipfcopd_rawadata.h5ad')
##################
#processing_adata#
##################
sc.pl.highest_expr_genes(umap_adata, n_top=20, )

sc.pp.filter_cells(umap_adata, min_genes=200)

sc.pp.filter_genes(umap_adata, min_cells=3)

umap_adata.var['mt'] = umap_adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(umap_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(umap_adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

sc.pl.scatter(umap_adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(umap_adata, x='total_counts', y='n_genes_by_counts')
###########
#pca_adata#
###########
sc.pp.normalize_total(umap_adata, target_sum=1e4)
umap_adata.raw = umap_adata
umap_adata = sc.pp.filter_genes_dispersion(umap_adata, subset = False, min_disp=.5, max_disp=None,
                                              min_mean=.0125, max_mean=10, n_bins=20, n_top_genes=None,
                                              log=True, copy=True)

sc.pp.log1p(umap_adata)

sc.pp.scale(umap_adata, max_value=10, zero_center=False)

sc.tl.pca(umap_adata, use_highly_variable=True)

sc.pl.pca_scatter(umap_adata, color=['project'])
sc.pl.pca_scatter(umap_adata, color=['diagnosis'])
sc.pl.pca_scatter(umap_adata, color=['group'])

sc.pl.pca_variance_ratio(umap_adata, log = True)

sc.pl.pca_loadings(umap_adata)
############
#umap_adata#
############
neighbor = 30
resolution = 1
pcs = 15
sce.pp.bbknn(umap_adata, 'sample_id')
sc.tl.leiden(umap_adata, resolution = resolution)#
sc.tl.umap(umap_adata)
sc.pl.umap(umap_adata, color=['leiden'])
#######################
#rank_gene_group_adata#
#######################
sc.tl.rank_genes_groups(umap_adata, 'leiden', method='wilcoxon')
deg_adata = pd.DataFrame(umap_adata.uns['rank_genes_groups']['names']).head(50)
deg_adata
######
#Save#
######
umap_adata.write('ipfcopd_umapadata.h5ad')