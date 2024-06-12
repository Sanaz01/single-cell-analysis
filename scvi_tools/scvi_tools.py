#!/usr/bin/env python

import scanpy as sc    ####
import scvi
import pandas as pd
# import matplotlib.pyplot as plt
sc._settings.ScanpyConfig.n_jobs = 32
from pathlib import Path
sc._settings.ScanpyConfig.figdir = Path('.')    # instead of './figures/'
import sys
import warnings
warnings.filterwarnings('ignore')
from sklearn.cluster import KMeans



def load_dataset(filename):
    data = sc.read_10x_h5(filename)
    data.var_names_make_unique()
    data.obs_names_make_unique()
    return data

def filter_qc(data):

    sc.pp.filter_cells(data, min_counts=100)
    sc.pp.filter_genes(data, min_cells=3)

    ## QC plots for total genes, counts and percentage of mitochondrial and ribosomal genes in cells
    data.var['mt'] = data.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    data.var['rb'] = data.var_names.str.startswith(('RPS','RPL'))  # annotate the group of ribosomal proteins as 'rb'

    ## Calculate QC Metrics
    sc.pp.calculate_qc_metrics(data, qc_vars=['mt', 'rb'], percent_top=None, log1p=False, inplace=True)

    # Filter by % of ribosomal and mitochondrial genes in each cell barcode
    data = data[(data.obs['pct_counts_mt'] < 15) & (data.obs['pct_counts_rb'] < 45)]
    
    # Load any doublet info, if available
    # sc.pp.scrublet(data, batch_key="batch")

    return data

def LogNormalize_anndata(data):
    data = data.copy()
    # Normalize for library size correction. Total-count normalizing data matrix to 10,000 reads per cell.
    sc.pp.normalize_total(data, target_sum=1e4)

    # Log transform the data
    sc.pp.log1p(data)

    return data

# Not adding Cell cycle scores or regress_out function | can add cell cycle score and provide it in setup_anndata
# beacuse scvi-tools works on raw count data

def find_hvg_scvi(data):
    data = data.copy()
    sc.pp.highly_variable_genes(data, subset=True, layer="counts", flavor="seurat_v3", batch_key = 'batch')
    # sc.pl.highly_variable_genes(adata)
    # print(data.var.highly_variable.value_counts())
    data = data[:, data.var.highly_variable]
    return data


if __name__ == "__main__":

    adata_path = sys.argv[1]
    annot_path = sys.argv[2]
    tag = sys.argv[3]

    # Load AnnData
    adata = sc.read_10x_h5(adata_path)
    # Load Annotation file
    annot = pd.read_csv(annot_path, delimiter='\t')

    # print('Data loaded')
    annot = annot.loc[(annot.Chr=='chr5') & annot.Start.between(171410537, 171410545)]
    annot['bar'] = annot.apply(lambda x: f"{x.cellBarcode}-1",axis=1)
    annot = annot[['bar', 'genotype']]
    annot = annot.drop_duplicates('bar')
    annot.set_index('bar', inplace=True)

    annot_dict = {'0/0': 'WT',
                '0/1': 'Heterozygous mutation',
                '1/1': 'Homozygous mutation'}
    annot['genotype'] = annot['genotype'].map(annot_dict)


    # Merge Anndata and annotation
    adata.obs['batch'] = annot
    adata.obs['batch'] = adata.obs['batch'].fillna('Unknown')

    # print('Annotation merged')


    adata  = filter_qc(adata)
    adata.layers["counts"] = adata.X.copy()         # store raw counts in 'counts' layers
    adata = LogNormalize_anndata(adata)             # Log-Normalize the data 
    adata.raw = adata                               # Save total RNA data as raw data

    # Keep only those batch category that have atleast 50 frequency
    value_counts = adata.obs['batch'].value_counts()
    batch_values_to_keep = value_counts[value_counts >= 50].index.tolist()
    adata=adata[adata.obs['batch'].isin(batch_values_to_keep)] 

    adata=find_hvg_scvi(adata)
    adata = adata.copy()

    # run setup_anndata(), which alerts scvi-tools to the locations of various matrices inside the anndata
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        # categorical_covariate_keys=["cell_source", "donor"],
        continuous_covariate_keys=["pct_counts_mt", "pct_counts_rb"],
        batch_key="batch"
    )

    ## Training scVI model
    model = scvi.model.SCVI(adata)
    model.train(max_epochs=300)
    # model.save('scVI-model', overwrite=True)

    # Load saved model
    # model = scvi.model.SCVI.load('scVI-model', adata=adata)

    # store latent representation in obsm
    SCVI_LATENT_KEY = "X_scVI"
    adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()         # latent representation of each cell in the dataset
    
    # store the normalized values back in the anndata layers
    SCVI_NORMALIZED_KEY = "scvi_normalized"
    adata.layers[SCVI_NORMALIZED_KEY] = model.get_normalized_expression(library_size=1e4)

    # Cluster using latent values
    sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)

    # sc.tl.leiden(adata, key_added="leiden_scvi", resolution=1)
    X_pca = adata.obsm[SCVI_LATENT_KEY]
    kmeans = KMeans(n_clusters=2, random_state=0).fit(X_pca)
    adata.obs['kmeans2'] = kmeans.labels_.astype(str)

    sc.tl.umap(adata, min_dist=0.3)
    
    sc.pl.umap(adata, layer='scvi_normalized', color=['kmeans2', 'batch'], title=[str(tag), 'batch'], save=f"_{tag}_scvi_cluster.png" )
    # sc.pl.umap(adata, color=['kmeans2'], legend_loc='on data', save=f"_{tag}_scvi_cluster2.png") 
    
    # Save Anndata
    # file_name = 
    # adata.write(filename=f'adata_scvi_{tag}.h5')


    ##############################
    ## Performing DEGs Analysis ##  using KMeans Clustering
    ##############################

    adata2 = adata.raw.to_adata()
    adata2.var_names_make_unique()
    adata2.obs_names_make_unique()

    
    sc.tl.rank_genes_groups(adata2, groupby='kmeans2', groups='all', reference="rest", use_raw=False, method='wilcoxon', key_added='rank_genes', rankby_abs=True, pts=True)
    
    # sc.tl.rank_genes_groups(adata2, groupby='leiden', groups='all', reference="rest", use_raw=False, method='wilcoxon', key_added='rank_genes', rankby_abs=True, pts=True)
    ## method = 'logreg', 't-test', 'wilcoxon', 't-test_overestim_var'
    
    sc.tl.filter_rank_genes_groups(adata2, key='rank_genes', key_added='filtered_rank_genes', min_in_group_fraction=0.25, min_fold_change=1, max_out_group_fraction=0.5, compare_abs=True)
    DEgenes = sc.get.rank_genes_groups_df(adata2, group=None, key='filtered_rank_genes', pval_cutoff=1e-5)
    DEgenes = DEgenes.loc[~DEgenes.names.isnull()]
    DEgenes.to_csv(f'scvi_DEG_{tag}.csv', header=True, index=False)
