#!/usr/bin/env python

import pandas as pd
import scanpy as sc
from pathlib import Path
sc._settings.ScanpyConfig.figdir = Path('.')    # instead of './figures/'
sc._settings.ScanpyConfig.n_jobs = 32
import warnings
warnings.filterwarnings('ignore')
import sys
# import functions as func


# print('Libraries loaded')

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

def concat_anndata(*args):
    import anndata as ad
    data = ad.concat(*args, join='outer', label='batch')
    return data

def LogNormalize_anndata(data):
    data = data.copy()
    # Normalize for library size correction. Total-count normalizing data matrix to 10,000 reads per cell.
    sc.pp.normalize_total(data, target_sum=1e4)

    # Log transform the data
    sc.pp.log1p(data)

    return data

def remove_CellCycle_effects(data):
    data = data.copy()
    ## General genes expressed in S and G2 phases of cell cycle
    s_genes = ['MCM5','PCNA','TYMS','FEN1','MCM7','MCM4','RRM1','UNG','GINS2','MCM6','CDCA7','DTL','PRIM1',
            'UHRF1','CENPU','HELLS','RFC2','POLR1B','NASP','RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7',
            'POLD3','MSH2','ATAD2','RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2',
            'USP1','CLSPN','POLA1','CHAF1B','MRPL36','E2F8']
    g2m_genes = ['HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2','NUF2','CKS1B',
                'MKI67','TMPO','CENPF','TACC3','PIMREG','SMC4','CCNB2','CKAP2L','CKAP2','AURKB','BUB1',
                'KIF11','ANP32E','TUBB4B','GTSE1','KIF20B','HJURP','CDCA3','JPT1','CDC20','TTK','CDC25C',
                'KIF2C','RANGAP1','NCAPD2','DLGAP5','CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA','PSRC1',
                'ANLN','LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA']

    s_genes = [x for x in s_genes if x in data.var_names]
    g2m_genes = [x for x in g2m_genes if x in data.var_names]

    ## Cell Cycle Score
    sc.tl.score_genes_cell_cycle(data, s_genes=s_genes, g2m_genes=g2m_genes)

    return data
    

def find_hvg(data):
    data = data.copy()
    sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='batch')
    # sc.pl.highly_variable_genes(adata)
    # print(data.var.highly_variable.value_counts())
    data = data[:, data.var.highly_variable]
    return data

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

    adata = filter_qc(adata)
    adata.layers["counts"] = adata.X.copy()                # store raw counts in 'counts' layers
    adata = LogNormalize_anndata(adata)                    # Log-Normalize the data 
    adata.raw = adata                                      # Save total RNA data as raw data

    # Keep only those batch category that have atleast 50 frequency
    value_counts = adata.obs['batch'].value_counts()
    batch_values_to_keep = value_counts[value_counts >= 50].index.tolist()
    adata=adata[adata.obs['batch'].isin(batch_values_to_keep)] 

    adata=find_hvg(adata)                                  # Find highly variable genes

    # adata = remove_CellCycle_effects(adata)         # Calculate cell-cycle scores
    regress_out_variables = ['total_counts', 'n_counts', 'pct_counts_mt', 'pct_counts_rb']
    if 'S_score' in adata.obs.columns:
        regress_out_variables = regress_out_variables + ['S_score', 'G2M_score'] 
    sc.pp.regress_out(adata, regress_out_variables)     # regressing out the confounding variables: total UMI count, MT percentage, and cell cycle scores - to eliminate their influence on downstream analysis
    sc.pp.scale(adata, max_value=10)                    # Scale each gene to zero mean and unit variance. equal weight in downstream analyses, so that highly-expressed genes do not dominate. Important before PCA

    sc.tl.pca(adata, svd_solver='auto')                    # Dimensionality reduction
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50)   # compute the neighborhood graph of cells using the PCA representation of the data matrix
    # sc.tl.leiden(adata, resolution=1)                                # Leiden graph-clustering method
    from sklearn.cluster import KMeans
    X_pca = adata.obsm['X_pca']
    kmeans = KMeans(n_clusters=2, random_state=0).fit(X_pca)
    adata.obs['kmeans2'] = kmeans.labels_.astype(str)
    sc.tl.umap(adata)                                  # Performing primary UMAP Clustering

    sc.pl.umap(adata, color=['kmeans2', 'batch'], title=[str(tag), 'batch'], save=f"_{tag}_scanpy_cluster.png")         ## Split plot by batch
    sc.pl.umap(adata, color=['kmeans2'], legend_loc='on data', save=f"_{tag}_scanpy_cluster2.png")         

    ## Save Anndata
    # file_name = f'adata_{tag}.h5'
    # adata.write(filename=file_name)
    # adata.obs.to_csv(f'{tag}_obs.csv')
 
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
    DEgenes.to_csv(f'scanpy_DEG_{tag}.csv', header=True, index=False)