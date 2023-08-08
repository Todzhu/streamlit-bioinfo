import os
import math
import diopy
import shutil
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import cosg as cosg
from PIL import Image
import seaborn as sns
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

import warnings
warnings.filterwarnings("ignore")

sc.settings.verbosity = 1
sc.logging.print_header()
sc.settings.set_figure_params(dpi=100, dpi_save=300, figsize=(5, 5), facecolor='white')

plt.style.use('seaborn-white')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=12)
plt.rcParams['pdf.fonttype'] = 42


@st.cache_data
def make_report_dir(out_dir):   
    os.makedirs(out_dir, exist_ok=True)
    return out_dir

def image_view(_container, png_file, caption=''):
    image = Image.open(png_file)
    _container.image(image, caption=caption, use_column_width=True)

def sub_cluster(_adata, n_pcs=10, n_neighbors=15, resolution=0.3, organism=None, out_dir='./'):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    sc.tl.pca(_adata, svd_solver='arpack')
    sc.pl.pca(_adata, color='sample', show=False)
    plt.savefig(os.path.join(out_dir, 'pca_sample.png'), facecolor='white', bbox_inches="tight", dpi=300)
    # sc.pl.pca(_adata, color='group', show=False)
    # plt.savefig(os.path.join(out_dir, 'pca_group.png'), facecolor='white', bbox_inches="tight", dpi=300)

    sc.pp.neighbors(_adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    # UMAP
    sc.tl.umap(_adata)
    sc.tl.leiden(_adata, resolution=resolution)
    sc.pl.umap(_adata, color='leiden', size=25, show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_cluster.png'), facecolor='white', bbox_inches="tight", dpi=300)

    # TSNE
    sc.tl.tsne(_adata)
    sc.pl.tsne(_adata, color='leiden', size=25, show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_tsne_cluster.png'), facecolor='white', bbox_inches="tight", dpi=300)

    # Cluster correlation
    sc.pl.correlation_matrix(_adata, 'leiden', show=False)
    plt.savefig(os.path.join(out_dir, 'cluster_correlation.png'), facecolor='white', bbox_inches="tight", dpi=300)

    # Marker
    cosg.cosg(_adata,
            groups='all',
            key_added='cosg',
            mu=1,
            n_genes_user=500,
            groupby='leiden')           
    colnames = ['names', 'scores']
    test = [pd.DataFrame(_adata.uns["cosg"][c]) for c in colnames]
    test = pd.concat(test, axis=1, names=[None, 'group'], keys=colnames)
    markers = {}
    cats = _adata.obs.leiden.cat.categories
    for i, c in enumerate(cats):
        cell_type_df = test.loc[:, 'names'][c]
        scores_df = test.loc[:, 'scores'][c]
        markers[c] = cell_type_df.values.tolist()[:4]
        all_markers = pd.DataFrame({'names': list(cell_type_df), 'scores': list(scores_df)})
        all_markers['cluster'] = c
        all_markers.to_csv(os.path.join(out_dir, f'cluster{c}.xls'), sep='\t', index=None)
    sc.pl.dotplot(_adata,
                var_names = markers, 
                groupby='leiden',
                cmap='Spectral_r',
                standard_scale='var',
                show=False)
    plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), facecolor='white', bbox_inches="tight", dpi=300)
    _adata.obs.to_csv(os.path.join(out_dir, 'all_cells_features.xls'), sep='\t')

    from matplotlib.pyplot import rc_context
    sc.set_figure_params(dpi=80, color_map = 'viridis_r')

    samples = _adata.obs['sample'].drop_duplicates().tolist()
    if len(samples) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
        fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 6))
        for i in range(len(samples)):
            if i == len(samples) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(_adata[_adata.obs["sample"]==samples[i]], color="leiden", size=25, legend_loc=legend, title=samples[i], ax=axs[i], show=False)
        plt.savefig(os.path.join(out_dir, 'umap_samples.png'), facecolor='white', bbox_inches="tight", dpi=300)

    pd.DataFrame(_adata.obsm['X_tsne'], index=_adata.obs_names, columns=['tSNE_1', 'tSNE_2']).to_csv(os.path.join(out_dir, 'data_tsne.xls'), sep='\t')
    pd.DataFrame(_adata.obsm['X_umap'], index=_adata.obs_names, columns=['UMAP_1', 'UMAP_2']).to_csv(os.path.join(out_dir, 'data_umap.xls'), sep='\t')

    if organism == 'mouse':
        _adata.var['gene'] = _adata.var.index
        _adata.var.index = [i.upper() for i in _adata.var.index]
    markers = dc.get_resource('PanglaoDB')
    # markers.to_csv('./source/PanglaoDB.xls', sep='\t', index=None)
    markers = markers[(markers[organism]=='True')&(markers['canonical_marker']=='True')&(markers['organ']!='Pancreas')]
    markers = markers[~markers.duplicated(['cell_type', 'genesymbol'])]
    dc.run_ora(mat=_adata, net=markers, source='cell_type', target='genesymbol', min_n=3, verbose=True,use_raw=False)
    acts = dc.get_acts(_adata, obsm_key='ora_estimate')
    mean_enr = dc.summarize_acts(acts, groupby="leiden", min_std=1)
    sns.clustermap(mean_enr, xticklabels=mean_enr.columns, cmap='viridis')
    annotation_dict = dc.assign_groups(mean_enr)
    _adata.obs['cell_type'] = [annotation_dict[clust] for clust in _adata.obs["leiden"]]

    if organism == 'mouse':
        _adata.var.index = _adata.var['gene']

    sc.pl.umap(_adata, color='cell_type', size=25, show=False, title='')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    sc.pl.tsne(_adata, color='cell_type', size=25, show=False, title='')
    plt.savefig(os.path.join(out_dir, 'dimplot_tsne_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)

    with rc_context({'figure.figsize': (6, 6)}):
        sc.pl.umap(_adata, color=["leiden","cell_type"], wspace=0.2, show=False)
        plt.savefig(os.path.join(out_dir, 'umap_leiden_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    
    with rc_context({'figure.figsize': (6, 6)}):
        sc.pl.tsne(_adata, color=["leiden","cell_type"], wspace=0.2, show=False)
        plt.savefig(os.path.join(out_dir, 'tsne_leiden_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)

    sc.pl.umap(_adata, color=["n_genes_by_counts", "total_counts"], wspace=0.2, show=False)
    plt.savefig(os.path.join(out_dir, 'umap_n_genes_by_counts.png'), facecolor='white', bbox_inches="tight", dpi=300)

    return _adata


@st.cache_data
def visua_gene_express(_adata, genes, out_dir):
    st.set_option('deprecation.showPyplotGlobalUse', False)
    sc.set_figure_params(dpi=100, color_map = 'viridis_r')
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.umap(_adata, color=genes, s=30, ncols=3, vmax='p99')
 
        return plt.show()







