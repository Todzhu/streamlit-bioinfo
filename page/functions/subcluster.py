import scanpy as sc
import streamlit as st
from utils.utils import *
from streamlit_pills import pills


def form(adata, method, context):
    params, response = context.columns([3, 7])
    with params.form('sub-cluster'):
        cluster = st.multiselect(label = '细胞Cluster/类型：', options = sorted(adata.obs[method].drop_duplicates().to_list()), placeholder='细胞Cluster/类型')
        n_pcs = st.number_input(label='n_pcs:', min_value=5, max_value=50, value=30, step=1)
        n_neighbors = st.number_input(label='n_neighbors:', min_value=10, max_value=50, value=15, step=1)
        resolution = st.number_input(label='resolution:', min_value=0.1, max_value=1.0, value=0.3, step=0.1)
        submitted = st.form_submit_button("🚀Subset", type='primary')

    if submitted:
        sub_adata = adata[adata.obs[method].isin(cluster) , :]
        st.write(sub_adata, cluster, n_pcs, n_neighbors, resolution)

        with response.expander(f"**:red[输出结果]** ", expanded=True):
            tabs = st.tabs(["🍒 UMAP", "🍒 TSNE",  "🍒 Cluster Correlation", "🍒 Cluster statistics"])        


def subcluster():
    _, context, _ = st.columns([0.5, 8, 1.2])
    _.button(label='↩ 返回', on_click=callback, type='primary')
    with context:
        st.subheader("🍓 Subclustering scRNA-seq datasets")
        st.write("After we have completed the scRNA-seq workflow and identified the various cell types present in our samples, we might decide that for a particular cell type, we would like to identify subtypes. For example, if we have a large cluster of CD4+ Helper T cells, we may want to identify subsets of Th1, Th2, Th17, Th9, and Tfh cells. To identify these cell subsets, we would subset the dataset to the cell type(s) of interest (e.g. CD4+ Helper T cells). ")
        uploaded_file = st.file_uploader(label="h5ad文件上传:", type='h5ad')
        if uploaded_file:
            adata = sc.read_h5ad(uploaded_file)
            st.write(adata)
            columns = adata.obs.columns
            types = [i for i in columns if i=='cell_type' or i.startswith('leiden_')]
            method = pills("", types, ["🎯", "🍊"], label_visibility="collapsed")
            form(adata, method, context)
        else:
            st.info('上传.h5ad文件, 选择细胞类型或细胞Cluster继续亚群细分。', icon="ℹ️")


   

