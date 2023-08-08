import os
import random
import datetime
import scanpy as sc
import streamlit as st
from utils.utils import *
from streamlit_pills import pills
from utils.module import make_report_dir, sub_cluster, image_view, visua_gene_express

out_dir = './tmp'
if 'out_dir' not in st.session_state:
    st.session_state['out_dir'] = make_report_dir(os.path.join(out_dir, 'sub-Cluster_'+''.join(random.sample('0123456789zyxwvutsrqponmlkjihgfedcba', 8))))
if 'n_pcs' not in st.session_state:
    st.session_state['n_pcs'] = 10
if 'n_neighbors' not in st.session_state:
    st.session_state['n_neighbors'] = 15
if 'resolution' not in st.session_state:
    st.session_state['resolution'] = 0.3
if 'sub_adata' not in st.session_state:
    st.session_state['sub_adata'] = []


def form(adata, method, context):
    params, response = context.columns([3, 7])
    with params.form('sub-cluster'):
        cluster = st.multiselect(label = 'cell cluster/cell cype:', options = sorted(adata.obs[method].drop_duplicates().to_list()), placeholder='细胞Cluster/类型')
        n_pcs = st.number_input(label='n_pcs:', min_value=5, max_value=50, value=10, step=1)
        n_neighbors = st.number_input(label='n_neighbors:', min_value=10, max_value=50, value=15, step=1)
        resolution = st.number_input(label='resolution:', min_value=0.1, max_value=1.0, value=0.3, step=0.1)
        organism = st.selectbox(label='organism:', options=['选择物种', 'human', 'mouse'], index=0)
        submitted = st.form_submit_button("🚀Subset", type='primary')
    if submitted:
        if not cluster:
            response.error('请至少选择1个细胞Cluster/类型！', icon="🚨")
        else:
            if organism == '选择物种':
                response.error('请选择物种信息！', icon="🚨")
            else:
                sub_adata = adata[adata.obs[method].isin(cluster) , :]
                sub_adata = sub_adata.raw.to_adata()
                sub_adata.raw = sub_adata
                with st.spinner('数据加载中...'):
                    st.session_state['sub_adata'] = sub_cluster(sub_adata, n_pcs, n_neighbors, resolution, organism, st.session_state['out_dir'])

    if st.session_state['sub_adata']:
        with response.expander(f"**:red[聚类结果]** ", expanded=True):
            tabs = st.tabs(["🍒 UMAP", "🍒 TSNE", "🍒 Cluster Correlation", "🍒 Cluster Markers"])        
            image_view(tabs[0], png_file=os.path.join(st.session_state['out_dir'], 'umap_leiden_celltype.png'))
            image_view(tabs[1], png_file=os.path.join(st.session_state['out_dir'], 'tsne_leiden_celltype.png'))
            image_view(tabs[2], png_file=os.path.join(st.session_state['out_dir'], 'cluster_correlation.png'))
            image_view(tabs[3], png_file=os.path.join(st.session_state['out_dir'], 'dotplot_marker.png'))

        with params.form('gene_exp'):
            genes = st.multiselect(label='visualization of gene expression:', options=st.session_state['sub_adata'].var_names)       
            submitted_view = st.form_submit_button("🚀Plot", type='primary')

        if submitted_view:
            with response.expander(f"**:red[基因表达]** ", expanded=True):
                tabs = st.tabs(["🍓 Color map"])
                if not genes:
                    st.error('Please select genes to visualization of gene expression.', icon="🔥")
                else:
                    fig = visua_gene_express(st.session_state['sub_adata'], genes, st.session_state['out_dir'])
                    tabs[0].pyplot(fig)



def subcluster():
    _, context, _ = st.columns([0.5, 8, 1.2])
    _.button(label='↩ 返回', on_click=callback, type='primary')
    with context:
        st.subheader("🍓 Subclustering scRNA-seq datasets")
        st.write("After we have completed the scRNA-seq workflow and identified the various cell types present in our samples, we might decide that for a particular cell type, we would like to identify subtypes. For example, if we have a large cluster of CD4+ Helper T cells, we may want to identify subsets of Th1, Th2, Th17, Th9, and Tfh cells. To identify these cell subsets, we would subset the dataset to the cell type(s) of interest (e.g. CD4+ Helper T cells). ")
        uploaded_file = st.file_uploader(label="h5ad文件上传:", type='h5ad')
        if uploaded_file:
            adata = sc.read_h5ad(uploaded_file)
            adata.uns['log1p']['base'] = None
            st.write(adata)
            columns = adata.obs.columns
            types = [i for i in columns if i=='cell_type' or i.startswith('leiden_')]
            method = pills("", types, ["🎯", "🍊"], label_visibility="collapsed")
            form(adata, method, context)
        else:
            st.info('上传.h5ad文件, 选择细胞类型或细胞Cluster继续亚群细分。', icon="ℹ️")


   

