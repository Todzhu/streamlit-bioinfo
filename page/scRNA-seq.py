import streamlit as st
from page.functions.styles import *
from page.functions.subcluster import *
from page.functions.pseudotime import *
from page.functions.cellproportion import *


if 'tool' not in st.session_state:
    st.session_state['tool'] = None


def selectTool(tool):
    st.session_state['tool'] = tool

if st.session_state['tool'] == '单细胞亚群细分':
    subcluster()
elif st.session_state['tool'] == '单细胞拟时序分析':
    pseudotime()
elif st.session_state['tool'] == '细胞比例差异分析':
    cellproportion()
else:
    st.subheader("🍊 scRNA-seq data visualization tools")
    tool_list = [{'name': '单细胞亚群细分', 'img': img_to_html('./img/subcluster.png'), 'description': '亚群细分是将细胞按照其基因表达模式和功能特征进行进一步的细分和分类，了解细胞群体之间的异质性，发现潜在的细胞亚群。'},
                 {'name': '单细胞拟时序分析', 'img': img_to_html('./img/pseudotime.png'), 'description': 'Monocle2 通过反向图嵌入(Reversed Graph Embedding)根据不同细胞亚群基因表达量随时间的变化情况，构建细胞谱系发育。'},
                 {'name': '细胞比例差异分析', 'img': img_to_html('./img/Roe.png'), 'description': '对单细胞数据进行亚群注释之后，使用Ro/e指标，统计观察到的细胞数与期望细胞数的比值，用于量化每个亚群对组织的偏好程度。'}            
                ]
    N_cards_per_col = 5
    for index, item in enumerate(tool_list):
        i = index%N_cards_per_col
        if i==0:
            st.write("")
            cols = st.columns(N_cards_per_col, gap="large")
        cols[i].button(label=item["name"], key=item["name"], on_click=selectTool, args=[item['name']], type='primary', use_container_width=True)
        cols[i].markdown(item["img"], unsafe_allow_html=True)
        cols[i].markdown(f"**:black[{item['description']}]**")

    add_color_to_cards()

