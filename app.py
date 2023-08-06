import os
import streamlit as st
from pathlib import Path
import streamlit_book as stb

from streamlit_option_menu import option_menu


st.set_page_config(
    page_title="BioTools",
    page_icon=":dna:",
    layout="wide",
)

# set page width
st.markdown(
    f"""
        <style>
        .main .block-container {{
            max-width: 100%;
        }}
        </style>
    """,
    unsafe_allow_html=True,
)

current_path = Path(__file__).parent.absolute()

stb.set_book_config(
                    menu_title="Tools",
                    menu_icon="bi bi-tools",
                    options=["Pipeline", "scRNA-seq", "TCGA"],
                    paths=[os.path.join(current_path, "page/pipeline.py"),
                           os.path.join(current_path, "page/scRNA-seq.py"),
                           os.path.join(current_path, "page/tcga.py")],
                    icons=['bi bi-pc-display-horizontal', "bi bi-flower1", "bi bi-puzzle"],
                    save_answers=False,
                   )












