#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main 
@File    ：genefam_exp_heatmap_generator.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2023/03/12 19:49 
@IDE     ：PyCharm 
-----------------------------
Description:
Input:
Output:

One-line usage note: python3 genefam_exp_heatmap_generator.py ,
- ARG1:

Environment requirement (if any): base
- Dependencies: None
-----------------------------
Notes:

-----------------------------
File Revisions:
    2023/03/12: Version 1 - Creation
"""
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import random
import os

# Used for passing pandas dataframe variable to R
pandas2ri.activate()  # Used for variable conversion (pd.DataFrame)


def genefam_heatmap_generator(exp_df, query_format):
    """
    Generate expression level heatmap for interested gene family members
    :param exp_df: dataframe with col names: ["gene_id", "sample_id", "exp_level", "tissue"],
        among which 'exp_level' is the real FPKM or TPM data.
    :param query_format: FPKM / TPM, used for plot labeling
    :return: folder_name
    """
    # Generate unique folder name and create the folder
    folder_name = "htmap" + str(random.randint(int(1e8), int(1e9) - 1))
    os.makedirs("Mainapp/static/Temp_R_genefam/" + folder_name, exist_ok=True)

    # --- Read R script ---
    rscript_fh = open("Mainapp/R_scripts/exp_heatmap_plotter.R")
    # readlines() generates a list, concat pending.
    rscript = "".join(rscript_fh.readlines())
    rscript_fh.close()
    robjects.r(rscript)
    # --- Finish reading R script ---

    robjects.globalenv['exp_df'] = exp_df
    robjects.globalenv['query_format'] = query_format
    robjects.globalenv['wkdir_name'] = folder_name

    # Execute function
    robjects.r['exp_heatmap_plotter'](
        exp_df, query_format
    )

    return folder_name
