#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main 
@File    ：genefam_phylogenetic_plt_generator.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2023/02/27 16:04 
@IDE     ：PyCharm 
-----------------------------
Description:
Input:
Output:

One-line usage note: python3 genefam_phylogenetic_plt_generator.py ,
- ARG1:

Environment requirement (if any): base
- Dependencies: None
-----------------------------
Notes:

-----------------------------
File Revisions:
    2023/02/27: Version 1 - Creation
"""
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import random
import os

# Used for passing pandas dataframe variable to R
pandas2ri.activate()  # Used for variable conversion (pd.DataFrame)


def genefam_phylogenetic_analysis(sequence_fasta_string,
                                  msa_method, pairwise_dist_method, tree_layout_method, bootstrap_rep_num,
                                  is_protein_sequence):
    # :param is_protein_sequence: True / False

    # Generate unique folder name and create the folder
    folder_name = "phylo" + str(random.randint(int(1e8), int(1e9) - 1))
    os.makedirs("Mainapp/static/Temp_R_genefam/" + folder_name, exist_ok=True)

    # Write sequence fasta string to seq_original.fasta, which can be identified by R script.
    with open("Mainapp/static/Temp_R_genefam/" + folder_name + "/seq_original.fasta", mode='w') as out_fh:
        out_fh.write(sequence_fasta_string)

    # --- Read R script ---
    if is_protein_sequence:
        rscript_fh = open("Mainapp/R_scripts/msa_tree_construction_protein.R")
    else:
        rscript_fh = open("Mainapp/R_scripts/msa_tree_construction_nucleotide.R")
    # readlines() generates a list, concat pending.
    rscript = "".join(rscript_fh.readlines())
    rscript_fh.close()
    robjects.r(rscript)
    # --- Finish reading R script ---

    robjects.globalenv['wkdir_name'] = folder_name
    robjects.globalenv['msa_method'] = msa_method
    robjects.globalenv['pairwisedist_method'] = pairwise_dist_method
    robjects.globalenv['phylotree_layout'] = tree_layout_method
    robjects.globalenv['bootstrap_rep'] = bootstrap_rep_num

    # Execute function
    robjects.r['genefam_phylogenetic_analysis'](
        folder_name, msa_method, pairwise_dist_method, tree_layout_method, bootstrap_rep_num,
    )

    return folder_name
