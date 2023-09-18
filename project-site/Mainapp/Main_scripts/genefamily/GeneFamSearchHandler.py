#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main 
@File    ：GeneFamSearchHandler.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2023/03/02 20:24 
@IDE     ：PyCharm 
-----------------------------
Description:
Input:
Output:

One-line usage note: python3 GeneFamSearchHandler.py ,
- ARG1:

Environment requirement (if any): base
- Dependencies: None
-----------------------------
Notes:

-----------------------------
File Revisions:
    2023/03/02: Version 1 - Creation
"""
import os
import random
import subprocess


def hmmsearch_deployer(species_filename, hmm_query_filename, seq_evalue_threshold, dom_evalue_threshold):

    if not os.path.exists("global-static/pfam_hmm/" + hmm_query_filename):
        return 102  # .hmm file not found

    # Generate unique folder name and create the folder
    output_file_name = "hmm" + str(random.randint(int(1e8), int(1e9) - 1)) + ".txt"

    # Generate full paths to input and output
    pep_input = "global-static/annot-sequences/%s.proteins.fasta" % species_filename
    hmm_input = "global-static/pfam_hmm/" + hmm_query_filename
    text_output = "Mainapp/static/Temp_R_genefam/" + output_file_name

    # hmmsearch [options] <hmmfile> <seqdb>
    hmmsearch_instructions = "hmmsearch -E %s --domE %s --cpu 1 %s %s > %s" % (
        seq_evalue_threshold, dom_evalue_threshold, hmm_input, pep_input, text_output
    )
    # print(hmmsearch_instructions)
    p = subprocess.Popen(str(hmmsearch_instructions), shell=True)
    return_code = p.wait()

    if return_code:
        return 103  # Program Failed

    return output_file_name
