#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main 
@File    ：PhyloAnalysisHandler.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2023/02/27 16:03 
@IDE     ：PyCharm 
-----------------------------
Description:
Input:
Output:

One-line usage note: python3 PhyloAnalysisHandler.py ,
- ARG1:

Environment requirement (if any): base
- Dependencies: None
-----------------------------
Notes:

-----------------------------
File Revisions:
    2023/02/27: Version 1 - Creation
"""
from Mainapp.Main_scripts.transcriptome import GeneSequenceObtainer
from Mainapp.Main_scripts import MainConfiguration


def genefam_data_validation(sequence_query, species_query, is_protein_sequence=False):
    if sequence_query.startswith(">"):
        # is fasta
        sequence_list_final = [indv_gene.rstrip("\r") for indv_gene in sequence_query.split("\n")]
        sequence_list_final = [indv_entry for indv_entry in sequence_list_final if indv_entry != '']
        if len(sequence_list_final) > 0:
            return "\n".join(sequence_list_final)
        else:
            return None
    else:
        # is gene IDs, pending query
        gene_list_final = [indv_gene.rstrip("\r") for indv_gene in sequence_query.split("\n")]
        gene_list_final = [indv_entry for indv_entry in gene_list_final if indv_entry != '']

        # Obtain filenames for fasta file containing gene sequences
        species_filename = MainConfiguration.query_tables(species_query + "_filename")

        if is_protein_sequence:
            return_fasta_list = GeneSequenceObtainer.gene_sequence_query(
                gene_list_final, species_filename, "protein", False, True
            )
        else:
            return_fasta_list = GeneSequenceObtainer.gene_sequence_query(
                gene_list_final, species_filename, "cds", False, True
            )
        # The last two parameters: write_fasta=False, primary_seq_req=False (default)
        # This sequence query is used for Multiple Sequence Alignment
        # so only the primary sequence should be displayed (primary_seq_req=True)
        # CDS is used for phylogenetic analysis.
        if return_fasta_list == 102:
            # Empty sequence list
            return None
        if len(return_fasta_list) > 0:
            return "".join(return_fasta_list)
        else:
            return None

