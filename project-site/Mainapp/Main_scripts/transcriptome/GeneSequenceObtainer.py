#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main 
@File    ：GeneSequenceObtainer.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2023/02/25 13:38 
@IDE     ：PyCharm 
-----------------------------
Description:
Input:
Output:

One-line usage note: python3 GeneSequenceObtainer.py ,
- ARG1:

Environment requirement (if any): HEMUdb_new
- Dependencies: None
-----------------------------
Notes:

-----------------------------
File Revisions:
    2023/02/25: Version 1 - Creation
"""
import re
import random
from Bio import SeqIO


def gene_sequence_query(gene_id_list, query_species, query_type, write_fasta=False, primary_seq_req=False):
    """

    :param primary_seq_req: Request primary sequences rather than all sequences, used for GO/KEGG enrichment or
        sequence phylogenetic analysis.
    :param gene_id_list:
    :param query_species:
    :param query_type:
    :param write_fasta: Boolean value (True/False), write query result to a single fasta file,
        returning its unique identifier (file in
    static/Temp_R_genefam/[identifier].fasta)
    :return:
    """

    # List containing all the fasta entries, with each element representing a single fasta line.
    output_record_list = []
    if primary_seq_req:
        primary_seq_filename_fill = "primary."
    else:
        primary_seq_filename_fill = ""

    if query_type == "gene":
        target_fasta = "global-static/annot-sequences/%s.gene.%sfasta" % (
            query_species, primary_seq_filename_fill)
    elif query_type == "transcript":
        target_fasta = "global-static/annot-sequences/%s.transcripts.%sfasta" % (
            query_species, primary_seq_filename_fill)
    elif query_type == "protein":
        target_fasta = "global-static/annot-sequences/%s.proteins.%sfasta" % (
            query_species, primary_seq_filename_fill)
    elif query_type == "cds":
        target_fasta = "global-static/annot-sequences/%s.cds.%sfasta" % (
            query_species, primary_seq_filename_fill)
    else:
        return 102  # Final list len = 0
    try:
        with open(target_fasta) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # Filter gene sequences from record ids, below is an example of a record id
                # Cl000212::Chr1:20200-28609(-)
                for indv_entry in gene_id_list:
                    if re.match(r'.*%s.*' % str(indv_entry), str(record.id)):
                        output_record_list.append(record)
                # if record.id.split("::")[0] in gene_id_list:
                #     output_record_list.append(record)
    except FileNotFoundError:
        return 101  # File not found

    if not len(output_record_list) == 0:
        final_return_list = []
        for entry in output_record_list:
            final_return_list.append(">" + str(entry.id) + "\n")
            final_return_list.append(str(entry.seq) + "\n")

        if not write_fasta:
            return final_return_list
        else:
            fasta_file_name = "sequence" + str(random.randint(int(1e8), int(1e9) - 1)) + ".fasta"
            # Write sequence fasta string and return the unique fasta_file_name
            with open("Mainapp/static/Temp_R_genefam/" + fasta_file_name, mode='w') as out_fh:
                out_fh.write("".join(final_return_list))
            return fasta_file_name

    else:
        return 102  # Final list len = 0

