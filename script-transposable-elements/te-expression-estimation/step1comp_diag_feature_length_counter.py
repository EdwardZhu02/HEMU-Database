#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_scripts 
@File    ：step1comp_diag_feature_length_counter.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2023/03/20 19:47 
@IDE     ：PyCharm 
-----------------------------
Description:
Input:
Output:

One-line usage note: python3 step1comp_diag_feature_length_counter.py ,
- ARG1:

Environment requirement (if any): base
- Dependencies: None
-----------------------------
Notes:

-----------------------------
File Revisions:
    2023/03/20: Version 1 - Creation
"""


# Count feature length using converted GTF files


def length_counter_revised(TE_GTF, genome_GTF, TEout_lenTable, Geneout_lenTable):
    """

    :param Geneout_lenTable:
    :param TEout_lenTable:
    :param TE_GTF:
    :param genome_GTF:
    :return:
    """
    fam_length_dict = {}  # fam_name -> total_len
    gene_length_dict = {}  # gene_name -> sum of all exon len
    print("Processing TE GTF")
    try:
        with open(TE_GTF, mode='r') as in_fh:
            for entry in in_fh:
                entry_split_list = entry.split("\t")

                if int(entry_split_list[4]) - int(entry_split_list[3]) <= 1:
                    print("Detected blank or single-nucleotide annotation, skipping.\n", entry)
                    continue

                gtf_entry_infoarea_split_list = entry_split_list[-1].rstrip("\n").split(";")
                gtf_entry_infoarea_dict = {}
                for indv_info in gtf_entry_infoarea_split_list:
                    if len(indv_info.split("\"")) == 3:
                        gtf_entry_infoarea_dict[indv_info.split('\"')[0].rstrip(' ').lstrip(' ')] = \
                            indv_info.split('\"')[1]
                # print(entry)
                # print(gtf_entry_infoarea_dict.keys())
                try:
                    if fam_length_dict[gtf_entry_infoarea_dict['family_id']]:
                        fam_length_dict[gtf_entry_infoarea_dict['family_id']] = \
                            int(fam_length_dict[gtf_entry_infoarea_dict['family_id']]) + \
                            (int(entry_split_list[4]) - int(entry_split_list[3]))
                except KeyError:
                    fam_length_dict[gtf_entry_infoarea_dict['family_id']] = int(entry_split_list[4]) - int(
                        entry_split_list[3])
    except FileNotFoundError:
        # File not found
        print("Exception when opening TE GTF")
        pass

    print("Processing gene GTF")
    try:
        with open(genome_GTF, mode='r') as in_fh:
            for entry in in_fh:
                entry_split_list = entry.split("\t")
                if entry_split_list[2] == 'exon':

                    if int(entry_split_list[4]) - int(entry_split_list[3]) <= 1:
                        print("Detected blank or single-nucleotide annotation, skipping.\n", entry)
                        continue

                    gtf_entry_infoarea_split_list = entry_split_list[-1].rstrip("\n").split(";")
                    gtf_entry_infoarea_dict = {}
                    for indv_info in gtf_entry_infoarea_split_list:
                        if len(indv_info.split("\"")) == 3:
                            gtf_entry_infoarea_dict[indv_info.split('\"')[0].rstrip(' ').lstrip(' ')] = \
                                indv_info.split('\"')[1]

                    try:
                        if gene_length_dict[gtf_entry_infoarea_dict['gene_id']]:
                            # print("added, %s %s + (%s -> %s)" % (gtf_entry_infoarea_dict['gene_id'],
                            # gene_length_dict[gtf_entry_infoarea_dict['gene_id']], entry_split_list[3],
                            # entry_split_list[4]))
                            gene_length_dict[gtf_entry_infoarea_dict['gene_id']] = \
                                int(gene_length_dict[gtf_entry_infoarea_dict['gene_id']]) + \
                                (int(entry_split_list[4]) - int(entry_split_list[3]))
                    except KeyError:
                        # print("created, %s, %s -> %s = %d" % (
                        # gtf_entry_infoarea_dict['gene_id'], entry_split_list[3], entry_split_list[4],
                        # (int(entry_split_list[4]) - int(entry_split_list[3]))))

                        gene_length_dict[gtf_entry_infoarea_dict['gene_id']] = (int(entry_split_list[4]) - int(
                            entry_split_list[3]))
    except FileNotFoundError:
        # File not found
        print("Exception when opening gene GTF")
        pass

    print("Writing output")
    with open(TEout_lenTable, mode='w') as out_fh:
        for ename, elen in fam_length_dict.items():
            out_fh.write(ename + "\t" + str(elen) + "\n")
    with open(Geneout_lenTable, mode='w') as out_fh:
        for ename, elen in gene_length_dict.items():
            out_fh.write(ename + "\t" + str(elen) + "\n")
    print("All done")


if __name__ == "__main__":
    length_counter_revised(
        "./Module6_TE_expression/TE_annot_converted.gtf",
        "./Module6_TE_expression/genome_annot_converted.gtf",
        "./Module6_TE_expression/TE_family.length",
        "./Module6_TE_expression/genome_gene.length",
    )
