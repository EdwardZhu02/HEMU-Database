#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main 
@File    ：ATACPeakDataCurator.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2023/04/10 23:28 
@IDE     ：PyCharm 
-----------------------------
Description:
Input:
Output:

One-line usage note: python3 ATACPeakDataCurator.py ,
- ARG1:

Environment requirement (if any): base
- Dependencies: None
-----------------------------
Notes:

-----------------------------
File Revisions:
    2023/04/10: Version 1 - Creation
"""

import pandas as pd
import pymysql
from Mainapp.Main_scripts.MainConfiguration import query_sql
import re

# Extract gene functional annotation from database
# Database configurations

# Used for pymysql services
dbhost = str(query_sql("host"))
dbuser = str(query_sql("user"))
dbpassword = str(query_sql("pwd"))
dbdatabase = str(query_sql("dbname"))


def atac_peak_df_builder(query_species, query_flankingregion_length, query_table_narrowpeak, query_table_gene, mainquery, querytype):
    """

    :param query_table_gene:
    :param query_table_narrowpeak:
    :param query_species:
    :param query_flankingregion_length:
    :param mainquery: geneid -> list, containing all gene IDs / range -> [Seq]:[start]..[end]
    :param querytype: geneid / range
    :return: gene_df, peak_df
    """

    db = pymysql.connect(host=dbhost, user=dbuser, password=dbpassword, database=dbdatabase)
    cursor = db.cursor()

    if querytype == "geneid":
        # [geneid1, geneid2, ..]

        # [["geneid", "seqid", "source", "ontology", "start", "end", "score", "strand", "phase", "attributes"], [], ...]
        gene_total_entry_list = []
        # [["chr", "start", "end", "name", "score", "strand", "signal_value", "-log10pvalue", "-log10qvalue", "peak"],.]
        peak_total_entry_list = []

        for indv_geneid in mainquery:

            sqlcmd_select_gene = "SELECT * FROM %s WHERE geneid='%s';" % (
                query_table_gene, indv_geneid
            )
            try:  # Execute SQL command
                cursor.execute(sqlcmd_select_gene)
                results_gene = cursor.fetchall()
            except:
                print("Exception occurred while querying database.")
                return RuntimeError
            # result_gene_list: A 10-col GFF-derived format, containing only the gene line.
            # ["geneid", "seqid", "source", "ontology", "start", "end", "score", "strand", "phase", "attributes"]
            result_gene_list = [list(indv[0:len(indv)]) for indv in results_gene][0]

            if not result_gene_list:
                # Query invalid, no gene found in the database sheet
                return [101]

            # Add gene entry to the final return list
            gene_total_entry_list.append(result_gene_list)

            # Query peak information regarding the target gene
            seqname = result_gene_list[1]
            seqstart = result_gene_list[4]
            seqend = result_gene_list[5]

            sqlcmd_select_overlap_peak = "SELECT * FROM %s WHERE seqid='%s' AND start <= %s AND end >= %s;" % (
                query_table_narrowpeak, seqname,
                str(int(seqend) + int(query_flankingregion_length)),  # End
                str(int(seqstart) - int(query_flankingregion_length)),  # Start
            )
            try:  # Execute SQL command
                cursor.execute(sqlcmd_select_overlap_peak)
                results_peak = cursor.fetchall()
            except:
                print("Exception occurred while querying database.")
                return RuntimeError
            # ["chr", "start", "end", "name", "score", "strand", "signal_value", "-log10pvalue", "-log10qvalue", "peak"]
            tmp_peak_list = [list(indv[0:len(indv)]) for indv in results_peak]
            for indv_peak_entry in tmp_peak_list:
                peak_total_entry_list.append(indv_peak_entry)

        # End of for loop
        return [gene_total_entry_list, peak_total_entry_list]

    elif querytype == "range":
        # [Seq]:[start]..[end]

        match_result = re.search(r'^(?P<seqname>.+):(?P<seqstart>\d+)..(?P<seqend>\d+)$', mainquery)
        if match_result:
            seqname = match_result.group("seqname")
            seqstart = match_result.group("seqstart")
            seqend = match_result.group("seqend")

            sqlcmd_select_overlap_peak = "SELECT * FROM %s WHERE seqid='%s' AND start <= %s AND end >= %s;" % (
                query_table_narrowpeak, seqname,
                str(int(seqend) + int(query_flankingregion_length)),  # End
                str(int(seqstart) - int(query_flankingregion_length)),  # Start
            )
            try:  # Execute SQL command
                cursor.execute(sqlcmd_select_overlap_peak)
                results = cursor.fetchall()
            except:
                print("Exception occurred while querying database.")
                return RuntimeError

            # ["chr", "start", "end", "name", "score", "strand", "signal_value", "-log10pvalue", "-log10qvalue", "peak"]
            return_list = [list(indv[0:len(indv)]) for indv in results]
            return return_list

        else:
            # Invalid entry format, not captured by regex
            return [101]
    else:
        # Invalid query type
        return RuntimeError
