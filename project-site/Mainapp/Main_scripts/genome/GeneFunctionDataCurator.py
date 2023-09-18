#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main 
@File    ：GeneFunctionDataCurator.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2023/04/02 15:09 
@IDE     ：PyCharm 
-----------------------------
Description:
Input:
Output:

One-line usage note: python3 GeneFunctionDataCurator.py ,
- ARG1:

Environment requirement (if any): base
- Dependencies: None
-----------------------------
Notes:

-----------------------------
File Revisions:
    2023/04/02: Version 1 - Creation
"""
import pandas as pd
import pymysql
from Mainapp.Main_scripts.MainConfiguration import query_sql

# Extract gene functional annotation from database
# Database configurations

# Used for pymysql services
dbhost = str(query_sql("host"))
dbuser = str(query_sql("user"))
dbpassword = str(query_sql("pwd"))
dbdatabase = str(query_sql("dbname"))


def gene_functional_annot_builder(gene_id, gokegg_sheet_name):

    db = pymysql.connect(host=dbhost, user=dbuser, password=dbpassword, database=dbdatabase)
    cursor = db.cursor()

    gene_id_query = "%" + gene_id + "%"  # ambiguous query
    sqlcmd_select_gokegg = "SELECT * FROM %s WHERE query LIKE '%s';" % (gokegg_sheet_name, gene_id_query)
    try:  # Execute SQL command
        cursor.execute(sqlcmd_select_gokegg)
        results = cursor.fetchall()
    except:
        print("Exception occurred while querying database.")
        return RuntimeError

    indv_gokegg_list = [[indv[0], indv[1], indv[2], indv[3], indv[4]] for indv in results][0]
    # gene_transcript_name, description, GOterm, KEGG_kterm, KEGG_pathway

    return indv_gokegg_list
