import pandas as pd
import numpy as np
import pymysql
from Mainapp.Main_scripts.MainConfiguration import query_sql

# Extract gene expression data from database and append tissue information regarding each sample
# Database configurations

# Used for pymysql services
dbhost = str(query_sql("host"))
dbuser = str(query_sql("user"))
dbpassword = str(query_sql("pwd"))
dbdatabase = str(query_sql("dbname"))


def gene_gokegg_df_builder(gene_id_list, gokegg_sheet_name):
    """

    :return: pd.DataFrame object, with column: [gene_transcript_name, description, GOterm, KEGG_kterm, KEGG_pathway]
    """

    global results, final_gokegg_df
    db = pymysql.connect(host=dbhost, user=dbuser, password=dbpassword, database=dbdatabase)
    cursor = db.cursor()
    df_init_flag = False
    final_gokegg_df = pd.DataFrame([])

    for gene_id in gene_id_list:
        gene_id_query = "%" + gene_id + "%"
        sqlcmd_select_gokegg = "SELECT * FROM %s WHERE query LIKE '%s';" % (gokegg_sheet_name, gene_id_query)
        try:  # Execute SQL command
            cursor.execute(sqlcmd_select_gokegg)
            results = cursor.fetchall()
        except:
            print("Exception occurred while querying database.")

        indv_gokegg_list = [[indv[0], indv[1], indv[2], indv[3], indv[4]] for indv in results]
        # gene_transcript_name, description, GOterm, KEGG_kterm, KEGG_pathway

        if not df_init_flag:
            final_gokegg_df = pd.DataFrame(indv_gokegg_list, columns=['query_name', 'Description', 'GOs',
                                                                      'KEGG_ko', 'KEGG_Pathway'])
            df_init_flag = True
        else:
            indv_gokegg_df = pd.DataFrame(indv_gokegg_list, columns=['query_name', 'Description', 'GOs',
                                                                     'KEGG_ko', 'KEGG_Pathway'])
            final_gokegg_df = pd.concat([final_gokegg_df, indv_gokegg_df])

    return final_gokegg_df
