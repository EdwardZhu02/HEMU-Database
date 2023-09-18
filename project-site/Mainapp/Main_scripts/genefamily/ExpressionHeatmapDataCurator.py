import pandas as pd
import numpy as np
import pymysql
from Mainapp.Main_scripts.MainConfiguration import query_sql

# Used for pymysql services
dbhost = str(query_sql("host"))
dbuser = str(query_sql("user"))
dbpassword = str(query_sql("pwd"))
dbdatabase = str(query_sql("dbname"))


def gene_expheatmap_df_builder(query_sheet_name, query_gene_list, query_sample_list, query_format):
    """

    :param query_sheet_name:
    :param query_gene_list:
    :param query_sample_list:
    :param query_format: FPKM / TPM
    :return: dataframe, [sample_id, sample_tissue, gene_id, tf_fam, tpm]
    """
    final_df_list = []

    db = pymysql.connect(host=dbhost, user=dbuser, password=dbpassword, database=dbdatabase)
    cursor = db.cursor()

    for indv_sample in query_sample_list:
        for indv_gene in query_gene_list:

            sqlcmd = "SELECT * FROM %s WHERE gene='%s' and sample_id='%s';" % \
                (query_sheet_name, indv_gene, indv_sample)
            
            try:  # Execute SQL command
                cursor.execute(sqlcmd)
                results = cursor.fetchall()
                
                if not results:
                    # No results generated from database query
                    print("Query error:", sqlcmd, sep=" ")

                if query_format == "FPKM":
                    # [[Gene_id, Sample_id, FPKM, tissue_type], [G2, S2, F2, t2]..]
                    tmp_result_list = [[indv[0], indv[1], indv[2], indv[4]] for indv in results]
                elif query_format == "TPM":
                    # [[Gene_id, Sample_id, TPM, tissue_type], [G2, S2, T2, t2]..]
                    tmp_result_list = [[indv[0], indv[1], indv[3], indv[4]] for indv in results]

                for indv_temp_result in tmp_result_list:
                    final_df_list.append(indv_temp_result)

            except:
                print("Exception occurred while querying database.")

    # Build dataframe
    init_df = pd.DataFrame(final_df_list,
                           columns=["gene_id", "sample_id", "exp_level", "tissue"])
    return init_df
