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


def gene_exp_df_builder(gene_id, exp_sheet_name):
    """

    :param gene_id:
    :param exp_sheet_name:
    :param sampleinfo_sheet_name:
    :return: pd.DataFrame object, with column index: [sample_id, FPKM, TPM, sample_tissue]
    """
    global colnames, results, init_df_list, result_tissue

    init_df_list = []  # List for generating initial dataframe (sample_id , FPKM, TPM)
    db = pymysql.connect(host=dbhost, user=dbuser, password=dbpassword, database=dbdatabase)
    cursor = db.cursor()

    sqlcmd_select_gene = "SELECT * FROM %s WHERE gene='%s';" % (exp_sheet_name, gene_id)

    try:  # Execute SQL command
        # print(sqlcmd_select_gene)
        cursor.execute(sqlcmd_select_gene)
        results = cursor.fetchall()
    except:
        print("Exception occurred while querying database.")
        return RuntimeError

    # Build dataframe
    init_df_list = [[indv[1], indv[2], indv[3], indv[4]] for indv in results]
    init_df = pd.DataFrame(init_df_list, columns=['sample_id', 'fpkm', 'tpm', 'tissue_type'])
    # SRR-ID FPKM TPM TISSUE INFO1 INFO2 INFO3

    return init_df


def fund_info_obtainer(gene_id, gene_exp_df, exp_format):
    """
    Obtain critical information regarding gene expression, Used in the right-side panel of RNA-seq search result.
    :param gene_id: gene ID
    :param gene_exp_df: sample_id FPKM TPM tissue_type
    :param exp_format: FPKM / TPM, used for extracting expression level from dataframe
    :return: gene_exp_detail = [gene_ID, expressed_samples, total_samples, max, min, median]
    """
    gene_exp_detail = []
    expressed_threshold = 1  # Modifiable, default with both TPM and FPKM.

    if exp_format == "TPM":
        exp_level_comp = gene_exp_df['tpm'].tolist()
        gene_exp_detail.append([
            gene_id,
            len([indv for indv in exp_level_comp if indv >= expressed_threshold]),
            len(exp_level_comp),
            max(exp_level_comp),
            min(exp_level_comp),
            np.median(exp_level_comp),
        ])
    elif exp_format == "FPKM":
        exp_level_comp = gene_exp_df['fpkm'].tolist()
        gene_exp_detail.append([
            gene_id,
            len([indv for indv in exp_level_comp if indv >= expressed_threshold]),
            len(exp_level_comp),
            max(exp_level_comp),
            min(exp_level_comp),
            np.median(exp_level_comp),
        ])

    return gene_exp_detail
