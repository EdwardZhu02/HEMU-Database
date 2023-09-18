import pandas as pd
import pymysql
from Mainapp.Main_scripts.MainConfiguration import query_sql

# Extract gene expression data from database, used for differential expression analysis
# Database configurations

# Used for pymysql services
dbhost = str(query_sql("host"))
dbuser = str(query_sql("user"))
dbpassword = str(query_sql("pwd"))
dbdatabase = str(query_sql("dbname"))


def wgcna_init_df_builder(exp_sheet_name, gene_list, sample_list, expression_format):
    """

    :param exp_sheet_name:
    :param gene_list:
    :param sample_list:
    :param expression_format: FPKM/TPM
    :return:
    """
    global results, final_df, indv_geneid
    # print(gene_list, sample_list)
    db = pymysql.connect(host=dbhost, user=dbuser, password=dbpassword, database=dbdatabase)
    cursor = db.cursor()

    total_df_dict = {}
    for indv_sampleid in sample_list:
        tmp_exp_list = []
        if expression_format == "FPKM":
            for indv_geneid in gene_list:
                sqlcmd = "SELECT fpkm FROM %s WHERE sample_id='%s' AND gene='%s'" % (
                    exp_sheet_name, indv_sampleid, indv_geneid)
                try:  # Execute SQL command
                    cursor.execute(sqlcmd)
                    results = cursor.fetchall()
                    tmp_exp_list.append(float(results[0][0]))
                except:
                    print("Exception occurred while querying database. \n%s" % sqlcmd)
        elif expression_format == "TPM":
            for indv_geneid in gene_list:
                sqlcmd = "SELECT tpm FROM %s WHERE sample_id='%s' AND gene='%s'" % (
                    exp_sheet_name, indv_sampleid, indv_geneid)
                try:  # Execute SQL command
                    cursor.execute(sqlcmd)
                    results = cursor.fetchall()
                    tmp_exp_list.append(float(results[0][0]))
                except:
                    print("Exception occurred while querying database. \n%s" % sqlcmd)
        else:
            continue

        total_df_dict[indv_sampleid] = tmp_exp_list

    final_df = pd.DataFrame(total_df_dict)
    final_df.index = gene_list

    return final_df
