import pandas as pd
import pymysql
from Mainapp.Main_scripts.MainConfiguration import query_sql

# Extract te expression data from database, by individual sample accession

# Used for pymysql services
dbhost = str(query_sql("host"))
dbuser = str(query_sql("user"))
dbpassword = str(query_sql("pwd"))
dbdatabase = str(query_sql("dbname"))


def TE_exp_df_builder_sample(sample_id, TE_sheet_name):
    """

    :return: pd.DataFrame object, with column index: [sample_id, FPKM, TPM, sample_tissue]
    """
    global colnames, results, init_df_list, result_tissue

    init_df_list = []  # List for generating initial dataframe (sample_id , FPKM, TPM)
    db = pymysql.connect(host=dbhost, user=dbuser, password=dbpassword, database=dbdatabase)
    cursor = db.cursor()

    sqlcmd_select = "SELECT * FROM %s WHERE sample_id='%s';" % (TE_sheet_name, sample_id)

    try:  # Execute SQL command
        cursor.execute(sqlcmd_select)
        results = cursor.fetchall()
    except:
        print("Exception occurred while querying database.")

    # Build dataframe, columns=['TE_id', 'fpkm', 'tpm', 'TE_class', 'TE_class_group']
    init_df_list = [[indv[0], indv[2], indv[3], indv[4], indv[5]] for indv in results]
    init_df = pd.DataFrame(init_df_list, columns=['TE_id', 'fpkm', 'tpm', 'TE_class', 'TE_class_group'])
    # print(init_df)
    return init_df
