import pandas as pd
import pymysql
from Mainapp.Main_scripts.MainConfiguration import query_sql

# Extract te expression data from database, by individual te family ID

# Used for pymysql services
dbhost = str(query_sql("host"))
dbuser = str(query_sql("user"))
dbpassword = str(query_sql("pwd"))
dbdatabase = str(query_sql("dbname"))


def TE_exp_df_builder_family(TE_family_id, TE_sheet_name):
    """

    :return: pd.DataFrame object, with column index: [sample_id, FPKM, TPM, sample_tissue]
    """
    global colnames, results, init_df_list, result_tissue

    init_df_list = []  # List for generating initial dataframe (sample_id , FPKM, TPM)
    db = pymysql.connect(host=dbhost, user=dbuser, password=dbpassword, database=dbdatabase)
    cursor = db.cursor()

    sqlcmd_select = "SELECT * FROM %s WHERE TE_id='%s';" % (TE_sheet_name, TE_family_id)

    try:  # Execute SQL command
        cursor.execute(sqlcmd_select)
        results = cursor.fetchall()
    except:
        print("Exception occurred while querying database.")

    try:
        # Obtain te classification information
        TE_class = results[0][3]
        TE_class_group = results[0][6]
    except IndexError:
        # Possible result of neglecting an "_INT" suffix when querying LTR-like TEs
        # re-do query
        sqlcmd_select = "SELECT * FROM %s WHERE TE_id='%s_INT';" % (TE_sheet_name, TE_family_id)

        try:  # Execute SQL command
            cursor.execute(sqlcmd_select)
            results = cursor.fetchall()
        except:
            print("Exception occurred while querying database.")

        try:
            # Obtain te classification information
            TE_class = results[0][3]
            TE_class_group = results[0][6]
        except IndexError:
            pass

    # Build dataframe, columns=['sample_id', 'fpkm', 'tpm', 'sample_tissue']
    init_df_list = [[indv[1], indv[2], indv[3], indv[6]] for indv in results]
    init_df = pd.DataFrame(init_df_list,
                           columns=['sample_id', 'fpkm', 'tpm', 'sample_tissue'])

    return TE_class, TE_class_group, init_df
