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


def gene_de_df_builder(exp_sheet_name,
                       group1_sampleid_list, group1_name,
                       group2_sampleid_list, group2_name):
    global results, final_df

    final_de_df_dict = {}
    total_sampleid_list = []  # [SRR1, SRR2, SRR3, ..]
    total_groupby_list = []  # [G1, G1, G1, G2, G2, G2]
    total_group_color_list = []  # ["#FFFFFF", "#FFFFFF"...]

    for indv_sample in group1_sampleid_list:
        total_sampleid_list.append(indv_sample)
        total_groupby_list.append(group1_name)
        total_group_color_list.append("#2878B5")  # Deep blue

    for indv_sample in group2_sampleid_list:
        total_sampleid_list.append(indv_sample)
        total_groupby_list.append(group2_name)
        total_group_color_list.append("#9AC9DB")  # Light blue

    db = pymysql.connect(host=dbhost, user=dbuser, password=dbpassword, database=dbdatabase)
    cursor = db.cursor()

    index_build_flag = 0
    sql_or_flag = 0

    sqlcmd = "SELECT * FROM %s WHERE " % exp_sheet_name
    for indv_sampleid in total_sampleid_list:
        if sql_or_flag == 0:
            sqlcmd += "sample_id='%s' " % indv_sampleid
            sql_or_flag = 1
        else:
            sqlcmd += "OR sample_id='%s' " % indv_sampleid
    sqlcmd += ";"

    try:  # Execute SQL command
        cursor.execute(sqlcmd)
        results = cursor.fetchall()
    except:
        print("Exception occurred while querying database.")

    # Build dataframe
    tmp_df_list = [[indv[0], indv[1], indv[3]] for indv in results]  # [[Gene_id, Sample_id, TPM], [G2, S2, T2]..]
    tmp_df_dict = {}  # sample_id => [[gene_id, TPM],[]]
    for indv_entry in tmp_df_list:
        if indv_entry[1] not in tmp_df_dict.keys():
            tmp_df_dict[indv_entry[1]] = []
            tmp_df_dict[indv_entry[1]].append([indv_entry[0], indv_entry[2]])
        else:
            tmp_df_dict[indv_entry[1]].append([indv_entry[0], indv_entry[2]])

    for sample_id in total_sampleid_list:
        # for sample_id, gene_inf_list in tmp_df_dict.items():  # Avoid re-ordering
        gene_inf_list = tmp_df_dict[sample_id]

        if index_build_flag == 0:
            final_df = pd.DataFrame(gene_inf_list, columns=['gene', sample_id])
            final_df = final_df.set_index('gene', drop=True)
            index_build_flag = 1
        elif index_build_flag == 1:
            tmp_df = pd.DataFrame(gene_inf_list, columns=['gene', sample_id])
            tmp_df = tmp_df.set_index('gene', drop=True)
            final_df = pd.concat([final_df, tmp_df], axis=1)

    return final_df, total_groupby_list, total_group_color_list

# if __name__ == "__main__":
# For testing
# gene_de_df_builder('zea_exp',
# ['SRR13587899','SRR13587900', 'SRR13587901', 'SRR13587902', 'SRR13587903', 'SRR13587904'], '',
# '', '')
