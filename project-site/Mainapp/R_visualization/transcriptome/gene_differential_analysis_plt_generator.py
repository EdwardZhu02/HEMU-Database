import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import pandas as pd  # for exporting raw expression value csv
import os
import random

# Used for passing pandas dataframe variable to R
pandas2ri.activate()  # Used for variable conversion (pd.DataFrame)


def GeneDifferentialAnalysis(DE_data_raw, DE_group_list, DE_group_color_list,
                             logfc_threshold, pvalue_threshold, heatmap_gene_count,
                             group1_name, group2_name):
    # Generate unique folder name and create the folder
    output_folder_name = "dge" + str(random.randint(int(1e8), int(1e9) - 1))

    os.makedirs('Mainapp/static/Temp_R_DEprojects/' + output_folder_name, exist_ok=True)
    try:
        DE_data_raw.to_csv('Mainapp/static/Temp_R_DEprojects/' + output_folder_name + '/raw_exp_values.csv')
    except:
        return 102  # Exception when writing raw expression value csv

    robjects.globalenv['DE_data_raw'] = DE_data_raw
    robjects.globalenv['DE_group_list'] = DE_group_list
    robjects.globalenv['DE_group_color_list'] = DE_group_color_list
    robjects.globalenv['logfc_threshold'] = logfc_threshold
    robjects.globalenv['pvalue_threshold'] = pvalue_threshold
    robjects.globalenv['heatmap_gene_count'] = heatmap_gene_count
    robjects.globalenv['group1_name'] = group1_name
    robjects.globalenv['group2_name'] = group2_name
    robjects.globalenv['output_folder_name'] = output_folder_name

    # --- Read R script ---
    rscript_fh = open("Mainapp/R_scripts/dge_report_plotter.R")
    # readlines() generates a list, concat pending.
    rscript = "".join(rscript_fh.readlines())
    rscript_fh.close()
    robjects.r(rscript)
    # --- Finish reading R script ---

    robjects.r['DE_analysis'](DE_data_raw, DE_group_list, DE_group_color_list, logfc_threshold, pvalue_threshold,
                              heatmap_gene_count, group1_name, group2_name, output_folder_name)  # Execute function

    return output_folder_name
