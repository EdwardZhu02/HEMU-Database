"""
这是最终生成大表达量矩阵的脚本文件。
目标物种为薏苡（Coix lacryma）.
使用方法
python Expression_Database_Generator [./directory_of_tab_storage] [./directory/Output.csv]
"""

import os
import re
import sys

import pandas as pd
import numpy as np
from pandas.testing import assert_frame_equal

tab_list = []


def obtain_tabs(directory):
    file_list = []
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            file_list.append(str(os.path.join(root, dir)))
        for file in files:
            file_list.append(str(os.path.join(root, file)))
    for file in file_list:
        if file.endswith(".tab"):
            tab_list.append(file)


def parse_tabs(output_dir):
    dataframe_all_comp = []
    flag1 = 0

    for tab_abs_dir in tab_list:

        flag2 = 0
        sample_indv_database = []

        (filepath, filename) = os.path.split(tab_abs_dir)
        tmp_match_obj = re.search(r'(SRR\d+)V4_Gene_Abundance.tab', filename)
        if tmp_match_obj:
            tab_SRR_name = tmp_match_obj.group(1)
        else:
            exit()
        # tab_SRR_name = str(filename.rstrip("V4_Gene_Abundance.tab"))
        print("Processing %s" % tab_SRR_name)

        try:
            with open(tab_abs_dir, "r") as tab_lines:
                for tab_line in tab_lines:
                    tab_line_sep = tab_line.split("\t")

                    match_object = re.search(r"(Cl\w+)", tab_line_sep[0])
                    if match_object:
                        sample_indv_database.append(
                            [match_object.group(1), tab_line_sep[7], tab_line_sep[8].rstrip("\n")])
                        # [7]-FPKM [8]-TPM
        except IOError:
            return ("IOError")

        if sample_indv_database:
            indv_df = pd.DataFrame(sample_indv_database,
                                   columns=["Genes", tab_SRR_name + "_FPKM", tab_SRR_name + "_TPM"])
            indv_df.set_index("Genes")
            indv_df = indv_df.sort_values(by=["Genes", tab_SRR_name + "_FPKM"], ascending=False)
            # indv_df.to_csv("%s.csv" % tab_SRR_name, index="")

            indv_df = indv_df.drop_duplicates(["Genes"])  # 去除含有重复Index的行
            # print(indv_df)

            if flag1 == 0:  # 生成一个含有全部不重复基因的DataFrame
                df_gene_cut = indv_df[["Genes"]]  # 单取基因的一行
                df_gene_cut = df_gene_cut.reset_index(drop=True)
                dataframe_all_comp.append(df_gene_cut)
                flag1 = 1

            if flag1 == 1:  # 比较是否与第一次生成的基因DataFrame相同，避免数据对齐错误
                df_gene_cut2 = indv_df[["Genes"]]  # 单取基因的一行
                df_gene_cut2 = df_gene_cut.reset_index(drop=True)

                as_res = assert_frame_equal(df_gene_cut, df_gene_cut2)
                if not as_res:
                    flag2 = 1

            if flag2 == 1:
                print("%s Gene id comp CHECKSUM COMPLETE" % tab_SRR_name)
                df_cut = indv_df[[tab_SRR_name + "_FPKM", tab_SRR_name + "_TPM"]]
                df_cut = df_cut.reset_index(drop=True)
                dataframe_all_comp.append(df_cut)
            else:
                print("%s Gene id comp CHECKSUM FAIL. SKIPPING." % tab_SRR_name)

    try:
        print("---\nGenerating Output CSV....\n---")
        df_all = pd.concat(dataframe_all_comp, axis=1)
        df_all.to_csv(output_dir, index="")
        print("---\nOutput CSV Generated: %s\n---" % output_dir)

    except PermissionError:
        print("CSV OUTPUT:Permission Error")
    except IOError:
        print("CSV OUTPUT:IO Error")


if __name__ == "__main__":
    obtain_tabs(sys.argv[1])
    # obtain_tabs("./Maize_V4_Tabs/")

    parse_tabs(sys.argv[2])
    # parse_tabs("./Main.csv")
