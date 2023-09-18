import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

# File deletion
import shutil
from pathlib import Path

pandas2ri.activate()  # Used for variable conversion (pd.DataFrame)


def TE_bysample_plt(_te_df, sample_id, query_format):
    # Detect and delete existing plots
    tmp_plot_path = "Mainapp/static/Temp_R_TE"
    for elm in Path(tmp_plot_path).glob(sample_id + '*'):
        elm.unlink() if elm.is_file() else shutil.rmtree(elm)  # delete folder and file

    robjects.globalenv['te_dataframe'] = pd.DataFrame(_te_df)
    robjects.globalenv['filename_id'] = sample_id

    # --- Read R script ---
    rscript_fh = open("Mainapp/R_scripts/te_expression_plotter.R")
    # readlines() generates a list, concat pending.
    rscript = "".join(rscript_fh.readlines())
    rscript_fh.close()
    robjects.r(rscript)
    # --- Finish reading R script ---
    if query_format == "TPM":
        robjects.r['TE_sample_plotter_tpm'](pd.DataFrame(_te_df), sample_id)  # Execute function
    elif query_format == "FPKM":
        robjects.r['TE_sample_plotter_fpkm'](pd.DataFrame(_te_df), sample_id)  # Execute function
    else:
        return RuntimeError
    return sample_id


def TE_byfamily_plt(_te_df, TE_id, query_format):
    # Detect and delete existing plots
    tmp_plot_path = "Mainapp/static/Temp_R_TE"
    for elm in Path(tmp_plot_path).glob(TE_id + '*'):
        elm.unlink() if elm.is_file() else shutil.rmtree(elm)  # delete folder and file

    robjects.globalenv['te_dataframe'] = pd.DataFrame(_te_df)
    robjects.globalenv['filename_id'] = TE_id

    # --- Read R script ---
    rscript_fh = open("Mainapp/R_scripts/te_expression_plotter.R")
    # readlines() generates a list, concat pending.
    rscript = "".join(rscript_fh.readlines())
    rscript_fh.close()
    robjects.r(rscript)
    # --- Finish reading R script ---
    if query_format == "TPM":
        result = robjects.r['TE_family_plotter_tpm'](pd.DataFrame(_te_df), TE_id).rx()  # Execute function
        final_list = [list(i)[0] for i in result]
        # [min_tpm, mean_tpm, max_tpm, sample_number, sample_number_expressed]
    elif query_format == "FPKM":
        result = robjects.r['TE_family_plotter_fpkm'](pd.DataFrame(_te_df), TE_id).rx()  # Execute function
        final_list = [list(i)[0] for i in result]
        # [min_fpkm, mean_fpkm, max_fpkm, sample_number, sample_number_expressed]
    else:
        return RuntimeError

    final_list = [float(i) for i in final_list]
    final_list.append(TE_id)
    # [min_tpm/fpkm, mean_tpm/fpkm, max_tpm/fpkm, sample_number, sample_number_expressed, TE_id]
    return final_list
