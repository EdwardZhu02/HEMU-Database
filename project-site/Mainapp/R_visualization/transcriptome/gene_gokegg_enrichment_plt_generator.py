import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import random

# Used for passing pandas dataframe variable to R
pandas2ri.activate()  # Used for variable conversion (pd.DataFrame)


def gokegg_enrich_plt(_gokegg_df, enrich_format, species_filename):
    filename_id = str(random.randint(int(1e7), int(1e8)))
    # For background gene profile generation
    bg_annot_path = "global-static/gokegg-bgannot/" + species_filename + ".emapper.annotations"

    robjects.globalenv['gokegg_dataframe'] = _gokegg_df
    robjects.globalenv['filename_id'] = filename_id
    robjects.globalenv['bg_annot'] = bg_annot_path

    # --- Read R script ---
    rscript_fh = open("Mainapp/R_scripts/gokegg_plotter.R")
    # readlines() generates a list, concat pending.
    rscript = "".join(rscript_fh.readlines())
    rscript_fh.close()
    robjects.r(rscript)
    # --- Finish reading R script ---

    if enrich_format == "GO":
        robjects.r['go_enricher'](_gokegg_df, filename_id, bg_annot_path)  # Execute function
    elif enrich_format == "KEGG":
        robjects.r['kegg_enricher'](_gokegg_df, filename_id, bg_annot_path)  # Execute function
    else:
        # Frontend exception
        return None

    return filename_id
