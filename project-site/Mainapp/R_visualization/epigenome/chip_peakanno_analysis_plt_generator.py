import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import os
import random

# Used for passing pandas dataframe variable to R
pandas2ri.activate()  # Used for variable conversion (pd.DataFrame)


def SingleChIPAnalysis(sample_id, tssRegion, species,
                       ignore_1st_exon, ignore_1st_intron, ignore_downstream, ignore_promoter_subcategory):

    def PeakFileValidation(_sample_id, _query_species):
        tmp_peakfile = '%sV4.bam_deduplicate_bam_peaks.narrowPeak' % _sample_id
        target_narrowPeak = "global-static/chip-narrowpeak/%s/%s" % (_query_species, tmp_peakfile)
        if not os.path.exists(target_narrowPeak):
            return 102  # narrowPeak file not found
        return target_narrowPeak

    narrowPeak_file = PeakFileValidation(sample_id, species)
    if narrowPeak_file == 102:
        return 102  # narrowPeak file not found

    # Generate unique folder name and create the folder
    output_folder_name = "peakanno" + str(random.randint(int(1e8), int(1e9) - 1))

    os.makedirs('Mainapp/static/Temp_R_epigenome/' + output_folder_name, exist_ok=True)
    # print("Folder made done, %s" % output_folder_name)

    robjects.globalenv['narrowPeak_file'] = narrowPeak_file
    robjects.globalenv['tssRegion'] = tssRegion
    robjects.globalenv['species'] = species
    robjects.globalenv['ignore_1st_exon'] = ignore_1st_exon
    robjects.globalenv['ignore_1st_intron'] = ignore_1st_intron
    robjects.globalenv['ignore_downstream'] = ignore_downstream
    robjects.globalenv['ignore_promoter_subcategory'] = ignore_promoter_subcategory
    robjects.globalenv['output_folder_name'] = output_folder_name

    # --- Read R script ---
    rscript_fh = open("Mainapp/R_scripts/chip_peakanno_plotter.R")
    # readlines() generates a list, concat pending.
    rscript = "".join(rscript_fh.readlines())
    rscript_fh.close()
    robjects.r(rscript)
    # --- Finish reading R script ---

    robjects.r['peak_analysis_single'](
        narrowPeak_file, tssRegion, species, ignore_1st_exon, ignore_1st_intron, ignore_downstream,
        ignore_promoter_subcategory, output_folder_name
    )  # Execute function

    return output_folder_name
