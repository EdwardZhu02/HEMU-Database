import os
import re

import rpy2.robjects as robjects


# Env req: /data5/condaenvs/rpy3dev

class IndividualAccession(object):
    def __init__(self):
        # Root folder for the project, containing all links.
        self._accession_name = os.getcwd().split("/")[-1]
        self._current_dir = os.getcwd() + "/"
        # bam-style alignment file for a certain species
        self._alignment_bam_dir_list = []

    def Count_Normalizar(self):
        # ============================================================================
        # Split and format raw TEcount files
        def cntTable_simplifier(TEcount_cntTable, TEout_cntTable, Geneout_cntTable):
            with open(TEcount_cntTable, mode='r') as in_fh:
                with open(TEout_cntTable, mode='w') as out_fh1:
                    with open(Geneout_cntTable, mode='w') as out_fh2:

                        te_name_dict = {}

                        for entry in in_fh:
                            if re.match(r'^".*"$', entry.split("\t")[0]):  # Gene
                                out_fh2.write("\t".join([entry.split("\t")[0].rstrip("\"").lstrip("\""),
                                                         entry.split("\t")[1]]))

                            elif re.match(r'^.*:.*:.*:.*:.*$', entry.split("\t")[0]):  # TE_longver
                                out_fh1.write("\t".join(
                                    [":".join([entry.split("\t")[0].split(":")[0], entry.split("\t")[0].split(":")[1]]),
                                     entry.split("\t")[1].rstrip("\n"),
                                     entry.split("\t")[0].split(":")[-1]]) + "\n")

                            elif re.match(r'^.*:.*:.*$', entry.split("\t")[0]):  # TE_shortver
                                out_fh1.write(
                                    "\t".join([entry.split("\t")[0].split(":")[0], entry.split("\t")[1].rstrip("\n"),
                                               entry.split("\t")[0].split(":")[2]]) + "\n")
                            else:
                                print("Skipped line, ", entry)
            print("Done processing cntTable")

        # ============================================================================

        # ============================================================================
        # Calculation of FPKM and TPM values
        def fpkm_tpm_calculation(_gene_count, _TE_count, _gene_lentable, _TE_lentable, _outdir):
            """
            All inputs are in dataframes.
            :param _outdir:
            :param _gene_count:
            :param _TE_count:
            :param _gene_lentable:
            :param _TE_lentable:
            :return:
            """

            robjects.globalenv['gene_count'] = _gene_count
            robjects.globalenv['te_count'] = _TE_count
            robjects.globalenv['gene_lentable'] = _gene_lentable
            robjects.globalenv['te_lentable'] = _TE_lentable
            robjects.globalenv['out_dir'] = _outdir

            # --- Read R script ---
            rscript_fh = open("./step2_count_normalize.R")
            # readlines() generates a list, concat pending.
            rscript = "".join(rscript_fh.readlines())
            rscript_fh.close()
            robjects.r(rscript)
            # --- Finish reading R script ---

            # Execute function
            robjects.r['gene_count_normalizer'](_gene_count, _gene_lentable, _outdir)
            # Output1: Normalized FPKM & TPM value - genes
            # Module6_TE_expression/TEcount_final/SRR111111_gene.tab
            robjects.r['te_count_normalizer'](_TE_count, _TE_lentable, _outdir)
            # Output2: Normalized FPKM & TPM value - TEs
            # Module6_TE_expression/TEcount_final/SRR111111_te.tab

        # ============================================================================

        # Search for TEcount-generated cntTable files
        print("Filtering TEcount tables and normalizing CPM to FPKM & TPM.")
        for indv_file in os.listdir(self._current_dir + "Module6_TE_expression/TEcount_raw"):  # ex. 1.cntTable
            if indv_file.endswith(".cntTable"):
                print("Splitting: ", indv_file)
                cntTable_simplifier(self._current_dir + "Module6_TE_expression/TEcount_raw/" + indv_file,
                                    self._current_dir + "Module6_TE_expression/TEcount_normalized/" + indv_file.rstrip(
                                        ".cntTable") + "_TE.count",
                                    self._current_dir + "Module6_TE_expression/TEcount_normalized/" + indv_file.rstrip(
                                        ".cntTable") + "_gene.count",
                                    )
                print("Normalizing Count: ", indv_file)
                fpkm_tpm_calculation(self._current_dir + "Module6_TE_expression/TEcount_normalized/" + indv_file.rstrip(
                    ".cntTable") + "_gene.count",
                                     self._current_dir + "Module6_TE_expression/TEcount_normalized/" + indv_file.rstrip(
                                         ".cntTable") + "_TE.count",
                                     self._current_dir + "Module6_TE_expression/genome_gene.length",
                                     self._current_dir + "Module6_TE_expression/TE_family.length",
                                     self._current_dir + "Module6_TE_expression/TEcount_final/" + indv_file.rstrip(
                                         ".cntTable") + "_"
                                     )
                print("Done:", indv_file)
                # break  # For diagnosis
        print("All done, check 'Module6_TE_expression/TEcount_final/' for TE and gene expression .csv files.")


mainAccession = IndividualAccession()
mainAccession.Count_Normalizar()
