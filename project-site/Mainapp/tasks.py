#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main
@File    ：tasks.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2022/10/2 18:00
@IDE     ：PyCharm
-----------------------------
Description: Celery-dependent async tasks
Input:
Output:

Environment requirement (if any): HEMUdb_new (Remote Server, Ubuntu 20.04)
- Dependencies: Django(v3.1), Project files, etc.
-----------------------------
Notes:

-----------------------------
File Revisions:
    2022/10/2: Version 1 - Creation
    2023/2/18: Version 2 - Finished async conversion for all transcriptome-derived analysis
"""

# ===================================================================================
# Global imports
# Async-dependent
from __future__ import absolute_import, unicode_literals
from celery import shared_task
# Function-dependent
import time
import pandas as pd
import random
import os
# ===================================================================================
# Module: Genome-derived analysis
# Gene Functional Annotation Query
from Mainapp.Main_scripts.genome import GeneFunctionDataCurator
# ===================================================================================

# ===================================================================================
# Module: Transcriptome-derived analysis
# Gene Expression Query
from Mainapp.Main_scripts.transcriptome import GeneExpressionDataCurator
from Mainapp.R_visualization.transcriptome import gene_expression_query_plt_generator
# Differential Gene Expression Analysis
from Mainapp.Main_scripts.transcriptome import GeneDEDataCurator
from Mainapp.R_visualization.transcriptome import gene_differential_analysis_plt_generator
# GO / KEGG Enrichment
from Mainapp.Main_scripts.transcriptome import GoKeggDataCurator
from Mainapp.R_visualization.transcriptome import gene_gokegg_enrichment_plt_generator
# WGCNA
from Mainapp.Main_scripts.transcriptome import WGCNASampleGeneDataCurator
from Mainapp.R_visualization.transcriptome import gene_wgcna_section_deployer
# Sequence Obtain
from Mainapp.Main_scripts.transcriptome import GeneSequenceObtainer
# ===================================================================================

# ===================================================================================
# Module: Transcriptome-derived analysis
# Phylogenetic Analysis
from Mainapp.Main_scripts.genefamily import PhyloAnalysisHandler
from Mainapp.R_visualization.genefamily import genefam_phylogenetic_plt_generator
# Family Identification - hmmsearch
from Mainapp.Main_scripts.genefamily import GeneFamSearchHandler
# Family expression heatmap
from Mainapp.Main_scripts.genefamily import ExpressionHeatmapDataCurator
from Mainapp.R_visualization.genefamily import genefam_exp_heatmap_generator
# ===================================================================================

# ===================================================================================
# Module: TE analysis
# TE Expression Query
from Mainapp.Main_scripts.te import TEBySampleDataCurator, TEByFamilyDataCurator
from Mainapp.R_visualization.te import te_expression_query_plt_generator
from Mainapp.Main_scripts.te import TELocationDataCurator
# ===================================================================================

# ===================================================================================
# Module: Epigenome analysis
# ChIP peak annotation analysis
from Mainapp.R_visualization.epigenome import chip_peakanno_analysis_plt_generator
# ATAC peak location search
from Mainapp.Main_scripts.epigenome import ATACPeakDataCurator


# ===================================================================================


# Sample task to validate if celery backend is working
@shared_task
def testing_task1():
    time.sleep(10)
    pass


# ===================================================================================
# Module: Genome-derived analysis
@shared_task
def genome_genefunction_query_handler(query_species_name, gene_list, gokegg_sheet_name):
    query_annotation_list_full = []  # [[gene_name, description, GOterm, KEGG_kterm, KEGG_pathway], [..]]
    for indv_gene in gene_list:
        try:
            # Add functional annotation information
            # [gene_name, primary_transcript_name, description, GOterm, KEGG_kterm, KEGG_pathway]
            indv_annotation_list = GeneFunctionDataCurator.gene_functional_annot_builder(indv_gene, gokegg_sheet_name)
            if not len(indv_annotation_list):
                return RuntimeError
            # insert gene name before transcript name
            indv_annotation_list.insert(0, indv_gene)
            query_annotation_list_full.append(indv_annotation_list)
        except:
            return RuntimeError

    gene_function_query_return_list = [gene_list, query_annotation_list_full, query_species_name]
    return gene_function_query_return_list


# ===================================================================================
# Module: Transcriptome-derived analysis
@shared_task
def gene_expression_query_handler(query_species_name, exp_sheet_name, gene_list, expression_format):
    """

    :param query_species_name:
    :param exp_sheet_name:
    :param gene_list:
    :param expression_format:
    :return:
    """

    query_list_full = []  # [[gene_ID, expressed_samples, total_samples, max, min, median], [..]]
    for indv_gene in gene_list:

        # Curate expression dataframe
        exp_df = GeneExpressionDataCurator.gene_exp_df_builder(indv_gene, exp_sheet_name)

        # Add basic expression statistics (expressed samples, avg, median, etc.)
        # [gene_ID, expressed_samples, total_samples, max, min, median]
        query_list_full.append(
            GeneExpressionDataCurator.fund_info_obtainer(indv_gene, exp_df, expression_format)
        )

        # Generate plots
        gene_expression_query_plt_generator.overview_barplot(
            exp_df, indv_gene, expression_format, query_species_name
        )
        gene_expression_query_plt_generator.tissue_specific_barplot(
            exp_df, indv_gene, expression_format, query_species_name
        )
        # Generate gene-specific expression table (CSV format)
        exp_df.to_csv("Mainapp/static/Temp_R_html/" + indv_gene + "_values_" + expression_format + ".csv")

    # Returning full gene expression information list to views function
    # [query_list, query_list_full, query_format, query_species, query_annotation_list_full]
    gene_exp_query_return_list = [gene_list, query_list_full, expression_format, query_species_name]
    return gene_exp_query_return_list


@shared_task
def gene_differential_analysis_handler(query_species_name, exp_sheet_name,
                                       group1_samples_list, group1_name, group2_samples_list, group2_name,
                                       logfc_threshold, pvalue_threshold, heatmap_gene_count):
    """

    :param query_species_name:
    :param exp_sheet_name:
    :param group1_samples_list:
    :param group1_name: [list]
    :param group2_samples_list:
    :param group2_name: [list]
    :param logfc_threshold:
    :param pvalue_threshold:
    :param heatmap_gene_count:
    :return:
    """
    DE_data_raw, DE_group_list, DE_group_color_list = GeneDEDataCurator.gene_de_df_builder(
        exp_sheet_name, group1_samples_list, group1_name, group2_samples_list, group2_name
    )
    task_destination_folder = gene_differential_analysis_plt_generator.GeneDifferentialAnalysis(
        DE_data_raw, DE_group_list, DE_group_color_list,
        logfc_threshold, pvalue_threshold, heatmap_gene_count,
        group1_name, group2_name
    )
    # Returning full differential gene expression information list to views function
    gene_differential_analysis_return_list = [task_destination_folder, query_species_name,
                                              group1_name, group1_samples_list, group2_name, group2_samples_list,
                                              logfc_threshold, pvalue_threshold, heatmap_gene_count]
    return gene_differential_analysis_return_list


@shared_task
def gene_gokegg_enrichment_handler(gene_list, gokegg_sheet_name, enrich_format, species_filename):

    gokegg_df = GoKeggDataCurator.gene_gokegg_df_builder(gene_list, gokegg_sheet_name)

    plt_filename_id = gene_gokegg_enrichment_plt_generator.gokegg_enrich_plt(
        gokegg_df, enrich_format, species_filename
    )
    gene_gokegg_enrichment_return_list = [plt_filename_id]
    return gene_gokegg_enrichment_return_list


@shared_task
def gene_wgcna_handler(args_dict, section_number):
    """

    :param args_dict:
    :param section_number:
    :return:
    """
    # try:
    if section_number == "1":

        # Generate unique folder name and create the folder
        dirname = "wgcna" + str(random.randint(int(1e8), int(1e9) - 1))
        os.makedirs("Mainapp/static/Temp_R_wgcna/" + dirname, exist_ok=True)

        # exp_sheet_name, gene_list, sample_list, expression_format
        sample_exp_df = WGCNASampleGeneDataCurator.wgcna_init_df_builder(
            args_dict["query_table"],
            args_dict["gene_list_final"],
            args_dict["accession_list_final"],
            args_dict["query_format"],
        )
        sampleGeneData = "Mainapp/static/Temp_R_wgcna/" + dirname + "/sample_expression_data.csv"
        sample_exp_df.to_csv(sampleGeneData)

        # sampleGeneData, RcCutoff, samplePerc, GeneNum, cutmethod, rscut, datatype, anamethod, dirname
        dirname_return = gene_wgcna_section_deployer.wgcna_section1_deployer(
            sampleGeneData,
            args_dict["exp_cut_threshold"],
            args_dict["sample_cut_threshold"],
            args_dict["soft_power_cutoff"],
            args_dict["query_format"],
            args_dict["transformation_method"],
            dirname,
        )

        return [dirname_return]
    elif section_number == "2":
        try:
            dirname_return = gene_wgcna_section_deployer.wgcna_section2_deployer(
                args_dict["query_sftpower"],
                args_dict["min_module_size"],
                args_dict["tree_cut_threshold"],
                args_dict["query_project_id"],
            )
        except FileNotFoundError:
            return [102]  # Project directory not found or previous section not successfully completed

        return [dirname_return]
    elif section_number == "3":
        dirname_return = gene_wgcna_section_deployer.wgcna_section3_deployer(
            args_dict["query_project_id"],
        )

        return [dirname_return]
    # except:
    #    return RuntimeError


@shared_task
def gene_sequence_obtain_handler(query_species_name_filename, gene_list, query_species_name, search_longest):
    try:
        # return a list containing all the fasta entries, with each element representing a single fasta line.
        if search_longest == "T":
            sequence_fasta_filename = GeneSequenceObtainer.gene_sequence_query(gene_list, query_species_name_filename,
                                                                               query_species_name, True, True)
            # The last two parameters: write_fasta=True, primary_seq_req=True
        else:
            sequence_fasta_filename = GeneSequenceObtainer.gene_sequence_query(gene_list, query_species_name_filename,
                                                                               query_species_name, True, False)
            # The last two parameters: write_fasta=True, primary_seq_req=False
        # This is the POST method for sequence query,
        # so all available sequences should be displayed (primary_seq_req=False)
    except:
        return RuntimeError
    gene_sequence_obtain_return_list = [sequence_fasta_filename]
    return gene_sequence_obtain_return_list


# ===================================================================================

# ===================================================================================
# Gene family analysis
@shared_task
def genefam_phylogenetic_analysis_handler(query_species_name, query_sequences,
                                          msa_method, pairwise_dist_method, tree_layout_method, bootstrap_rep_num,
                                          is_protein_sequence_raw="dna"):
    """

    :param query_species_name:
    :param query_sequences:
    :param msa_method:
    :param pairwise_dist_method:
    :param tree_layout_method:
    :param bootstrap_rep_num:
    :param is_protein_sequence_raw: "protein" / "dna"
    :return:
    """
    try:
        # Validate sequence types
        if is_protein_sequence_raw == "protein":
            is_protein_sequence = True
        elif is_protein_sequence_raw == "dna":
            is_protein_sequence = False
        else:
            return RuntimeError

        sequence_fasta_string = PhyloAnalysisHandler.genefam_data_validation(query_sequences, query_species_name,
                                                                             is_protein_sequence)
        result_folder_name = genefam_phylogenetic_plt_generator.genefam_phylogenetic_analysis(
            sequence_fasta_string,
            msa_method,
            pairwise_dist_method,
            tree_layout_method,
            bootstrap_rep_num,
            is_protein_sequence,
        )
    except:
        return RuntimeError
    genefam_phylo_analysis_return_list = [result_folder_name]
    return genefam_phylo_analysis_return_list


@shared_task
def genefam_identification_hmmsearch_handler(species_filename, hmm_query_filename,
                                             seq_evalue_threshold, dom_evalue_threshold):
    try:
        result_folder_name = GeneFamSearchHandler.hmmsearch_deployer(species_filename, hmm_query_filename,
                                                                     seq_evalue_threshold, dom_evalue_threshold)
        return [result_folder_name]
    except:
        return RuntimeError


@shared_task
def genefam_expression_heatmap_handler(query_sheet_name, query_gene_list, query_sample_list, query_format):
    try:
        try:
            tmp_df = ExpressionHeatmapDataCurator.gene_expheatmap_df_builder(
                query_sheet_name, query_gene_list, query_sample_list, query_format
            )
        except LookupError:
            return [102]  # Database query error
        result_folder_name = genefam_exp_heatmap_generator.genefam_heatmap_generator(
            tmp_df, query_format
        )
        return [result_folder_name, query_format]
    except:
        return [100]


# ===================================================================================

# ===================================================================================
# Transposable Element (TE) analysis
@shared_task
def te_expression_query_handler(query_species_name, exp_sheet_name, query_list, query_method, query_format):
    """

    :param query_format: TPM / FPKM
    :param query_species_name:
    :param exp_sheet_name: ex. coix_te
    :param query_list: sample accession list or te family ID list
    :param query_method: sample / teid
    :return:
    """
    final_return_list = [[]]

    if query_method == "sample":
        for indv_query in query_list:
            tmp_exp_df = TEBySampleDataCurator.TE_exp_df_builder_sample(
                indv_query, exp_sheet_name
            )
            sample_id = te_expression_query_plt_generator.TE_bysample_plt(
                tmp_exp_df, indv_query, query_format
            )
            final_return_list[0].append(sample_id)

    elif query_method == "teid":
        for indv_query in query_list:
            te_class, te_class_group, tmp_exp_df = TEByFamilyDataCurator.TE_exp_df_builder_family(
                indv_query, exp_sheet_name
            )
            result_list = te_expression_query_plt_generator.TE_byfamily_plt(
                tmp_exp_df, indv_query, query_format
            )
            # Information regarding a single TE family ID query
            # [min_tpm/fpkm, mean_tpm/fpkm, max_tpm/fpkm, sample_number, sample_number_expressed, TE_id]
            final_return_list[0].append(result_list)
    else:
        return RuntimeError
    # [0], result_list, nested list containing information regarding each query
    final_return_list.append(query_species_name)  # [1], ex. zea
    final_return_list.append(query_method)  # [2], sample / teid
    final_return_list.append(query_format)  # [3], FPKM / TPM
    return final_return_list


@shared_task
def te_insertion_location_search_handler(query_species, query_flankingregion_length,
                                         query_table_tecoord, query_table_gene, mainquery, querytype):
    """

    :param query_table_gene:
    :param query_table_tecoord:
    :param query_species:
    :param query_flankingregion_length:
    :param mainquery: geneid -> list, containing all gene IDs / range -> [Seq]:[start]..[end]
    :param querytype: geneid / range
    :return:
    """
    try:
        task_return_list = TELocationDataCurator.te_location_df_builder(
            query_species, query_flankingregion_length, query_table_tecoord, query_table_gene, mainquery, querytype
        )
        if task_return_list[0] == 101:
            # Range query, input not in the correct format
            return 101
        else:
            # task return list:
            # if task identified by range, then the list comprises [[peak_entry1], [peak_entry2],...]
            # if task identified by gene id, then the list comprises [[gene_entries], [peak_entries]],
            #   and each [entries] is a nested list containing all the individual entries.
            return [query_species, querytype, task_return_list]

    except:
        return RuntimeError


# ===================================================================================

# ===================================================================================
# Module: Epigenome analysis
@shared_task
def chip_peakanno_analysis_handler(sample_id, tssRegion, species,
                                   ignore_1st_exon, ignore_1st_intron, ignore_downstream, ignore_promoter_subcategory):
    """

    :param sample_id:
    :param tssRegion:
    :param species:
    :param ignore_1st_exon: T / F
    :param ignore_1st_intron: T / F
    :param ignore_downstream: T / F
    :param ignore_promoter_subcategory: T / F
    :return:
    """
    try:
        # Generate plots and files
        output_folder_name = chip_peakanno_analysis_plt_generator.SingleChIPAnalysis(
            sample_id, tssRegion, species, ignore_1st_exon='F', ignore_1st_intron='F', ignore_downstream='F',
            ignore_promoter_subcategory='F'
        )
        if output_folder_name == 102:
            return 102  # narrowPeak file not found
    except:
        return RuntimeError

    ChIP_analysis_return_list = [sample_id, tssRegion, species,
                                 ignore_1st_exon, ignore_1st_intron, ignore_downstream, ignore_promoter_subcategory,
                                 output_folder_name]
    return ChIP_analysis_return_list


@shared_task
def atac_peak_location_search_handler(query_species, query_flankingregion_length,
                                      query_table_narrowpeak, query_table_gene, mainquery, querytype):
    """

    :param query_table_gene:
    :param query_table_narrowpeak:
    :param query_species:
    :param query_flankingregion_length:
    :param mainquery: geneid -> list, containing all gene IDs / range -> [Seq]:[start]..[end]
    :param querytype: geneid / range
    :return:
    """
    try:
        task_return_list = ATACPeakDataCurator.atac_peak_df_builder(
            query_species, query_flankingregion_length, query_table_narrowpeak, query_table_gene, mainquery, querytype
        )
        if task_return_list[0] == 101:
            # Range query, input not in the correct format
            return 101
        else:
            # task return list:
            # if task identified by range, then the list comprises [[peak_entry1], [peak_entry2],...]
            # if task identified by gene id, then the list comprises [[gene_entries], [peak_entries]],
            #   and each [entries] is a nested list containing all the individual entries.
            return [query_species, querytype, task_return_list]

    except:
        return RuntimeError
