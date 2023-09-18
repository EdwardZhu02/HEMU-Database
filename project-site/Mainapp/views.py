#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main
@File    ：views.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2022/8/27 18:00
@IDE     ：PyCharm
-----------------------------
Description: The main views function for HEMU database.
Input:
Output:

Environment requirement (if any): HEMUdb_new (Remote Server, Ubuntu 20.04)
- Dependencies: Django(v3.1), Project files, etc.
-----------------------------
Notes:

-----------------------------
File Revisions:
    2022/8/27: Version 1 - Creation
"""
import json
import os
from django.shortcuts import render, HttpResponse, redirect

# Server Configurations
from Mainapp.Main_scripts import MainConfiguration
from Mainapp.Main_scripts.site_manage import TmpFileCleaner

# Jbrowse2 Module
from Mainapp.Main_scripts.jbrowse import jbrowse_species_catalog

# Celery applications
import Mainapp.tasks as celery_tasks  # Async task list
from celery.result import AsyncResult  # Interface for task status querying


# ==============ASYNC FRAMEWORK==============
def get_async_task_progress(request):
    """
    Global view function for querying task progeess.
    :param request: AJAX request, GET, ?taskid=[celery_task_id]
    :return: json, {'state': 'PENDING'} / {'state': 'SUCCESS'} / {'state': 'FAILURE'}
    """
    if request.is_ajax():
        taskid = request.GET.get('taskid')
        if taskid:
            result = AsyncResult(taskid)
            # Another validation method, print(result.successful())
            response_data = {'state': result.state, }
            return HttpResponse(json.dumps(response_data), content_type='application/json')
    return render(request, 'error.html', {'error_message': 'INVALID_REQUEST'})


def render_asyc_onhold(request):
    task_id = request.GET.get('taskid')
    redirect_url = request.GET.get('dir')
    is_var_transfer = request.GET.get('var')
    if task_id and redirect_url:
        if is_var_transfer == "1":
            # Transfer vars in the redirect url
            return render(request, 'async_task_onhold_vartransfer.html',
                          {'task_id': task_id, 'redirect_url': redirect_url})
        # No var transfer in the redirect url
        return render(request, 'async_task_onhold.html',
                      {'task_id': task_id, 'redirect_url': redirect_url})
    return render(request, 'error.html', {'error_message': 'INVALID_URL_QUERY'})


def async_test(request):
    """
    Testing task deployer, submit a 10-second time consuming task to celery queue.
    :param request:
    :return:
    """
    if request.method == "GET":
        if request.GET.get('task') == "task1":
            async_result = celery_tasks.testing_task1.delay()
            return render(request, 'async_test/test_progressbar.html', {'task_id': async_result.task_id})
        elif request.GET.get('success') == "1":
            return render(request, 'error.html', {'error_message': 'TASK_IS_COMPLETED'})

        return render(request, 'async_test/test_progressbar.html')


# ==============HOME PAGE AND STATIC PAGES==============
# HOME PAGE RENDERER
def init_scr(request):
    return render(request, 'home.html')


# SITE MANAGER VALIDATOR & RENDERER
def site_manage(request):
    if request.method == "GET":
        # Common GET, request main page
        return render(request, 'site_manage/site_manage_auth.html')
    elif request.method == "POST":

        if request.POST.get('clear_tmp_files') == "1":
            identity_number_deleted = int(TmpFileCleaner.clean_tmp_files())
            tmp_message = "Cleaning completed, removed %d identities." % identity_number_deleted
            return render(request, 'site_manage/site_manage_dashboard.html',
                          {'clear_tmp_files_message': tmp_message})

        # Authentication
        if str(request.POST.get('password')) == MainConfiguration.query_admin_pwd():
            return render(request, 'site_manage/site_manage_dashboard.html')
        else:
            return render(request, 'site_manage/site_manage_auth.html',
                          {'error_message': 'Password Incorrect!'})


# AUX DATA RENDERER
def aux_data_render(request, aux_data_name):
    return render(request, 'user_guide.html')


# ==============MODULE1: GENOME AND GENE INFORMATION==============
# 1.1 GENE FUNCTION QUERY
def genome_gene_function_query_async(request):
    if request.method == 'GET':
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return [gene_list, query_annotation_list_full, query_species_name]
                    # query_annotation_list_full = \
                    # [[gene_name, primary_transcript_name, description, GOterm, KEGG_kterm, KEGG_pathway], [..]]
                    context = {
                        'query_list': task_return[0],
                        'query_annotation_list_full': task_return[1],
                        'query_species': task_return[2],
                        'last_query': ";".join(task_return[0]),
                        'flag_showtable': 1,
                    }
                    return render(request, 'genome/gene_annot_search_main.html', context)

                # task_return == None means no database queue hit, probably due to faulty input.
                return render(request, 'genome/gene_annot_search_main.html',
                              {'error_message': 'Invalid entry. Please check gene nomenclature and submit the query '
                                                'again.'})

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE'})
        # Render default view
        return render(request, 'genome/gene_annot_search_main.html')

    elif request.method == "POST":
        # Task submission
        # Obtain query from frontend form submission
        query_list = request.POST.get('main_query').split(';')  # Original query data
        query_species = request.POST.get('query_species')  # Species to query

        query_gokegg_table = MainConfiguration.query_tables(str(query_species) + "_gokegg")

        # Deploy task
        async_result = celery_tasks.genome_genefunction_query_handler.delay(
            query_species,
            query_list,
            query_gokegg_table,
        )
        # Return page containing task id, awaiting frontend recognition
        return render(request, 'genome/gene_annot_search_main.html', {'task_id': async_result.task_id})


# 1.2 GENE STRUCTURE VISUALIZER
def genome_gene_sturcture_shinyapprender(request):
    if request.method == "GET":
        return render(request, 'genome/gene_structure_visualization.html')


# 1.3 GENOME SYNTENY VIEWER SELECTOR
def genome_gene_synteny_render(request):
    if request.method == "GET":
        # Synteny session ID, gene/synteny/?id=B73v4_BTx623
        full_name = jbrowse_species_catalog.jbrowse_syntenysession_query(str(request.GET.get('id')))
        if full_name:  # Determine whether client is sending a query
            if not full_name == -1:  # Capture key error
                return render(request, 'Jbrowse/jbrowse_synteny_insession.html',  # Sending a valid query
                              {'species_full_name': full_name})
        return render(request, 'genome/gene_synteny_nav.html')  # Sending an invalid query or not sending a quer


# ==============MODULE2: TRANSCRIPTOME-DERIVED ANALYSIS==============
# 2.1 GENE EXPRESSION PROFILE QUERY HANDLER
def gene_expression_query_async(request):
    if request.method == 'GET':
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return [query_list, query_list_full, query_format, query_species]
                    context = {
                        'query_list': task_return[0],
                        'query_list_full': task_return[1],
                        # 'query_annotation_list_full': task_return[4],
                        'query_format': task_return[2],
                        'query_species': task_return[3],
                        'last_query': ";".join(task_return[0]),
                    }
                    return render(request, 'gene_expression/gene_search_main.html', context)

                # task_return == None means no database queue hit, probably due to faulty input.
                return render(request, 'gene_expression/gene_search_main.html',
                              {'error_message': 'Invalid entry. Please check gene nomenclature and submit the query '
                                                'again.'})

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE'})
        # Render default view
        return render(request, 'gene_expression/gene_search_main.html')

    elif request.method == "POST":
        # Task submission
        # Obtain query from frontend form submission
        query_list = request.POST.get('main_query').split(';')  # Original query data
        query_species = request.POST.get('query_species')  # Species to query
        query_format = request.POST.get('optionsRadios')  # FPKM / TPM

        # Obtain sql tables
        query_exp_table = MainConfiguration.query_tables(str(query_species) + "_exp")

        # Deploy task
        async_result = celery_tasks.gene_expression_query_handler.delay(
            query_species,
            query_exp_table,
            query_list,
            query_format,
        )
        # Return page containing task id, awaiting frontend recognition
        return render(request, 'gene_expression/gene_search_main.html', {'task_id': async_result.task_id})


# 2.2 SEQUENCE ACQUISITION: Gene, Transcript, Protein
def gene_sequence_obtain_async(request):
    if request.method == "GET":

        # Submit GET-format task, used for small amount of querying.
        query_species = request.GET.get("sp")
        query_raw = request.GET.get("gene")
        query_format = request.GET.get("format")
        if query_species and query_raw and query_format:
            if len(query_raw.split(";")) > 10:
                return render(request, 'error.html', {'error_message': 'QUERY_TOO_LONG'})

            # Obtain filenames for fasta file containing gene sequences
            species_filename = MainConfiguration.query_tables(query_species + "_filename")
            # Deploy task
            # "F" means do not only search the longest sequence (Ref. tasks)
            async_result = celery_tasks.gene_sequence_obtain_handler.delay(
                species_filename, query_raw.split(";"), query_format, "F"
            )
            # Return async task manager, awaiting frontend recognition
            return render(request, 'gene_expression/gene_sequence_acquire.html', {'task_id': async_result.task_id})

        # Process task return status
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    # [sequence_fasta_filename]
                    if not task_return[0].startswith("sequence"):
                        # Empty result list
                        # Nomenclature: "sequence" + str(random.randint(int(1e8), int(1e9) - 1)) + ".fasta"
                        return render(request, 'gene_expression/gene_sequence_acquire.html',
                                      {'error_message': 'Invalid entry. Please check gene nomenclature and submit the '
                                                        'query again.'})

                    # Read file content and render to frontend
                    try:
                        with open("Mainapp/static/Temp_R_genefam/" + str(task_return[0]), mode='r') as fasta_fh:
                            fasta_content = fasta_fh.read()

                    except FileNotFoundError:
                        return render(request, 'gene_expression/gene_sequence_acquire.html',
                                      {'error_message': 'Invalid entry. Please check gene nomenclature and submit the '
                                                        'query again.'})
                    context = {
                        'result_list_full': [fasta_content],
                        'fasta_filename': task_return[0],
                    }
                    return render(request, 'gene_expression/gene_sequence_acquire.html', context)
                else:
                    return render(request, 'error.html', {'error_message': 'EMPTY_TASK_RESPONSE'})

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE'})
        # Render default view
        return render(request, 'gene_expression/gene_sequence_acquire.html')

    elif request.method == "POST":

        # Redirect for BLAST instance
        seq_transfer = request.POST.get("sequence_raw")
        if seq_transfer:
            return render(request, 'BLAST/sequenceserver_display.html', {'seq_transfer': seq_transfer, })

        # Task submission
        # Obtain query from frontend form submission
        query_list = request.POST.get('main_query').split(';')  # Original query data
        query_species = request.POST.get('query_species')  # Species to query
        query_format = request.POST.get('optionsRadios')  # gene / transcript / protein
        search_longest = request.POST.get('search_longest')  # T / F
        if not search_longest == "T":
            search_longest = "F"

        # Obtain filenames for fasta file containing gene sequences
        species_filename = MainConfiguration.query_tables(query_species + "_filename")
        # Deploy task (write .fasta file to disk in order to prevent celery overflow)
        async_result = celery_tasks.gene_sequence_obtain_handler.delay(
            species_filename, query_list, query_format, search_longest
        )
        # Return async task manager, awaiting frontend recognition
        return render(request, 'gene_expression/gene_sequence_acquire.html', {'task_id': async_result.task_id})


# 2.3 DGE QUERY HANDLER
def gene_differential_analysis_async(request):
    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    # [task_destination_folder, query_species_name, group1_name, group1_samples_list,
                    # group2_name, group2_samples_list, logfc_threshold, pvalue_threshold, heatmap_gene_count]
                    context = {
                        'task_destination_folder': task_return[0],
                        'species_query': task_return[1],
                        'group1_name': task_return[2],
                        'group1_samples_list': task_return[3],
                        'group2_name': task_return[4],
                        'group2_samples_list': task_return[5],
                        'logfc_threshold': task_return[6],
                        'pvalue_threshold': task_return[7],
                        'heatmap_gene_count': task_return[8],
                    }
                    return render(request, 'gene_expression/gene_DE_report_display.html', context)

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view
        return render(request, 'gene_expression/gene_DE_main.html')

    elif request.method == "POST":
        # Task submission
        # Obtain query from frontend form submission
        query_species = request.POST.get('species_query')
        logfc_threshold = request.POST.get('logfc_threshold')
        pvalue_threshold = request.POST.get('pvalue_threshold')
        heatmap_gene_count = request.POST.get('heatmap_gene_count')

        group1_samples_list = request.POST.get('group1_samples').split(";")
        group1_name = request.POST.get('group1_name')
        group2_samples_list = request.POST.get('group2_samples').split(";")
        group2_name = request.POST.get('group2_name')

        # Obtain sql tables
        query_exp_table = MainConfiguration.query_tables(str(query_species) + "_exp")

        # Deploy task
        async_result = celery_tasks.gene_differential_analysis_handler.delay(
            query_species, query_exp_table,
            group1_samples_list, group1_name, group2_samples_list, group2_name,
            logfc_threshold, pvalue_threshold, heatmap_gene_count
        )
        # Return async task manager, awaiting frontend recognition
        return redirect("/HEMUdb/asynctasks/onhold?taskid=%s&dir=%s" %
                        (async_result.task_id, "/HEMUdb/gene/dge"))


# 2.3 DGE STATIC FILE LOAD HANDLER
def load_DE_staticfile(request, identifier_name, file_name, frame_height, frame_width):
    return render(request, 'static_html_display.html',
                  {
                      'identifier_name': identifier_name,
                      'file_name': file_name,
                      'frame_height': frame_height,
                      'frame_width': frame_width,
                  })


# 2.4 GO/KEGG ENRICHMENT
def gene_gokegg_enrichment_async(request):
    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    # [plt_filename_id]
                    context = {
                        'filename_id': task_return[0],
                    }
                    return render(request, 'gene_expression/gokegg_enrich_result.html', context)

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view
        return render(request, 'gene_expression/gokegg_enrich.html')

    elif request.method == "POST":
        # Task submission
        # Obtain query from frontend form submission
        species_query = request.POST.get("query_species")
        enrich_format = request.POST.get('optionsRadios')  # GO / KEGG
        # Separate query and remove blank entries
        gene_list_final = [indv_gene.rstrip("\r") for indv_gene in request.POST.get("query_gene_list").split("\n")]
        gene_list_final = [indv_entry for indv_entry in gene_list_final if indv_entry != '']

        # Obtain sql tables
        query_table = MainConfiguration.query_tables(str(species_query) + "_gokegg")
        # Obtain filename for accessing background annotation
        species_filename = MainConfiguration.query_tables(str(species_query) + "_filename")

        # Deploy task
        # params: gene_list, gokegg_sheet_name, enrich_format
        async_result = celery_tasks.gene_gokegg_enrichment_handler.delay(
            gene_list_final, query_table, enrich_format, species_filename
        )
        # Return async task manager, awaiting frontend recognition
        return redirect("/HEMUdb/asynctasks/onhold?taskid=%s&dir=%s" %
                        (async_result.task_id, "/HEMUdb/gene/gokegg"))


# 2.5 WGCNA ANALYSIS
def gene_wgcna_analysis_async(request):
    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            section_redirect = request.GET.get('section')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Section 1 completed, display sft information for user to assess
                    if section_redirect == "1":
                        context = {'result_folder_name': task_return[0], }
                        return render(request, 'gene_expression/wgcna_result_section1.html', context)
                    elif section_redirect == "2":
                        if task_return[0] == "102":
                            # Project directory not found or previous section not successfully completed
                            return render(request, 'error.html', {'error_message': 'PROJECT_ID_ERROR'})
                        context = {'result_folder_name': task_return[0], }
                        return render(request, 'gene_expression/wgcna_result_section2.html', context)
                    elif section_redirect == "3":
                        context = {'result_folder_name': task_return[0], }
                        return render(request, 'gene_expression/wgcna_result_section3.html', context)
                    else:
                        return render(request, 'error.html', {'error_message': 'REQUESTED_SECTION_NOT_AVAILABLE'})

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view (data submission)
        return render(request, 'gene_expression/wgcna_ana_dashboard.html')

    elif request.method == "POST":
        # Task submission
        # Create a new instance, obtain query from frontend form submission
        form_section_id = request.POST.get("form_section_id")

        if form_section_id == "1":
            query_species = request.POST.get("query_species")
            query_format = request.POST.get("query_exp_met")  # FPKM / TPM
            transformation_method = request.POST.get("transformation_method")  # rawFPKM / logFPKM
            exp_cut_threshold = request.POST.get("exp_cut_threshold")
            sample_cut_threshold = request.POST.get("sample_cut_threshold")
            soft_power_cutoff = request.POST.get("soft_power_cutoff")

            # Obtain gene list
            gene_list_final = [indv_gene.rstrip("\r") for indv_gene in request.POST.get("query_gene_list").split("\n")]
            gene_list_final = [indv_entry for indv_entry in gene_list_final if indv_entry != '']  # Remove blank entries
            # Obtain accession list
            accession_list_final = [indv_accession.rstrip("\r") for indv_accession in
                                    request.POST.get("query_accession_list").split("\n")]
            accession_list_final = [indv_entry for indv_entry in
                                    accession_list_final if indv_entry != '']  # Remove blank entries
            # Check gene list & accession list length
            if not 1 < len(gene_list_final) < 4001:
                return render(request, 'gene_expression/wgcna_ana_dashboard.html', {
                    'error_message': 'Invalid gene list length (1<N<4001), please re-submit the query.'
                })
            if not 1 < len(accession_list_final) < 61:
                return render(request, 'gene_expression/wgcna_ana_dashboard.html', {
                    'error_message': 'Invalid sample accession list length (1<N<61), please re-submit the query.'
                })

            # Obtain sql tables for expression data
            query_table = MainConfiguration.query_tables(str(query_species) + "_exp")
            # Pack function context
            task_context = {
                'query_table': query_table,
                'query_species': query_species,
                'query_format': query_format,
                'transformation_method': transformation_method,
                'exp_cut_threshold': exp_cut_threshold,
                'sample_cut_threshold': sample_cut_threshold,
                'soft_power_cutoff': soft_power_cutoff,
                'gene_list_final': gene_list_final,
                'accession_list_final': accession_list_final,
            }
            # Submit task, section 1
            async_result = celery_tasks.gene_wgcna_handler.delay(task_context, "1")
            return redirect("/HEMUdb/asynctasks/onhold?var=1&taskid=%s&dir=%s" %
                            (async_result.task_id, "/HEMUdb/gene/wgcna?section=1"))
        elif form_section_id == "2":
            query_project_id = request.POST.get("query_projectid")
            query_sftpower = request.POST.get("query_sftpower")
            min_module_size = request.POST.get("min_module_size")
            tree_cut_threshold = request.POST.get("tree_cut_threshold")

            # Pack function context
            task_context = {
                'query_project_id': query_project_id,
                'query_sftpower': query_sftpower,
                'min_module_size': min_module_size,
                'tree_cut_threshold': tree_cut_threshold,
            }
            # Submit task, section 2
            async_result = celery_tasks.gene_wgcna_handler.delay(task_context, "2")
            return redirect("/HEMUdb/asynctasks/onhold?var=1&taskid=%s&dir=%s" %
                            (async_result.task_id, "/HEMUdb/gene/wgcna?section=2"))
        elif form_section_id == "3":
            query_project_id = request.POST.get("query_projectid")
            sampletrait_tbl_file = request.FILES.get("query_sampletrait_tbl_file", None)

            # Filter file name
            if not sampletrait_tbl_file.name.endswith(".csv"):
                return render(request, 'error.html', {'error_message': 'FILE_FORMAT_ERROR'})

            # Detect folder and previous project files
            if not os.path.exists('Mainapp/static/Temp_R_wgcna/' + query_project_id):
                return render(request, 'error.html', {'error_message': 'PROJECT_ID_ERROR'})
            if not os.path.exists('Mainapp/static/Temp_R_wgcna/' + query_project_id + '/step1-vars.RData'):
                return render(request, 'error.html', {'error_message': 'PROJECT_TMPFILE_ERROR'})
            if not os.path.exists('Mainapp/static/Temp_R_wgcna/' + query_project_id + '/step1-vars.RData'):
                return render(request, 'error.html', {'error_message': 'PROJECT_TMPFILE_ERROR'})

            destination = open('Mainapp/static/Temp_R_wgcna/' + query_project_id + '/sample_trait_data.csv', 'wb+')
            for chunk in sampletrait_tbl_file.chunks():
                destination.write(chunk)
            destination.close()

            # Pack function context
            task_context = {
                'query_project_id': query_project_id,
            }
            # Submit task, section 2
            async_result = celery_tasks.gene_wgcna_handler.delay(task_context, "3")
            return redirect("/HEMUdb/asynctasks/onhold?var=1&taskid=%s&dir=%s" %
                            (async_result.task_id, "/HEMUdb/gene/wgcna?section=3"))
        elif form_section_id == "4":
            # Access results
            query_project_id = request.POST.get("query_projectid")
            # Detect folder and previous project files
            if not os.path.exists('Mainapp/static/Temp_R_wgcna/' + query_project_id):
                return render(request, 'error.html', {'error_message': 'PROJECT_ID_NOT_DETECTED'})
            if os.path.exists(
                    'Mainapp/static/Temp_R_wgcna/' + query_project_id + '/step3-module-trait-correlation.csv'):
                # Section 3 completed
                return render(request, 'gene_expression/wgcna_result_section3.html',
                              {'result_folder_name': query_project_id, })
            if os.path.exists('Mainapp/static/Temp_R_wgcna/' + query_project_id + '/step2-vars.RData'):
                # Section 2 completed
                return render(request, 'gene_expression/wgcna_result_section2.html',
                              {'result_folder_name': query_project_id, })
            if os.path.exists('Mainapp/static/Temp_R_wgcna/' + query_project_id + '/step1-vars.RData'):
                # Section 1 completed
                return render(request, 'gene_expression/wgcna_result_section1.html',
                              {'result_folder_name': query_project_id, })
            return render(request, 'error.html', {'error_message': 'SECTION_DAMAGED'})
        else:
            return render(request, 'error.html', {'error_message': 'FORM_ID_ERROR'})


# ============MODULE3: GENE FAMILY ANALYSIS=============
# 3.1 GENE PHYLOGENETIC ANALYSIS
def genefam_phylogenetic_analysis_async(request):
    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    # [plt_filename_id]
                    context = {
                        'result_folder_name': task_return[0],
                    }
                    return render(request, 'gene_family/gene_phylogenic_analysis_result.html', context)

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view
        return render(request, "gene_family/gene_phylogenetic_analysis.html")

    elif request.method == "POST":
        # Task submission
        # Obtain query from frontend form submission
        species_query = request.POST.get("query_species")  # by HTML tag name
        sequence_query = request.POST.get("query_sequence")
        msa_method = request.POST.get("query_msa_method")
        pairwise_dist_method = request.POST.get("query_pairwise_dist_method")
        tree_layout_method = request.POST.get("query_tree_layout_method")
        bootstrap_rep_num = request.POST.get("query_bootstrap_rep_num")
        is_protein_sequence_raw = request.POST.get("query_seqtype")  # protein / dna

        try:
            # Validate if variable bootstrap_rep_num is a number
            int(bootstrap_rep_num) / 2
        except ValueError:
            return render(request, 'error.html', {'error_message': 'INVALID_QUERY'})

        # Deploy task
        # params: gene_list, gokegg_sheet_name, enrich_format
        async_result = celery_tasks.genefam_phylogenetic_analysis_handler.delay(
            species_query, sequence_query,
            msa_method, pairwise_dist_method, tree_layout_method, bootstrap_rep_num,
            is_protein_sequence_raw,
        )
        # Return async task manager, awaiting frontend recognition
        return redirect("/HEMUdb/asynctasks/onhold?taskid=%s&dir=%s" %
                        (async_result.task_id, "/HEMUdb/genefam/phylo"))


# 3.2 HMM-BASED GENE FAMILY IDENTIFICATION
def genefam_identification_hmm_async(request):
    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    # [plt_filename_id]
                    if task_return[0] == 102:
                        return render(request, 'gene_family/gene_hmmsearch_analysis.html',
                                      {'error_message': 'HMM file not found on the server. Please submit again.'})
                    if task_return[0] == 103:
                        return render(request, 'error.html', {'error_message': 'HMMSEARCH_RUNTIME_ERROR'})
                    context = {
                        'output_file_name': task_return[0],
                    }
                    return render(request, 'gene_family/gene_fam_identification_result.html', context)

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view
        return render(request, 'gene_family/gene_hmmsearch_analysis.html')

    elif request.method == "POST":
        # Task submission
        # Obtain query from frontend form submission
        species_query = request.POST.get("query_species")  # by HTML tag name
        hmm_query = request.POST.get("hmm_query")
        sequence_evalue_threshold = request.POST.get("sequence_evalue")
        domain_evalue_threshold = request.POST.get("domain_evalue")
        try:
            # Validate if variable sequence_evalue_threshold and domain_evalue_threshold are numbers
            int(float(sequence_evalue_threshold)) / 2
            int(float(domain_evalue_threshold)) / 2
            # Make sure both numbers are positive
            if int(float(sequence_evalue_threshold)) < 1 or int(float(domain_evalue_threshold)) < 1:
                raise ValueError
        except ValueError:
            return render(request, 'error.html', {'error_message': 'INVALID_QUERY'})

        # Obtain filename for accessing background annotation
        species_filename = MainConfiguration.query_tables(str(species_query) + "_filename")

        # Deploy task
        # params: species_filename, hmm_query_filename, seq_evalue_threshold, dom_evalue_threshold
        async_result = celery_tasks.genefam_identification_hmmsearch_handler.delay(
            species_filename, hmm_query, int(float(sequence_evalue_threshold)), int(float(domain_evalue_threshold))
        )
        # Return async task manager, awaiting frontend recognition
        return redirect("/HEMUdb/asynctasks/onhold?taskid=%s&dir=%s" %
                        (async_result.task_id, "/HEMUdb/genefam/hmm"))


# 3.3 HMM-BASED GENE FAMILY IDENTIFICATION
def genefam_identification_blastp_async(request):
    if request.method == "GET":
        return render(request, 'gene_family/gene_blastpsearch_analysis.html')


# 3.4 GENE FAMILY EXPRESSION HEATMAP GENERATOR
def genefam_expheatmapgen_async(request):
    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    # [plt_filename_id]
                    if len(task_return) < 2:
                        if task_return[0] == 102:  # Database query error
                            return render(request, 'gene_family/gene_expheatmap_gen_main.html',
                                          {'error_message': 'DATABASE QUERY ERROR: No database hit in some of the '
                                                            'entries. Please check input and submit the query again.'})
                        else:  # Database query error
                            return render(request, 'gene_family/gene_expheatmap_gen_main.html',
                                          {'error_message': 'GENERAL ERROR: No database hit in some of the '
                                                            'entries. Please check input and submit the query again.'})

                    context = {
                        'result_folder_name': task_return[0],
                        'query_format': task_return[1],
                    }
                    return render(request, 'gene_family/gene_expheatmap_gen_result.html', context)

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view
        return render(request, 'gene_family/gene_expheatmap_gen_main.html')

    elif request.method == "POST":
        # Task submission

        # Obtain query from frontend form submission
        query_species = request.POST.get("query_species")
        query_format = request.POST.get("query_exp_met")  # FPKM / TPM

        # Obtain gene list
        gene_list_final = [indv_gene.rstrip("\r") for indv_gene in request.POST.get("query_gene_list").split("\n")]
        gene_list_final = [indv_entry for indv_entry in gene_list_final if indv_entry != '']  # Remove blank entries
        # Obtain accession list
        accession_list_final = [indv_accession.rstrip("\r") for indv_accession in
                                request.POST.get("query_accession_list").split("\n")]
        accession_list_final = [indv_entry for indv_entry in
                                accession_list_final if indv_entry != '']  # Remove blank entries
        # Check gene list & accession list length
        if not 1 < len(gene_list_final) < 101:
            return render(request, 'gene_expression/wgcna_ana_dashboard.html', {
                'error_message': 'Invalid gene list length (1<N<101), please re-submit the query.'
            })
        if not 1 < len(accession_list_final) < 101:
            return render(request, 'gene_expression/wgcna_ana_dashboard.html', {
                'error_message': 'Invalid sample accession list length (1<N<101), please re-submit the query.'
            })

        # Obtain sql tables for expression data
        query_table = MainConfiguration.query_tables(str(query_species) + "_exp")

        # Deploy task
        # params: query_sheet_name, query_gene_list, query_sample_list, query_format
        async_result = celery_tasks.genefam_expression_heatmap_handler.delay(
            query_table, gene_list_final, accession_list_final, query_format
        )

        # Return async task manager, awaiting frontend recognition
        return redirect("/HEMUdb/asynctasks/onhold?taskid=%s&dir=%s" %
                        (async_result.task_id, "/HEMUdb/genefam/htmap"))


# ==============MODULE4: TE-DERIVED ANALYSIS==============
# 4.1 TE EXPRESSION PROFILE SEARCH
def te_expression_query_async(request):
    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    # [[query1, query2, ...], query_species_name, query_method, query_format]
                    if task_return[2] == "sample":
                        context = {
                            'query_sample_result_list': task_return[0],
                            'query_species_name': task_return[1],
                            'query_method': task_return[2],
                            'query_format': task_return[3],
                        }
                    elif task_return[2] == "teid":
                        context = {
                            'query_teid_result_list': task_return[0],
                            # nested list containing information regarding each query
                            # example of a single element:
                            # [min_tpm/fpkm, mean_tpm/fpkm, max_tpm/fpkm, sample_number, sample_number_expressed, TE_id]

                            'query_species_name': task_return[1],
                            'query_method': task_return[2],
                            'query_format': task_return[3],
                        }
                    else:
                        context = {'error_message': 'Error: empty query method'}
                    return render(request, 'te/te_search_main.html', context)

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view
        return render(request, 'te/te_search_main.html')

    elif request.method == "POST":
        # Task submission

        # Obtain query from frontend form submission
        form_section_id = request.POST.get('form_section_id')  # 1-sample; 2-TEid
        query_species = request.POST.get('species_query')
        query_format = request.POST.get('optionsRadios')  # FPKM / TPM

        # Obtain sql tables and any static files (if needed)
        query_table = MainConfiguration.query_tables(str(query_species) + "_te")

        if form_section_id == "1":  # by-sample-accession search
            # Obtain sample list
            sample_list_final = [indv_sample for indv_sample in request.POST.get('te_query_accession').split(";")]
            sample_list_final = [indv_sample for indv_sample in sample_list_final if indv_sample != '']

            # Deploy task
            async_result = celery_tasks.te_expression_query_handler.delay(
                query_species,
                query_table,
                sample_list_final,
                "sample",
                query_format,
            )

        elif form_section_id == "2":  # by-te family ID search
            # Obtain TE id list
            teid_list_final = [indv_te for indv_te in request.POST.get('te_query_teid').split(";")]
            teid_list_final = [indv_te for indv_te in teid_list_final if indv_te != '']

            # Deploy task
            async_result = celery_tasks.te_expression_query_handler.delay(
                query_species,
                query_table,
                teid_list_final,
                "teid",
                query_format,
            )
        else:
            return render(request, 'error.html', {'error_message': 'EXCEPTION'})

        # Return async task manager, awaiting frontend recognition
        return render(request, 'te/te_search_main.html', {'task_id': async_result.task_id})


# 4.2 TE INSERTION LOCATION PROFILE SEARCH
def te_insertion_shinyapprender(request):
    if request.method == "GET":
        return render(request, 'te/te_insertion_main.html')


def te_location_search_async(request):
    # return HttpResponse("Feature under maintenance.")

    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    if task_return == 101:
                        # Range query, input not in the correct format
                        return render(request, 'te/te_insertion_loc_search.html', {
                            'error_message': '[QUERY ERROR] Query invalid. Please check and submit again.'
                        })
                    # [query_species, querytype, task_return_list]
                    # if task identified by range, then the list comprises [[peak_entry1], [peak_entry2],...]
                    # if task identified by gene id, then the list comprises [[gene_entries], [peak_entries]],
                    #   and each [entries] is a nested list containing all the individual entries.
                    if task_return[1] == "range":
                        context = {
                            'query_species': task_return[0],
                            'query_type': task_return[1],
                            'peak_list_full': task_return[2],
                            'flag_showpeaktable': 1,
                        }
                        return render(request, 'te/te_insertion_loc_search.html', context)
                    elif task_return[1] == "geneid":
                        context = {
                            'query_species': task_return[0],
                            'query_type': task_return[1],
                            'gene_list_full': task_return[2][0],
                            'peak_list_full': task_return[2][1],
                            'flag_showgenetable': 1,
                            'flag_showpeaktable': 1,
                        }
                        return render(request, 'te/te_insertion_loc_search.html', context)

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view
        return render(request, 'te/te_insertion_loc_search.html')

    elif request.method == "POST":
        # Task submission

        # Obtain query from frontend form submission
        query_species = request.POST.get("query_species")
        query_flankingregion_length = request.POST.get("query_flanking_length")
        query_gene_id = request.POST.get("gene_query")
        query_range = request.POST.get("range_query")

        # Validate query, make sure only one condition is submitted at a time
        if query_gene_id and query_range:
            return render(request, 'te/te_insertion_loc_search.html', {
                'error_message': '[QUERY ERROR] Please do not submit gene ID and range query at the same time.'
            })
        if not query_gene_id and not query_range:
            return render(request, 'te/te_insertion_loc_search.html', {
                'error_message': '[QUERY ERROR] Please submit a valid query.'
            })

        # Obtain sql tables and any static files (if needed)
        query_table_tecoord = MainConfiguration.query_tables(str(query_species) + "_tecoord")
        query_table_gene = MainConfiguration.query_tables(str(query_species) + "_genecoord")

        if query_gene_id:
            # Process query, split into unique elements
            geneid_list_final = [indv_gene for indv_gene in query_gene_id.split(";")]
            geneid_list_final = [indv_gene for indv_gene in geneid_list_final if indv_gene != '']
            # Deploy task
            async_result = celery_tasks.te_insertion_location_search_handler.delay(
                query_species, query_flankingregion_length, query_table_tecoord, query_table_gene, geneid_list_final,
                'geneid'
            )
        elif query_range:
            # Deploy task
            async_result = celery_tasks.te_insertion_location_search_handler.delay(
                query_species, query_flankingregion_length, query_table_tecoord, query_table_gene, query_range,
                'range'
            )
        else:
            return render(request, 'error.html', {'error_message': 'INVALID_QUERY'})

        # Return async task manager, awaiting frontend recognition
        return render(request, 'te/te_insertion_loc_search.html', {'task_id': async_result.task_id})


# ==============MODULE5: EPIGENOME ANALYSIS==============
# ChIP PEAK ANNOTATION ANALYSIS (SINGLE-SAMPLE)
def chip_peakanno_analysis_async(request):
    if request.method == 'GET':
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:

                    if task_return == 102:
                        # File not found error
                        return render(request, 'epigenome/chip_peakanno_main.html',
                                      {
                                          'error_message': 'Sample ID not found for peak file. Please check '
                                                           'nomenclature and species type, then submit the query again.'
                                      })
                    # Task completed, return list
                    # [sample_id, tssRegion, species, ignore_1st_exon, ignore_1st_intron,
                    # ignore_downstream,ignore_promoter_subcategory, output_folder_name]
                    context = {
                        'sample_id': task_return[0],
                        'tssRegion': task_return[1],
                        'query_species': task_return[2],
                        'ignore_1st_exon': task_return[3],
                        'ignore_1st_intron': task_return[4],
                        'ignore_downstream': task_return[5],
                        'ignore_promoter_subcategory': task_return[6],
                        'result_folder_name': task_return[7]
                    }
                    return render(request, 'epigenome/chip_peakanno_result.html', context)

                # task_return == None means no database queue hit, probably due to faulty input.
                return render(request, 'epigenome/chip_peakanno_main.html',
                              {'error_message': 'Invalid entry. Please check gene nomenclature and submit the query '
                                                'again.'})
        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE'})
        # Render default view
        return render(request, 'epigenome/chip_peakanno_main.html')

    elif request.method == "POST":
        # Task submission
        # Obtain query from frontend form submission
        query_list = request.POST.get('main_query')  # Original query data
        query_species = request.POST.get('query_species')  # Species to query
        query_tss = request.POST.get('query_tss_length')  # Tss Region setting
        ignore_1st_exon = request.POST.get('ignore_1st_exon')
        if not ignore_1st_exon == "T":
            ignore_1st_exon = "F"
        ignore_1st_intron = request.POST.get('ignore_1st_intron')
        if not ignore_1st_intron == "T":
            ignore_1st_intron = "F"
        ignore_downstream = request.POST.get('ignore_downstream')
        if not ignore_downstream == "T":
            ignore_downstream = "F"
        ignore_promoter_subcategory = request.POST.get('ignore_promoter_subcategory')
        if not ignore_promoter_subcategory == "T":
            ignore_promoter_subcategory = "F"

        # Deploy task
        async_result = celery_tasks.chip_peakanno_analysis_handler.delay(
            query_list, query_tss, query_species,
            ignore_1st_exon, ignore_1st_intron, ignore_downstream, ignore_promoter_subcategory
        )
        # Return page containing task id, awaiting frontend recognition
        return redirect("/HEMUdb/asynctasks/onhold?taskid=%s&dir=%s" %
                        (async_result.task_id, "/HEMUdb/epigenome/chip"))


def atac_peak_location_search_async(request):
    if request.method == "GET":
        taskstate = request.GET.get("success")
        if taskstate == "1":
            # async task completed
            taskid = request.GET.get('taskid')
            if taskid:
                task_result = AsyncResult(taskid)
                if not task_result.state == "SUCCESS":
                    # Detect illegal task submissions
                    return render(request, 'error.html', {'error_message': 'REQUESTED_TASK_NOT_COMPLETED'})
                task_return = task_result.get()
                if task_return:
                    # Task completed, return list
                    if task_return == 101:
                        # Range query, input not in the correct format
                        return render(request, 'epigenome/atac_peakloc_search.html', {
                            'error_message': '[QUERY ERROR] Query invalid. Please check and submit again.'
                        })
                    # [query_species, querytype, task_return_list]
                    # if task identified by range, then the list comprises [[peak_entry1], [peak_entry2],...]
                    # if task identified by gene id, then the list comprises [[gene_entries], [peak_entries]],
                    #   and each [entries] is a nested list containing all the individual entries.
                    if task_return[1] == "range":
                        context = {
                            'query_species': task_return[0],
                            'query_type': task_return[1],
                            'peak_list_full': task_return[2],
                            'flag_showpeaktable': 1,
                        }
                        return render(request, 'epigenome/atac_peakloc_search.html', context)
                    elif task_return[1] == "geneid":
                        context = {
                            'query_species': task_return[0],
                            'query_type': task_return[1],
                            'gene_list_full': task_return[2][0],
                            'peak_list_full': task_return[2][1],
                            'flag_showgenetable': 1,
                            'flag_showpeaktable': 1,
                        }
                        return render(request, 'epigenome/atac_peakloc_search.html', context)

        elif taskstate == "0":
            # async task failed
            return render(request, 'error.html', {'error_message': 'CELERY_ASYNC_TASK_FAILURE_EXCEPTION'})
        # Render default view
        return render(request, 'epigenome/atac_peakloc_search.html')

    elif request.method == "POST":
        # Task submission

        # Obtain query from frontend form submission
        query_species = request.POST.get("query_species")
        query_flankingregion_length = request.POST.get("query_flanking_length")
        query_gene_id = request.POST.get("gene_query")
        query_range = request.POST.get("range_query")

        # Validate query, make sure only one condition is submitted at a time
        if query_gene_id and query_range:
            return render(request, 'epigenome/atac_peakloc_search.html', {
                'error_message': '[QUERY ERROR] Please do not submit gene ID and range query at the same time.'
            })
        if not query_gene_id and not query_range:
            return render(request, 'epigenome/atac_peakloc_search.html', {
                'error_message': '[QUERY ERROR] Please submit a valid query.'
            })

        # Obtain sql tables and any static files (if needed)
        query_table_narrowpeak = MainConfiguration.query_tables(str(query_species) + "_narrowpeak")
        query_table_gene = MainConfiguration.query_tables(str(query_species) + "_genecoord")

        if query_gene_id:
            # Process query, split into unique elements
            geneid_list_final = [indv_gene for indv_gene in query_gene_id.split(";")]
            geneid_list_final = [indv_gene for indv_gene in geneid_list_final if indv_gene != '']
            # Deploy task
            async_result = celery_tasks.atac_peak_location_search_handler.delay(
                query_species, query_flankingregion_length, query_table_narrowpeak, query_table_gene, geneid_list_final,
                'geneid'
            )
        elif query_range:
            # Deploy task
            async_result = celery_tasks.atac_peak_location_search_handler.delay(
                query_species, query_flankingregion_length, query_table_narrowpeak, query_table_gene, query_range,
                'range'
            )
        else:
            return render(request, 'error.html', {'error_message': 'INVALID_QUERY'})

        # Return async task manager, awaiting frontend recognition
        return render(request, 'epigenome/atac_peakloc_search.html', {'task_id': async_result.task_id})


# ==============MODULE6: MISCELLANEOUS==============
# 6.1 JBROWSE2 SPECIES SELECTOR
def jbrowse_catalog_render(request):
    if request.method == "GET":
        # Accession ID, jbrowse/?id=GCF_000005005.1
        full_name = jbrowse_species_catalog.jbrowse_config_query(str(request.GET.get('id')))
        if full_name:  # Determine whether client is sending a query
            if not full_name == -1:  # Capture key error
                return render(request, 'Jbrowse/jbrowse_insession.html',  # Sending a valid query
                              {'species_full_name': full_name})
        return render(request, 'Jbrowse/jbrowse_main_catalog.html')  # Sending an invalid query or not sending a query


# 6.2 BLAST SERVER
def seqserver(request):
    if request.method == "GET":
        return render(request, 'BLAST/sequenceserver_display.html', {'seq_transfer': ""})  # No autofill
