#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main
@File    ：views.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2023/7/22 14:00
@IDE     ：PyCharm
-----------------------------
Description: The main views function for the sesame information resource, a submodule for HEMU database.
Input:
Output:

Environment requirement (if any): HEMUdb_new (Remote Server, Ubuntu 20.04)
- Dependencies: Django(v3.1), Project files, etc.
-----------------------------
Notes:

-----------------------------
File Revisions:
    2023/7/22: Version 1 - Creation
"""
from django.shortcuts import render, HttpResponse, redirect

# App configuration
from Sesamedb.Main_scripts import MainConfiguration


# ==============HOME PAGE AND STATIC PAGES==============
# HOME PAGE RENDERER
def init_scr(request):
    return render(request, 'sdhome.html')


def init_scr_redirect(request):
    return redirect('/hezhi/home')


# ==============SECTION 1: GENOME==============
# Genome Jbrowse2 instance renderer
def genome_jbr_handler(request):
    if request.method == "GET":

        query_species = request.GET.get('id')
        if query_species:
            if not MainConfiguration.jbrowse_species_validation(query_species) == 102:
                return render(request, 'sdgenome/sdgenome_display.html',
                              {'sp_identifier': MainConfiguration.jbrowse_species_validation(query_species)})
            else:
                return
        return render(request, 'sdgenome/sdgenome_selector.html')


