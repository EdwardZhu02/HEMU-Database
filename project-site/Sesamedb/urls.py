#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU-Extended-Main 
@File    ：urls.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2023/07/22 14:28 
@IDE     ：PyCharm 
-----------------------------
Description:
Input:
Output:

One-line usage note: python3 urls.py ,
- ARG1:

Environment requirement (if any): base
- Dependencies: None
-----------------------------
Notes:

-----------------------------
File Revisions:
    2023/07/22: Version 1 - Creation
"""
from django.urls import path, include, re_path
from Sesamedb import views

urlpatterns = [

    # Home Page
    path('', views.init_scr_redirect),
    path('home', views.init_scr, name='sd_home_page'),

    # Section1: Genome
    path('genome', views.genome_jbr_handler, name='sd_genome_selector'),

    ]