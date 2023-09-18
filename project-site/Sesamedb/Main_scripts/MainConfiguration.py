#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU-Extended-Main 
@File    ：MainConfiguration.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2023/08/08 23:27 
@IDE     ：PyCharm 
-----------------------------
Description:
Input:
Output:

One-line usage note: python3 MainConfiguration.py ,
- ARG1:

Environment requirement (if any): base
- Dependencies: None
-----------------------------
Notes:

-----------------------------
File Revisions:
    2023/08/08: Version 1 - Creation
"""


def jbrowse_species_validation(query_species):
    valid_species_list = ["bl123", "wh153"]
    if str(query_species) in valid_species_list:
        return query_species
    else:
        return 102  # species not found
