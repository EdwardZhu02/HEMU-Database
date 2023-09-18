#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project ：HEMU_Database_Main 
@File    ：jbrowse_species_catalog.py
@Author  ：Edward Zhu, zhuyzh37@mail2.sysu.edu.cn
@Date    ：2023/1/10 12:27 
@IDE     ：PyCharm 
-----------------------------
Description: Return redirect links for jbrowse
Input:
Output:

Environment requirement (if any): HEMUdb_new (Remote Server, Ubuntu 20.04)
- Dependencies: None
-----------------------------
Notes:

-----------------------------
File Revisions:
    2023/1/10: Version 1 - Creation
"""
jbrowse_config_pos_dict = {
    # Zea
    "GCF_000005005.2": "GCF_000005005.2-Zea_mays_B73_RefGen_v4",
    "GCF_902167145.1": "GCF_902167145.1-Zea_mays_REFERENCE_NAM_5.0",

    # Sorghum
    "GCF_000003195.3": "GCF_000003195.3-Sorghum_bicolor_BTx623",
    "GCA_903166325.1": "GCA_903166325.1-Sorghum_bicolor_TX436",
    "GCA_903166285.1": "GCA_903166285.1-Sorghum_bicolor_TX2783",
    "GCA_003482435.1": "GCA_003482435.1-Sorghum_bicolor_Tx430",
    "GCA_015952705.1": "GCA_015952705.1-Sorghum_bicolor_Rio",

    # Saccharum
    "GCA_020631735.1": "GCA_020631735.1-Saccharum_officinarum",
    "GCA_022457205.1": "GCA_022457205.1-Saccharum_spontaneum",

    # Coix
    "GCA_009763385.1": "GCA_009763385.1-Coix_lacryma_jobi_var_lacrymajobi",
    "GCA_009758035.1": "GCA_009758035.1-Coix_lacryma_jobi_var_mayuen",
    "GCA_009725075.1": "GCA_009725075.1-Coix_aquatica",

    # Miscanthus
    "GCA_904845875.1": "GCA_904845875.1-Miscanthus_lutarioriparius",

    # Bothriochloa
    "GCA_023333625.1": "GCA_023333625.1-Bothriochloa_decipiens",
    "GCA_022036555.1": "GCA_022036555.1-Microstegium_vimineum",
    "GCA_018135685.1": "GCA_018135685.1-Themeda_triandra",
}

jbrowse_syntenysession_name_dict = {
    # Zea - Sorghum
    "B73v4_BTx623": "B73v4_BTx623",

}


def jbrowse_syntenysession_query(species_query):
    try:
        return str(jbrowse_syntenysession_name_dict[str(species_query)])
    except KeyError:
        return -1  # Key error


def jbrowse_config_query(species_query):
    try:
        return str(jbrowse_config_pos_dict[str(species_query)])
    except KeyError:
        return -1  # Key error

