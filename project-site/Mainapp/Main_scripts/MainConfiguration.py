"""
MainConfiguration.py
Stores critical information and table names used in the project
"""
tables_dict = {  # Table names
    'coix_exp': 'coix_exp',
    'coix_samp': 'coix_samp',
    'coix_te': 'coix_te',
    'coix_gokegg': 'clacr_gokegg',
    'coix_filename': 'Clacr',
    'coix_genecoord': 'clacr_genecoord',
    'coix_tecoord': 'clacr_tecoord',

    'zea_exp': 'zea_exp',
    'zea_samp': 'zea_samp',
    'zea_te': 'zea_te',
    'zea_gokegg': 'zmaysv4_gokegg',
    'zea_filename': 'Zmaysv4',
    'zea_narrowpeak': 'zmaysv4_narrowpeak',
    'zea_genecoord': 'zmaysv4_genecoord',
    'zea_tecoord': 'zmaysv4_tecoord',

    'sorghum_exp': 'sorghum_exp',
    'sorghum_samp': 'sorghum_samp',
    'sorghum_te': 'sbico_te',
    'sorghum_gokegg': 'sbico_gokegg',
    'sorghum_filename': 'Sbico',
    'sorghum_narrowpeak': 'sbico_narrowpeak',
    'sorghum_genecoord': 'sbico_genecoord',
    'sorghum_tecoord': 'sbico_tecoord',

    'saccharum_exp': 'saccharum_exp',
    'saccharum_samp': 'saccharum_samp',
    'saccharum_gene_raw': 'sspon_gene_raw',
    'saccharum_gokegg': 'sspon_gokegg',
    'saccharum_filename': 'Sspon',
    'saccharum_genecoord': 'sspon_genecoord',
    'saccharum_tecoord': 'sspon_tecoord',

    'miscanthus_exp': 'mluta_exp',
    'miscanthus_samp': 'mluta_samp',
    'miscanthus_te': 'mluta_te',
    'miscanthus_gokegg': 'mluta_gokegg',
    'miscanthus_filename': 'Mluta',
    'miscanthus_genecoord': 'mluta_genecoord',
    'miscanthus_tecoord': 'mluta_tecoord',

    'miscanthus_sine_exp': 'msine_exp',
    'miscanthus_sine_samp': 'msine_samp',
    'miscanthus_sine_te': 'msine_te',
    'miscanthus_sine_gokegg': 'msine_gokegg',
    'miscanthus_sine_filename': 'Msine',

    'chrysopogon_gokegg': 'cserr_gokegg',
    'chrysopogon_filename': 'Cserr',

    'hyparrhenia_gokegg': 'hdipl_gokegg',
    'hyparrhenia_filename': 'Hdipl',

    'themeda_exp': 'ttria_exp',
    'themeda_samp': 'ttria_samp',
    'themeda_gokegg': 'ttria_gokegg',
    'themeda_filename': 'Ttria',
}

sql_user = ""
sql_pwd = ""
sql_host = ""
sql_port = ""
sql_dbname = ""

admin_pwd = ""


def query_tables(tbl2query):
    try:
        return tables_dict[tbl2query]
    except NameError:
        return None


def query_admin_pwd():
    return admin_pwd


def query_sql(identity):
    if identity == "host":
        return sql_host
    elif identity == "user":
        return sql_user
    elif identity == "pwd":
        return sql_pwd
    elif identity == "port":
        return sql_port
    elif identity == "dbname":
        return sql_dbname
    else:
        return None
