import shutil
from pathlib import Path


def clean_tmp_files():

    dirs_to_del = [
        'Mainapp/static/Temp_R_DEprojects',
        'Mainapp/static/Temp_R_epigenome',
        'Mainapp/static/Temp_R_genefam',
        'Mainapp/static/Temp_R_gokegg',
        'Mainapp/static/Temp_R_html',
        'Mainapp/static/Temp_R_TE',
        'Mainapp/static/Temp_R_wgcna',
    ]

    def del_all(path):
        file_count = 0
        for elm in Path(path).glob('*'):
            file_count += 1
            elm.unlink() if elm.is_file() else shutil.rmtree(elm)  # delete folder and file
        return file_count

    file_count_total = 0
    for indv_dir in dirs_to_del:
        file_count_total += del_all(indv_dir)
    return file_count_total
