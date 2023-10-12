import pandas as pd
import numpy as np


def filter_extra_GOlines(enricher_df, go2term_df):
    total_go_term_comp = []
    lines_to_del = []

    enricher_list = np.array(enricher_df).tolist()
    goterm_list = np.array(go2term_df).tolist()

    #enricher_list.remove(enricher_list[0]) # Remove Index Lines
    #goterm_list.remove(goterm_list[0])
    enricher_list = list(enricher_list)
    goterm_list = list(goterm_list)

    for term in goterm_list:
        if not term[0] in total_go_term_comp:
            total_go_term_comp.append(term[0])
    for index,entry in enumerate(enricher_list):
        if entry[0] not in total_go_term_comp:
            lines_to_del.append(index)
            entry.clear()

    enricher_list = [x for x in enricher_list if len(x) > 0]
    enricher_df_new = pd.DataFrame(enricher_list)

    return enricher_df_new