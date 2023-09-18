import os

import pandas as pd
from sqlalchemy import create_engine

# Used for sqlalchemy engine services
sql_engine = create_engine("")

# Final database sheet
final_TE_datasheet_df = pd.DataFrame()
# Counter for samples
samples_completed = 1

with open("/mnt/f/Sbico-te-exp-sheet-assembly/step3_TEexpDB_construct.log", mode='w') as out_fh:
    # Constructing dict, sample_id -> sample_tissue
    sampleinfo_data = pd.read_csv("/mnt/f/Sbico-te-exp-sheet-assembly/Sbico_sampleinfo.csv", low_memory=False, sep=",")
    sampleinfo_df = pd.DataFrame(sampleinfo_data)
    sample_tissue_dict = {}  # sample_id -> sample_tissue
    for sample_entry in sampleinfo_df.index:
        sample_tissue_dict[sampleinfo_df.loc[sample_entry].sample_id] = str(
            sampleinfo_df.loc[sample_entry].sample_tissue).lower()
    out_fh.write("Sample-tissue dict construction completed. Length: %d" % len(sample_tissue_dict.keys()))
    print("Sample-tissue dict construction completed. Length: %d" % len(sample_tissue_dict.keys()))

    for indv_file in os.listdir("/mnt/f/Sbico-te-exp-sheet-assembly/TEcount_final"):  # ex. SRR0001_TE.csv

        if indv_file.endswith("_TE.csv"):
            TE_data = pd.read_csv("/mnt/f/Sbico-te-exp-sheet-assembly/TEcount_final/" + indv_file, low_memory=False, sep=",")
            TE_df = pd.DataFrame(TE_data)

            # Delete single-copy TEs
            TE_df.drop(TE_df[TE_df['te_id'].str.contains(pat='Seq', regex=False)].index, inplace=True)

            TE_df_final = TE_df
            # Insert sample_id at the second column
            TE_df_final.insert(loc=1, column='sample_id', value=indv_file.rstrip("_TE.csv"))

            # Insert tissue information at the last column
            try:
                tissue_type = str(sample_tissue_dict[indv_file.rstrip("_TE.csv")])
            except KeyError:
                out_fh.write("Index not found for %s, considering tissue unknown.\n" % indv_file.rstrip("_TE.csv"))
                tissue_type = "unknown"
            TE_df_final.insert(loc=6, column='tissue_type', value=tissue_type)

            final_TE_datasheet_df = pd.concat([final_TE_datasheet_df, TE_df_final])
            out_fh.write("Completed %d samples\n" % samples_completed)
            print("Completed %d samples\n" % samples_completed)
            samples_completed += 1
            # break  # Diag purposes
    out_fh.write("\n\nCreating csv backup..")
    print("\n\nCreating csv backup..")
    final_TE_datasheet_df.to_csv("/mnt/f/Sbico-te-exp-sheet-assembly/backup_Sbico_sheet.csv")
    out_fh.write("\n\nCreating database sheet..")
    print("\n\nCreating database sheet..")
    final_TE_datasheet_df.to_sql(name='sbico_te', con=sql_engine, if_exists="replace", index=False, chunksize=100)
    out_fh.write("All done")
    print("All done")
