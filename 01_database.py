import pandas as pd

## Obtain TFs of HUMAN DBD FILE
dbd_df = pd.read_csv(
    "/home/shey/Documentos/seqlab/TF/Human_DBDs_CisBP_2.0.txt", sep="\t"
)
dbd_df = dbd_df.drop_duplicates(subset=["Protein_ID", "From", "To"])
tf_ids = set(dbd_df["Gene_ID"].dropna().astype(str))
tf_syms = set(dbd_df["Protein_ID"].dropna().astype(str))
print(f"Total TFs únicos por Gene_ID: {len(tf_ids)}")
print(f"Total TFs únicos por Protein_ID: {len(tf_syms)}")

# Obtain APPRIS information of TFs
df_apris_pri = pd.read_csv(
    "/home/shey/Documentos/seqlab/TF/appris_data.appris.txt", sep="\t", low_memory=False
)
mask_id = df_apris_pri["Ensembl Gene ID"].astype(str).isin(tf_ids)
mask_sym = df_apris_pri["Translation ID"].astype(str).isin(tf_syms)
db_TF_appris = df_apris_pri[mask_id | mask_sym].copy()
print(f"Summary db_TF_appris: {db_TF_appris.nunique()}")
db_TF_appris.to_csv("/home/shey/Documentos/inflamation/TF/appris_TF.csv")
