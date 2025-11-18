import pandas as pd

## Obtain TF DBD of HUMAN DBD FILE
df = pd.read_csv("/home/shey/Documentos/seqlab/TF/HumanTFs.txt", sep="\t")
df_tf = df[df["Is TF?"].str.strip().str.lower() == "yes"]
tf_ids = set(df_tf["Ensembl ID"].dropna().astype(str))
tf_syms = set(df_tf["HGNC symbol"].dropna().astype(str))
print(f"Total TFs únicos extraídos: {len(tf_ids)}")
print(f"Total TFs únicos extraídos: {len(tf_syms)}")    

#Obtein APPRIS information of TF DBD
df_apris_pri = pd.read_csv("/home/shey/Documentos/seqlab/TF/appris_data.appris.txt", sep="\t", low_memory=False)
mask_id  = df_apris_pri["Ensembl Gene ID"].astype(str).isin(tf_ids)
mask_sym = df_apris_pri["Gene name (HGNC)"].astype(str).isin(tf_syms)
db_TF_appris = df_apris_pri[mask_id | mask_sym].copy()
db_TF_appris.to_csv("/home/shey/Documentos/seqlab/TF/appris_TF.csv")






