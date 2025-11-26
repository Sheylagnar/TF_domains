import pandas as pd

# ================= Obtain DBD information of TFs =====================
# 1) Load DBD database and TF list
dbd_df = pd.read_csv(
    "/home/shey/Documentos/seqlab/TF/Human_DBDs_CisBP_2.0.txt", sep="\t"
)  # 21277 rows
dbd_df = dbd_df.drop_duplicates(
    subset=["Protein_ID", "Aas", "From", "To"]
)  # 21257 (ENSP00000480026 tenia 11 filas duplicadas y 9 isoformas)
dbd_min = (
    dbd_df[["Protein_ID", "Pfam", "From", "To"]]
    .dropna(subset=["From", "To"])
    .astype({"From": int, "To": int})
)
tf_df = pd.read_csv("/home/shey/Documentos/seqlab/TF/TF_completeIDs.csv")

# 2) merge: duplica filas cuando hay varios DBD por Translation ID
df_long = tf_df.merge(
    dbd_min, how="left", left_on="Translation ID", right_on="Protein_ID"
).drop(columns="Protein_ID")
df_long["From"] = pd.to_numeric(df_long["From"], errors="coerce").astype("Int64")
df_long["To"] = pd.to_numeric(df_long["To"], errors="coerce").astype("Int64")

# 2b) Verificar que las coordenadas no exceden la longitud de la secuencia
df_long["seq_len"] = df_long["Sequence aa"].fillna("").str.len()

mask_fuera_rango = (
    df_long["Sequence aa"].notna()
    & df_long["To"].notna()
    & (df_long["To"] >= df_long["seq_len"])  # To no puede ser >= len(seq)
)

print("Filas con To >= len(Sequence aa):", mask_fuera_rango.sum())
print(
    df_long.loc[
        mask_fuera_rango, ["Ensembl Gene ID", "Translation ID", "From", "To", "seq_len"]
    ].head()
)


# 3) Extrar la secuencia en base a la coordenada 0-based INCLUSIVO
def extraer_dbd_seq(row):
    if (
        pd.notna(row["From"])
        and pd.notna(row["To"])
        and isinstance(row.get("Sequence aa"), str)
    ):
        f = int(row["From"])
        t = int(row["To"])
        seq = row["Sequence aa"]
        if f > t or t >= len(seq):
            return pd.NA
        return seq[f : t + 1]
    return pd.NA


df_long["DBD_seq"] = df_long.apply(extraer_dbd_seq, axis=1)

# 4) Create index for DBD per Translation ID
df_long["dbd_idx"] = (
    df_long.groupby("Translation ID")["From"]
    .transform(lambda s: s.notna().cumsum().where(s.notna()))
    .astype("Int64")
)
df_long.to_csv("/home/shey/Documentos/seqlab/TF/TF_dbd_data .csv", index=False)
