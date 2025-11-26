import pandas as pd

# Cargar tabla curada (ya contiene la columna 'TF name')
dbd_df = pd.read_csv(
    "/home/shey/Documentos/seqlab/TF/Human_DBDs_CisBP_2.0.txt", sep="\t"
)  # 21277 rows
dbd_df = dbd_df.drop_duplicates(
    subset=["Protein_ID", "Aas", "From", "To"]
)  # 21257 (ENSP00000480026 tenia 11 filas duplicadas y 9 isoformas)

tf_df = pd.read_csv("/home/shey/Documentos/seqlab/TF/TF_completeIDs.csv")
len(tf_df["Ensembl Gene ID"].unique())
len(tf_df["Translation ID"].unique())

## create new db with dbd
#
# df_long = (
#     tf_df
#     .merge(dbd_df[["Protein_ID","Pfam","From","To"]],
#            how="left",
#            left_on="Translation ID",
#            right_on="Protein_ID")
#     .drop(columns="Protein_ID")                      # <- elimina la redundante
#     .dropna(subset=["From","To"])
#     .sort_values(["Translation ID","From","To"])
#     .reset_index(drop=True)
# )
# df_long["dbd_idx"] = df_long.groupby("Translation ID").cumcount() + 1
#
# # incluye el aminoácido en la posición To
# df_long["DBD_seq"] = df_long.apply(
#     lambda r: r["Sequence aa"][int(r["From"]):int(r["To"]) + 1],
#     axis=1
# )
# # sort
# df_long["_iso_n"] = pd.to_numeric(
#     df_long["Isoform name"].str.extract(r"(\d+)$", expand=False),
#     errors="coerce"
# )
#
# df_long = df_long.sort_values(
#     ["Gene name (HGNC)", "_iso_n", "Isoform name", "Transcript ID", "From", "To"],
#     na_position="last"
# ).drop(columns="_iso_n").reset_index(drop=True)

# df_long.to_csv("/home/shey/Documentos/inflamation/TF/dbd_tf.csv", index=False)
# len(df_long["Ensembl Gene ID"].unique())


# 1) (opcional) quedarnos con lo necesario y asegurar ints
dbd_min = (
    dbd_df[["Protein_ID", "Pfam", "From", "To"]]
    .dropna(subset=["From", "To"])
    .astype({"From": int, "To": int})
)

# 2) merge: duplica filas cuando hay varios DBD por Translation ID
df_long = tf_df.merge(
    dbd_min, how="left", left_on="Translation ID", right_on="Protein_ID"
).drop(columns="Protein_ID")

# 3) secuencia del dominio (0-based INCLUSIVO: [From, To])
# Asegura enteros (permite NA con Int64)
df_long["From"] = pd.to_numeric(df_long["From"], errors="coerce").astype("Int64")
df_long["To"] = pd.to_numeric(df_long["To"], errors="coerce").astype("Int64")

# 0-based INCLUSIVO: [From, To] -> seq[f:t+1]
df_long["DBD_seq"] = df_long.apply(
    lambda r: (
        r["Sequence aa"][int(r["From"]) : int(r["To"]) + 1]
        if (
            pd.notna(r["From"])
            and pd.notna(r["To"])
            and isinstance(r["Sequence aa"], str)
        )
        else pd.NA
    ),
    axis=1,
)

# 4) (opcional) ordenar bonito y numerar dominios
# df_long = df_long.sort_values(["Gene name (HGNC)", "Isoform name", "Translation ID", "From", "To"]).reset_index(drop=True)
df_long["dbd_idx"] = df_long.groupby("Translation ID")["From"].transform(
    lambda s: s.notna().cumsum().where(s.notna())
)
df_long.to_csv("/home/shey/Documentos/inflamation/TF/TF_dbd_clean.csv", index=False)


################################

##SI TIENES UNA LISTA DE ISOFORMAS QUE NO DEBEN USARSE (POR EJEMPLO QUE LOS DOMINIOS ERAN MAS GRANDE QUE LA SECUENCIA TOTAL DE LA ISOFORMA)
isoformas_fuera_por_To = []
import pandas as pd
import numpy as np

# --- 0) Set de isoformas "malas" (las 27) ---
bad_isoforms = set(isoformas_fuera_por_To["Isoform name"])

# Mapear esas isoformas a sus Translation ID (porque el merge es por Translation ID)
bad_tids = set(
    tf_df.loc[tf_df["Isoform name"].isin(bad_isoforms), "Translation ID"]
    .dropna()
    .astype(str)
    .unique()
)

# --- 1) Base de DBD (solo columnas necesarias, ints) ---
dbd_min = (
    dbd_df[["Protein_ID", "Pfam", "From", "To"]].dropna(subset=["From", "To"]).copy()
)
dbd_min["From"] = pd.to_numeric(dbd_min["From"], errors="coerce").astype("Int64")
dbd_min["To"] = pd.to_numeric(dbd_min["To"], errors="coerce").astype("Int64")

# EXCLUIR los Translation ID que están en las isoformas que se pasan
dbd_min = dbd_min[~dbd_min["Protein_ID"].astype(str).isin(bad_tids)]

# --- 2) Merge: duplica filas cuando hay varios DBD por Translation ID ---
df_long = tf_df.merge(
    dbd_min, how="left", left_on="Translation ID", right_on="Protein_ID"
).drop(columns="Protein_ID")

# --- 3) Secuencia del dominio (0-based INCLUSIVO: [From, To]) ---
# (Si tus coords fuente fueran 1-based, ajusta a: f = r["From"]-1; t = r["To"]-1)
df_long["From"] = pd.to_numeric(df_long["From"], errors="coerce").astype("Int64")
df_long["To"] = pd.to_numeric(df_long["To"], errors="coerce").astype("Int64")


def slice_or_na(r):
    f, t, seq = r["From"], r["To"], r["Sequence aa"]
    if pd.isna(f) or pd.isna(t) or not isinstance(seq, str):
        return pd.NA
    try:
        return seq[int(f) : int(t) + 1]  # 0-based inclusivo
    except Exception:
        return pd.NA


df_long["DBD_seq"] = df_long.apply(slice_or_na, axis=1)

# --- 4) Numerar dominios por Translation ID (solo los que sí quedaron) ---
df_long["dbd_idx"] = df_long.groupby("Translation ID")["From"].transform(
    lambda s: s.notna().cumsum().where(s.notna())
)
