import pandas as pd

ptm_cols = ["Protein", "UniProt", "Position", "PTM", "PMID", "Peptide sequence"]
df_ptm = pd.read_csv(
    "/home/shey/Documentos/seqlab/TF/dbptm/Acetylation",
    sep="\t",
    names=ptm_cols,
)
tf2 = pd.read_csv("/home/shey/Documentos/seqlab/TF/TF_dbd_data.csv", sep=",")
tf = pd.read_csv("/home/shey/Documentos/seqlab/TF/appris_TF_complete.csv", sep=",")

df_ptm = df_ptm[df_ptm["Protein"].str.endswith("_HUMAN", na=False)].copy()
df_ptm["UniProt"] = (
    df_ptm["UniProt"].astype(str).str.strip().str.split("-").str[0].str.upper()
)

tf["UniProt_ID"] = (
    tf["UniProt_ID_norm"].astype(str).str.strip().str.split("-").str[0].str.upper()
)

mask = df["UniProt"].isin(tf["UniProt_ID"])
df_match = df[mask].copy()

##SECOND WAY

# 3) (A) Mapa UniProt -> Ensembl Gene ID (pueden ser varios por UniProt)
map_u2g = (
    tf[["UniProt_ID_norm", "Ensembl Gene ID"]]
    .drop_duplicates()
    .rename(columns={"UniProt_ID_norm": "UniProt"})
)

# 3) (B) Anexar Ensembl Gene ID(s) a cada fila de df_match según su UniProt
m1 = df_match.merge(map_u2g, on="UniProt", how="left")

# 4) Traer TODAS las isoformas de esos Ensembl Gene ID (sin restringir por UniProt)
cols_tf = [
    "Ensembl Gene ID",
    "HGNC symbol",
    "Isoform name",
    "Uniprot_isoform",
    "Transcript ID",
    "Translation ID",
    "Sequence aa",
    "UniProt_ID_norm",
]
m = m1.merge(tf[cols_tf], on="Ensembl Gene ID", how="left")

# 5) Buscar el péptido (21 aa) en la secuencia de la isoforma
m["Peptide_no_dash"] = (
    m["Peptide sequence"].astype(str).str.replace("-", "", regex=False)
)
m = m[
    m.apply(
        lambda r: isinstance(r["Sequence aa"], str)
        and r["Peptide_no_dash"] in r["Sequence aa"],
        axis=1,
    )
].copy()

# 6) Posiciones y AA modificado
m["pep_start_in_isoform"] = m.apply(
    lambda r: r["Sequence aa"].find(r["Peptide_no_dash"]) + 1, axis=1
)
m["Modified aa"] = m["Peptide_no_dash"].str[10]
m["Real Position"] = m["pep_start_in_isoform"] + 10
m = m.rename(columns={"Position": "Reference position"})

# 7) Resultado base
result = (
    m[
        [
            "PTM",
            "Ensembl Gene ID",
            "Transcript ID",
            "Translation ID",
            "HGNC symbol",
            "Protein",
            "UniProt",
            "Reference position",
            "Real Position",
            "Isoform name",
            "Uniprot_isoform",
            "Modified aa",
            "Peptide sequence",
            "PMID",
        ]
    ]
    .sort_values(["Protein", "Reference position", "Uniprot_isoform"])
    .reset_index(drop=True)
)


# 8) Index por sitio dentro de cada proteína (1..n)
def add_ac_idx(g):
    keys = (
        g[["Reference position", "Peptide sequence"]]
        .drop_duplicates()
        .sort_values(["Reference position", "Peptide sequence"])
        .reset_index(drop=True)
    )
    mapping = {
        (r["Reference position"], r["Peptide sequence"]): i + 1
        for i, r in keys.iterrows()
    }
    g["Index"] = g.apply(
        lambda r: mapping[(r["Reference position"], r["Peptide sequence"])], axis=1
    )
    return g


result = (
    result.groupby("Protein", group_keys=False)
    .apply(add_ac_idx)
    .sort_values(["HGNC symbol", "Index", "Uniprot_isoform"])
    .reset_index(drop=True)
)

# 9) Poner "Index" al lado de "HGNC symbol"
col = result.pop("Index")
result.insert(result.columns.get_loc("HGNC symbol") + 1, "Index", col)

result.to_csv("/home/shey/Documentos/seqlab/TF/ACET_result.csv", index=False)


import os
import pandas as pd

# --- Rutas ---
PTM_DIR = "/home/shey/Documentos/seqlab/TF/dbptm"
TF_CSV = "/home/shey/Documentos/seqlab/TF/TF_completeIDs.csv"
OUT_CSV = "/home/shey/Documentos/seqlab/TF/PTM_aggregated.csv"

# --- Carga TF + maps básicos ---
tf = pd.read_csv(TF_CSV)
map_u2g = (
    tf[["UniProt_ID", "Ensembl Gene ID"]]
    .drop_duplicates()
    .rename(columns={"UniProt_ID": "UniProt"})
)
cols_tf = [
    "Ensembl Gene ID",
    "Gene name (HGNC)",
    "Isoform name",
    "Transcript ID",
    "Translation ID",
    "Sequence aa",
    "UniProt_ID",
]


# --- Helper de un archivo DBPTM (tsv sin header) ---
def process_ptm_file(path: str) -> pd.DataFrame:
    ptm_cols = ["Protein", "UniProt", "Position", "PTM", "PMID", "Peptide sequence"]
    df = pd.read_csv(
        path, sep="\t", names=ptm_cols, dtype=str, engine="python", on_bad_lines="skip"
    )
    if df.empty:
        return pd.DataFrame()

    # humanos + normalización
    df = df[df["Protein"].str.endswith("_HUMAN", na=False)].copy()
    if df.empty:
        return pd.DataFrame()
    df["UniProt"] = (
        df["UniProt"].astype(str).str.strip().str.split("-").str[0].str.upper()
    )
    df["PTM"] = df["PTM"].fillna("").replace("", os.path.basename(path))

    # match por UniProt presente en TF
    df = df[df["UniProt"].isin(tf["UniProt_ID"])].copy()
    if df.empty:
        return pd.DataFrame()

    # Ensembl y todas las isoformas
    m = df.merge(map_u2g, on="UniProt", how="left").merge(
        tf[cols_tf], on="Ensembl Gene ID", how="left"
    )
    if m.empty:
        return pd.DataFrame()

    # buscar el 21-mer en la isoforma
    m["pep21"] = m["Peptide sequence"].astype(str).str.replace("-", "", regex=False)
    m = m[
        m.apply(
            lambda r: isinstance(r["Sequence aa"], str)
            and r["pep21"] in r["Sequence aa"],
            axis=1,
        )
    ]
    if m.empty:
        return pd.DataFrame()

    # posición real y AA modificado (centro del 21-mer)
    m["pep_start"] = m.apply(lambda r: r["Sequence aa"].find(r["pep21"]) + 1, axis=1)
    m["Modified aa"] = m["pep21"].str[10]
    m["Real Position"] = m["pep_start"] + 10
    m = m.rename(columns={"Position": "Reference position"})

    # columnas finales + orden básico
    # columnas finales + orden por HGNC + Reference position numérica
    res = m[
        [
            "PTM",
            "Ensembl Gene ID",
            "Transcript ID",
            "Translation ID",
            "Protein",
            "UniProt",
            "Gene name (HGNC)",
            "Isoform name",
            "Reference position",
            "Real Position",
            "Modified aa",
            "Peptide sequence",
            "PMID",
        ]
    ].copy()

    # ordenar por posición como número (los no numéricos quedan al final)
    res["_refnum"] = pd.to_numeric(res["Reference position"], errors="coerce").fillna(
        10**9
    )
    res = res.sort_values(["Gene name (HGNC)", "_refnum"]).reset_index(drop=True)

    if res.empty:
        return res

    # Index por sitio dentro de Protein (Reference position + Peptide)
    # Index por sitio dentro de Protein (Reference position + Peptide)
    res["site_key"] = (
        res["Reference position"].astype(str) + "|" + res["Peptide sequence"]
    )
    res = res.sort_values(["Gene name (HGNC)", "_refnum", "Peptide sequence"])
    res["Index"] = res.groupby("Gene name (HGNC)")["site_key"].transform(
        lambda s: (s != s.shift()).cumsum()
    )
    res = res.drop(columns=["site_key", "_refnum"])

    # mover Index al lado de Gene name (HGNC)
    col = res.pop("Index")
    res.insert(res.columns.get_loc("Isoform name") + 1, "Index", col)
    return res


# --- Limpia salida y procesa todo en orden alfabético ---
if os.path.exists(OUT_CSV):
    os.remove(OUT_CSV)

files = sorted(
    os.path.join(PTM_DIR, f)
    for f in os.listdir(PTM_DIR)
    if os.path.isfile(os.path.join(PTM_DIR, f))
)

total = 0
for i, fp in enumerate(files, 1):
    name = os.path.basename(fp)
    print(f"[{i}/{len(files)}] {name} …", flush=True)
    try:
        df_res = process_ptm_file(fp)
    except Exception as e:
        print(f"   [warn] {name}: {e}", flush=True)
        continue
    if df_res.empty:
        print(f"   {name}: 0 filas.", flush=True)
        continue

    header = not os.path.exists(OUT_CSV)
    df_res.to_csv(OUT_CSV, mode="a", index=False, header=header)
    total += len(df_res)
    print(f"   {name}: +{len(df_res)} (acum {total})", flush=True)

print(f"[OK] Guardado: {OUT_CSV} | Filas totales: {total}")


#
# cols_tf = [
#     "UniProt_ID_norm", "HGNC symbol", "Uniprot_isoform", "Ensembl Gene ID",
#     "Transcript ID", "Translation ID", "Sequence aa", "Isoform name"
# ]
#
# # 1) Merge 1:N por UniProt base (sin tocar tus IDs)
# m = m1.merge(
#     tf[cols_tf],
#     left_on="UniProt",
#     right_on="UniProt_ID_norm",
#     how="left"
# )
#
# # 2) Quitar solo los guiones del péptido para poder buscarlo en la secuencia
# m["Peptide_no_dash"] = m["Peptide sequence"].astype(str).str.replace("-", "", regex=False)
#
# # 3) Quedarnos con las isoformas cuya 'Sequence aa' contiene ese péptido (21 aa)
# mask = m.apply(
#     lambda r: isinstance(r["Sequence aa"], str) and r["Peptide_no_dash"] in r["Sequence aa"],
#     axis=1
# )
# m = m[mask].copy()
#
# # 4) (opcional) posición 1-based donde inicia el match en la isoforma
# m["pep_start_in_isoform"] = m.apply(
#     lambda r: r["Sequence aa"].find(r["Peptide_no_dash"]) + 1, axis=1
# )
#
# # 5) Aminoácido modificado = 11° de la ventana (asumiendo 21 aa)
# m["Modified aa"] = m["Peptide_no_dash"].str[10]
#
# # 6) Real Position (en la isoforma) = inicio + 10 (porque el aa modificado es el #11)
# m["Real Position"] = m["pep_start_in_isoform"] + 10
#
# # 7) Renombrar 'Position' a 'Reference position' (posición del archivo original)
# m = m.rename(columns={"Position": "Reference position"})
#
# # 8) Salida final (si varias isoformas contienen el péptido, verás filas duplicadas como deseas)
# result = m[
#     [
#         "PTM", "Ensembl Gene ID", "Transcript ID", "Translation ID", "HGNC symbol", "Protein", "UniProt",
#         "Reference position", "Real Position", "Isoform name",
#         "Uniprot_isoform", "Modified aa", "Peptide sequence", "PMID"
#     ]
# ].sort_values(["Protein", "Reference position", "Uniprot_isoform"]).reset_index(drop=True)
#
# def add_ac_idx(g):
#     keys = (g[["Reference position","Peptide sequence"]]
#             .drop_duplicates()
#             .sort_values(["Reference position","Peptide sequence"])
#             .reset_index(drop=True))
#     mapping = { (r["Reference position"], r["Peptide sequence"]): i+1 for i, r in keys.iterrows() }
#     g["Index"] = g.apply(lambda r: mapping[(r["Reference position"], r["Peptide sequence"])], axis=1)
#     return g
#
# result = (
#     result
#     .groupby("Protein", group_keys=False)
#     .apply(add_ac_idx)
#     .sort_values(["HGNC symbol", "Index", "Uniprot_isoform"])
#     .reset_index(drop=True)
# )
# col = result.pop("Index")
# result.insert(result.columns.get_loc("HGNC symbol") + 1, "Index", col)
#
# result.to_csv("/home/shey/Documentos/seqlab/TF/partial_result.csv", index=False)
