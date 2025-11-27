import os
import pandas as pd

# Paths
PTM_DIR = "/home/shey/Documentos/seqlab/TF/dbptm"
TF_CSV = "/home/shey/Documentos/seqlab/TF/TF_completeIDs.csv"
OUT_CSV = "/home/shey/Documentos/seqlab/TF/PTM_aggregated.csv"

# Load files
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


# --- Helper de un archivo DBPTM  ---
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
    m["pep_start"] = m.apply(lambda r: r["Sequence aa"].find(r["pep21"]), axis=1)
    m["Modified aa"] = m["pep21"].str[10]
    m["Real Position"] = m["pep_start"] + 10
    m["Reference position"] = (
        pd.to_numeric(m["Position"], errors="coerce") - 1
    ).astype("Int64")

    # columnas finales + orden básico
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
