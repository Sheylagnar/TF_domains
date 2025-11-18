import time
import pandas as pd
import requests, sys

# Paths
INPUT_CSV = r"/home/shey/Documentos/inflamation/TF/appris_TF.csv"
OUTPUT_CSV = r"/home/shey/Documentos/inflamation/TF/TF_completeIDs.csv"

# Ensembl headers
BASE_URL = "https://rest.ensembl.org"

session = requests.Session()
session.headers.update(
    {
        "Accept": "application/json",
        "Content-Type": "application/json",
        "User-Agent": "IsoformUniProtMapper/1.0 (mailto:sheyla.carmen.2@gmail.com)",
    }
)


# Helper functions
def clean_ids(ids):
    """Removes NaN, spaces, and empty spaces; returns a single list"""
    s = pd.Series(ids).dropna().astype(str).str.strip()
    s = s[s != ""]
    return s.unique().tolist()


def get_sequences(translation_ids, chunk=50, pause=0.2):
    """Returns dict TranslationID -> sequence aa"""
    ids = clean_ids(translation_ids)
    out = {}
    for i in range(0, len(ids), chunk):
        batch = ids[i : i + chunk]
        r = session.post(
            f"{BASE_URL}/sequence/id",
            json={"ids": batch, "type": "protein"},
            timeout=60,
        )
        if not r.ok:
            raise RuntimeError(f"HTTP {r.status_code} from Ensembl: {r.text[:300]}")
        try:
            data = r.json()
        except Exception:
            ct = r.headers.get("Content-Type")
            txt = r.text[:500]
            raise RuntimeError(f"Respuesta no-JSON (Content-Type={ct}):\n{txt}")
        for item in data:
            out[item.get("id")] = item.get("seq")
        time.sleep(pause)
    return out


def fetch_transcript_names(ids, chunk=100, pause=0.2):
    ids = clean_ids(ids)
    out = {}
    for i in range(0, len(ids), chunk):
        batch = ids[i : i + chunk]
        r = session.post(
            f"{BASE_URL}/lookup/id",
            params={"expand": 0},
            json={"ids": batch},
            timeout=60,
        )
        r.raise_for_status()
        data = r.json()
        if isinstance(data, dict):
            it = data.items()
        else:
            it = []
            for d in data:
                if not isinstance(d, dict):
                    continue
                _id = d.get("id")
                if _id is None:
                    continue
                it.append((_id, d))
        for _id, obj in it:
            if not isinstance(obj, dict):
                continue
            name = obj.get("display_name") or obj.get("external_name")
            out[_id] = name if name is not None else _id
        time.sleep(pause)
    return out


def uni_from_ensps(ensps):
    ids = clean_ids(ensps)
    out = {}
    for ensp in ensps:
        r = session.get(f"{BASE_URL}/xrefs/id/{ensp}", timeout=30)
        if not r.ok:
            out[ensp] = None
            continue
        xs = r.json()
        sp = next(
            (
                x
                for x in xs
                if x.get("dbname") in ("UniProtKB/Swiss-Prot", "Uniprot/SWISSPROT")
            ),
            None,
        )
        tr = next(
            (
                x
                for x in xs
                if x.get("dbname") in ("UniProtKB/TrEMBL", "Uniprot/SPTREMBL")
            ),
            None,
        )
        pick = sp or tr
        out[ensp] = (
            (pick.get("primary_id") or pick.get("display_id") or pick.get("id"))
            if pick
            else None
        )
    return out


## Pipeline
df = pd.read_csv(INPUT_CSV, sep=",")

# Sequences aa
seq_map = get_sequences(df["Translation ID"])
df["Sequence aa"] = df["Translation ID"].astype(str).str.strip().map(seq_map)
print(df[["Translation ID", "Sequence aa"]].head())

# Isoform name
name_map = fetch_transcript_names(df["Transcript ID"])
df["Isoform name"] = df["Transcript ID"].astype(str).str.strip().map(name_map)

# UniProt ID
mask = df["Translation ID"].notna() & (
    df["Translation ID"].astype(str).str.strip() != ""
)
ensps = df.loc[mask, "Translation ID"]
ensp2uni = uni_from_ensps(ensps)
df["UniProt_ID"] = df["Translation ID"].astype(str).str.strip().map(ensp2uni)

df.to_csv("/home/shey/Documentos/inflamation/TF/TF_completeIDs.csv")
