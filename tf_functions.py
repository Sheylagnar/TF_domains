import time
import pandas as pd
import requests, sys
# Cargar tabla curada (ya contiene la columna 'TF name')
df = pd.read_csv("/home/shey/Documentos/inflamation/TF/appris_TF.csv", sep=",")

#BASE
server = "https://rest.ensembl.org"
ext = "/sequence/id"
headers = {"Content-Type": "application/json", "Accept": "application/json"}
r = requests.post(server + ext, headers=headers, data='{ "ids" : ["ENST00000696671", "ENST00000696670" ] }')

if not r.ok:
    r.raise_for_status()
    sys.exit()

decoded = r.json()
print(repr(decoded))


ENSEMBL_URL = "https://rest.ensembl.org/sequence/id"
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json"
}

def get_sequences(transcript_ids, chunk=50, pause=0.2):
    """Devuelve dict TranscriptID -> secuencia (DNA por defecto del endpoint)."""
    ids = pd.Series(transcript_ids).dropna().astype(str).str.strip()
    ids = ids[ids != ""].unique().tolist()
    out = {}
    for i in range(0, len(ids), chunk):
        batch = ids[i:i+chunk]
        r = requests.post(ENSEMBL_URL, headers=HEADERS, json={"ids": batch, "type":"protein"}, timeout=60)
        if not r.ok:
            # Verás el código y un pedazo de respuesta (útil si es HTML de error)
            raise RuntimeError(f"HTTP {r.status_code} from Ensembl: {r.text[:300]}")
        try:
            data = r.json()
        except Exception:
            # Cuando no es JSON (p. ej., HTML de rate limit/maintenance)
            ct = r.headers.get("Content-Type")
            txt = r.text[:500]
            raise RuntimeError(f"Respuesta no-JSON (Content-Type={ct}):\n{txt}")
        for item in data:
            out[item.get("id")] = item.get("seq")
        time.sleep(pause)  # sé amable con la API
    return out

df["Translation ID"] = df["Translation ID"].astype(str).str.strip()
seq_map = get_sequences(df["Translation ID"].tolist())
df["Sequence aa"] = df["Translation ID"].map(seq_map)
print(df[["Translation ID", "Sequence aa"]].head())



##ISOFORM NAME
LOOKUP_URL = "https://rest.ensembl.org/lookup/id"
HEADERS = {
    "Accept": "application/json",
    "Content-Type": "application/json",
    "User-Agent": "IsoformFetcher/1.0 (mailto:tucorreo@ejemplo.com)"
}

def fetch_transcript_names(ids, chunk=100, pause=0.2):
    ids = pd.Series(ids).dropna().astype(str).str.strip().tolist()
    out = {}
    for i in range(0, len(ids), chunk):
        batch = ids[i:i+chunk]
        r = requests.post(LOOKUP_URL, headers=HEADERS,
                          params={"expand": 0},
                          json={"ids": batch}, timeout=60)
        r.raise_for_status()
        data = r.json()
        # Soporta respuesta tipo lista o dict (según versión)
        if isinstance(data, dict):
            it = data.items()
        else:
            it = [(d.get("id"), d) for d in data]
        for _id, obj in it:
            out[_id] = obj.get("display_name") or obj.get("external_name")
        time.sleep(pause)
    return out

# Añadir la columna al DataFrame
df["Transcript ID"] = df["Transcript ID"].astype(str).str.strip()
name_map = fetch_transcript_names(df["Transcript ID"])
df["Isoform name"] = df["Transcript ID"].map(name_map)


#OBTEIN UNIPROT
import requests, pandas as pd

s = requests.Session()
s.headers.update({
    "Accept": "application/json",
    "User-Agent": "UniProtMapper/1.0 (mailto:tucorreo@ejemplo.com)"
})

def uni_from_ensps(ensps):
    out = {}
    for ensp in ensps:
        r = s.get(f"https://rest.ensembl.org/xrefs/id/{ensp}", timeout=30)
        if not r.ok:
            out[ensp] = None
            continue
        xs = r.json()
        sp = next((x for x in xs if x.get("dbname") in ("UniProtKB/Swiss-Prot","Uniprot/SWISSPROT")), None)
        tr = next((x for x in xs if x.get("dbname") in ("UniProtKB/TrEMBL","Uniprot/SPTREMBL")), None)
        pick = sp or tr
        out[ensp] = (pick.get("primary_id") or pick.get("display_id") or pick.get("id")) if pick else None
    return out

# aplicar sobre tu df
mask = df["Translation ID"].notna() & (df["Translation ID"].astype(str).str.strip() != "")
ensps = df.loc[mask, "Translation ID"].astype(str).str.strip().unique().tolist()

ensp2uni = uni_from_ensps(ensps)
df["UniProt_ID"] = df["Translation ID"].astype(str).str.strip().map(ensp2uni)

df.to_csv("/home/shey/Documentos/inflamation/TF/TF_completeIDs.csv")