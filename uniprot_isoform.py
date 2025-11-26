import pandas as pd
import time, re, requests
from functools import lru_cache


# 1) Funciones para uniprot

S = requests.Session()
ENSEMBL = "https://rest.ensembl.org"
UNI = "https://rest.uniprot.org"
HJSON = {"Content-Type":"application/json"}

def _get(url, **kwargs):
    timeout = kwargs.pop("timeout", 30)
    r = S.get(url, timeout=timeout, **kwargs)
    r.raise_for_status()
    return r

@lru_cache(maxsize=10000)
def get_ensp_from_enst(enst: str):
    if not enst or not str(enst).startswith("ENST"):
        return None
    r = _get(f"{ENSEMBL}/lookup/id/{enst}", params={"expand":1}, headers=HJSON)
    trans = (r.json().get("Translation") or {})
    return trans.get("id")

@lru_cache(maxsize=10000)
def ensembl_xrefs_isoform(ensp: str):
    if not ensp: return None
    r = _get(f"{ENSEMBL}/xrefs/id/{ensp}", headers=HJSON)
    xs = r.json()
    # buscar display_id o primary_id tipo Qxxxx-#
    for x in xs:
        db = (x.get("dbname") or "").lower()
        disp = x.get("display_id") or x.get("primary_id") or ""
        if "uniprot" in db and "-" in disp:
            return disp
    return None

@lru_cache(maxsize=10000)
def get_ensp_seq(ensp: str):
    if not ensp: return None
    r = _get(f"{ENSEMBL}/sequence/id/{ensp}", headers={"Content-Type":"text/plain"})
    return r.text.strip()

@lru_cache(maxsize=10000)
def uniprot_base_from_ensp(ensp: str):
    if not ensp: return None
    job = S.post(f"{UNI}/idmapping/run",
                 data={"from":"Ensembl_Protein","to":"UniProtKB","ids":ensp}, timeout=30)
    job.raise_for_status()
    job_id = job.json()["jobId"]
    while True:
        st = _get(f"{UNI}/idmapping/status/{job_id}").json()
        if st.get("jobStatus") in (None, "FINISHED"): break
        time.sleep(0.4)
    res = _get(f"{UNI}/idmapping/results/{job_id}", params={"format":"json"}).json()
    for h in res.get("results", []):
        to = h.get("to")
        if isinstance(to, str):
            return to.split("-")[0]
        if isinstance(to, dict):
            acc = to.get("primaryAccession") or to.get("uniProtkbId")
            if acc: return acc.split("-")[0]
    return None

@lru_cache(maxsize=10000)
def uniprot_isoforms(base_acc: str):
    if not base_acc: return tuple()
    r = _get(f"{UNI}/uniprotkb/{base_acc}.json")
    entry = r.json()
    ids = set()
    for c in entry.get("comments", []):
        if (c.get("commentType") or c.get("type")) in ("ALTERNATIVE PRODUCTS","ALTERNATIVE_PRODUCTS"):
            for iso in c.get("isoforms", []):
                for iid in iso.get("isoformIds", []):
                    if re.match(r"^[A-NR-Z0-9]{6,10}-\d+$", iid):
                        ids.add(iid)
    if not ids and base_acc:
        ids.add(f"{base_acc}-1")
    return tuple(sorted(ids))

@lru_cache(maxsize=100000)
def fetch_uniprot_fasta(acc_or_iso: str):
    if not acc_or_iso: return None
    r = S.get(f"{UNI}/uniprotkb/{acc_or_iso}.fasta", timeout=30)
    if r.status_code != 200: return None
    return "".join(ln.strip() for ln in r.text.splitlines() if ln and not ln.startswith(">"))

def enst_to_uniprot_isoform(enst: str, empty_label: str = "") -> str:
    """
    Devuelve 'Qxxxx-#' si encuentra isoforma; si no, devuelve empty_label (por defecto "").
    """
    if not isinstance(enst, str) or not enst.startswith("ENST"):
        return empty_label
    ensp = get_ensp_from_enst(enst)
    if not ensp:
        return empty_label  # transcript no-protein
    # 1) prioriza lo que diga Ensembl 
    iso = ensembl_xrefs_isoform(ensp)
    if iso:
        return iso
    # 2) match x secuencia
    ensp_seq = get_ensp_seq(ensp)
    base = uniprot_base_from_ensp(ensp)
    if not base or not ensp_seq:
        return empty_label
    for iid in uniprot_isoforms(base):
        seq = fetch_uniprot_fasta(iid)
        if seq and seq == ensp_seq:
            return iid
    return empty_label  # si no hubo match exacto, se deja en blanco


# 2) Aplicar al DataFrame 

def add_uniprot_isoform_with_progress(df: pd.DataFrame,
                                      col_in: str = "Transcript ID",
                                      col_out: str = "Uniprot_isoform",
                                      empty_label: str = "") -> pd.DataFrame:
    n = len(df)
    if n == 0:
        df[col_out] = []
        return df

    results = []
    printed = set()
    start = time.time()

    values = df[col_in].astype(str).tolist()
    for i, enst in enumerate(values, 1):
        results.append(enst_to_uniprot_isoform(enst, empty_label=empty_label))

        # progreso en múltiplos de 10%
        pct = int(i * 100 / n)
        step = (pct // 10) * 10
        if step in {10,20,30,40,50,60,70,80,90,100} and step not in printed:
            elapsed = time.time() - start
            eta = elapsed / i * (n - i)
            def _fmt(sec):
                m, s = divmod(int(sec), 60)
                return f"{m:02d}:{s:02d}"
            print(f"{step}% ({i}/{n}) - elapsed { _fmt(elapsed) } - ETA { _fmt(eta) }")
            printed.add(step)

        # opcional: ligera pausa cada 200 requests para ser amable con APIs
        if i % 200 == 0:
            time.sleep(0.2)

    df[col_out] = results
    # pequeño resumen
    found = sum(bool(x) for x in results)
    print(f"Listo. Total filas: {n} | con isoforma: {found} | en blanco: {n-found}")
    return df

# =========================
# 3) Ejecutar sobre tu tf_2
# =========================
tf_2 = pd.read_csv("/home/shey/Documentos/seqlab/TF/TF_completeIDs.csv")
tf_2 = add_uniprot_isoform_with_progress(tf_2, col_in="Transcript ID", col_out="Uniprot_isoform", empty_label="")

# tf_2.to_csv("home/shey/Documentos/seqlab/TF/tf_2_con_uniprot.csv", index=False)


