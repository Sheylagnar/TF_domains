import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

# 1) Aligner ---------------------------------------------------------------
MATRIX = substitution_matrices.load("BLOSUM62")


def make_aligner(gap_open: float, gap_extend: float, mode: str = "global"):
    al = Align.PairwiseAligner()
    al.substitution_matrix = MATRIX
    al.mode = mode
    al.open_gap_score = -gap_open
    al.extend_gap_score = -gap_extend
    # solo 1 alineamiento óptimo (si la versión lo soporta)
    try:
        al.maximum_number_of_alignments = 1
    except AttributeError:
        pass
    return al


# 2) Concatenar DBDs por isoforma **y por Pfam** ---------------------------
def concat_dbd_per_isoform_by_pfam(df_refs: pd.DataFrame) -> dict:
    """
    Devuelve { (isoform_name, pfam): (DBD_concat, n_segmentos) }.
    Respeta dbd_idx y filtra NaN/''.
    """
    pres = ~(
        df_refs["DBD_seq"].isna()
        | df_refs["DBD_seq"].astype(str).str.strip().isin(["", "nan", "NaN", "-"])
        | df_refs["Pfam"].isna()
        | df_refs["Pfam"].astype(str).str.strip().isin(["", "nan", "NaN", "-"])
    )
    out = {}
    if not pres.any():
        return out

    dfv = df_refs.loc[pres].copy()
    for (iso, pfam), sub in dfv.groupby(["Isoform name", "Pfam"]):
        sub = sub.sort_values("dbd_idx", na_position="last")
        segs = [
            s
            for s in sub["DBD_seq"].astype(str).str.strip()
            if s not in ("", "nan", "NaN", "-")
        ]
        if segs:
            out[(iso, pfam)] = ("".join(segs), len(segs))
    return out


# 3) Alinear y calcular métricas (+ FASTA oneline siempre) -----------------
def align_and_metrics(ref_seq: str, tgt_seq: str, aligner):
    """
    Devuelve:
      score, pct_identity, pct_positives, coverage_ref_pct,
      target_start, target_end, alignment_fasta_oneline.
    """
    it = aligner.align(ref_seq, tgt_seq)
    try:
        aln = next(it)  # solo el mejor
    except StopIteration:
        return None

    ref_blocks, tgt_blocks = aln.aligned
    t_start = int(tgt_blocks[0][0] + 1)
    t_end = int(tgt_blocks[-1][1])

    matches = positives = ref_no_gaps = 0
    for (rs, re), (ts, te) in zip(ref_blocks, tgt_blocks):
        r_seg = ref_seq[rs:re]
        t_seg = tgt_seq[ts:te]
        ref_no_gaps += re - rs
        for r, t in zip(r_seg, t_seg):
            if r == t:
                matches += 1
                positives += 1
            elif MATRIX.get((r, t), MATRIX.get((t, r), -1)) > 0:
                positives += 1

    aln_fasta = aln.format("fasta").replace("\n", "\\n")

    return {
        "score": aln.score,
        "pct_identity": round(100.0 * matches / ref_no_gaps, 2) if ref_no_gaps else 0.0,
        "pct_positives": (
            round(100.0 * positives / ref_no_gaps, 2) if ref_no_gaps else 0.0
        ),
        "coverage_ref_pct": (
            round(100.0 * ref_no_gaps / len(ref_seq), 2) if ref_seq else 0.0
        ),
        "target_start": t_start,
        "target_end": t_end,
        "alignment_fasta_oneline": aln_fasta,
    }


# 4) Orquestador por gen (targets y refs = todas las isoformas con Sequence aa) ----
def align_dbds_samegene_stream(
    df: pd.DataFrame,
    *,
    mode: str = "global",
    single_gaps=(10.0, 0.5),
    multi_gaps=(6.0, 0.2),
    min_ref_len: int = 8,
    include_seqs: bool = True,
) -> pd.DataFrame:
    """
    Por gen ('Gene name (HGNC)'):
      - Targets: **todas** las isoformas con 'Sequence aa' no vacío.
      - Refs por familia Pfam: para **todas** las isoformas con 'Sequence aa' no vacío.
        * Si la ref tiene DBDs de esa familia -> usa su concatenado.
        * Si no tiene -> usa un fallback (el DBD concatenado más largo de esa familia en el gen).
    Devuelve una fila por (Pfam, target_isoform, ref_isoform), incluso si la ref no
    posee esa familia (gracias al fallback).
    """
    al_single = make_aligner(*single_gaps, mode=mode)
    al_multi = make_aligner(*multi_gaps, mode=mode)

    rows = []
    for gene, df_gene in df.groupby("Gene name (HGNC)"):
        if pd.isna(gene) or str(gene).strip() == "":
            continue

        # --- isoformas con Sequence aa: serán targets y refs
        seq_ok_mask = ~(
            df_gene["Sequence aa"].isna()
            | df_gene["Sequence aa"]
            .astype(str)
            .str.strip()
            .isin(["", "nan", "NaN", "-"])
        )
        df_with_seq = (
            df_gene[seq_ok_mask]
            .drop_duplicates(subset=["Isoform name"])[
                ["Isoform name", "Sequence aa", "Translation ID"]
            ]
            .rename(
                columns={
                    "Isoform name": "isoform",
                    "Sequence aa": "seq",
                    "Translation ID": "translation_id",
                }
            )
        )
        if df_with_seq.empty:
            continue

        # targets y refs = todas las isoformas con secuencia
        targets_df = df_with_seq.copy()
        all_ref_isoforms = sorted(df_with_seq["isoform"].astype(str).unique().tolist())

        # --- construir mapa de DBDs concatenados por (isoforma, Pfam) dentro del gen
        pres_dbd = ~(
            df_gene["DBD_seq"].isna()
            | df_gene["DBD_seq"].astype(str).str.strip().isin(["", "nan", "NaN", "-"])
            | df_gene["Pfam"].isna()
            | df_gene["Pfam"].astype(str).str.strip().isin(["", "nan", "NaN", "-"])
        )
        if not pres_dbd.any():
            # el gen no tiene DBDs anotados en ninguna isoforma -> no hay con qué alinear
            continue

        # (iso, pfam) -> (concat, n_seg)
        fam_map = {}
        for (iso, pfam), sub in df_gene.loc[pres_dbd].groupby(["Isoform name", "Pfam"]):
            sub = sub.sort_values("dbd_idx", na_position="last")
            segs = [
                s
                for s in sub["DBD_seq"].astype(str).str.strip()
                if s not in ("", "nan", "NaN", "-")
            ]
            if segs:
                fam_map[(str(iso), str(pfam))] = ("".join(segs), len(segs))

        if not fam_map:
            continue

        # familias presentes en el gen
        pfams_in_gene = sorted({pfam for (_, pfam) in fam_map.keys()})

        # precomputar fallback por familia (más largo)
        fallback_by_pfam = {}
        for pf in pfams_in_gene:
            candidates = [
                (iso, seq_nseg) for (iso, p), seq_nseg in fam_map.items() if p == pf
            ]
            if not candidates:
                continue
            fb_iso, (fb_seq, fb_nseg) = max(candidates, key=lambda kv: len(kv[1][0]))
            fallback_by_pfam[pf] = (fb_iso, fb_seq, fb_nseg)

        # --- crear filas
        for _, trow in targets_df.iterrows():
            t_iso = str(trow["isoform"])
            t_seq = str(trow["seq"])
            t_tid = trow["translation_id"]

            for pf in pfams_in_gene:
                # si no hay fallback para esta familia (raro), saltamos
                if pf not in fallback_by_pfam:
                    continue
                fb_iso, fb_seq, fb_nseg = fallback_by_pfam[pf]

                for r_iso in all_ref_isoforms:
                    key = (str(r_iso), pf)
                    if key in fam_map:
                        r_seq, n_seg = fam_map[key]
                        used_iso = r_iso
                        ref_has_family = True
                    else:
                        r_seq, n_seg = fb_seq, fb_nseg
                        used_iso = fb_iso
                        ref_has_family = False

                    if len(r_seq) < min_ref_len:
                        continue

                    aligner = al_single if n_seg == 1 else al_multi
                    met = align_and_metrics(r_seq, t_seq, aligner)

                    row = {
                        "Gene name (HGNC)": gene,
                        "pfam_family": pf,
                        "target_isoform": t_iso,
                        "ref_isoform": r_iso,  # par solicitado
                        "ref_family_source_isoform": used_iso,  # de dónde salió el DBD usado
                        "ref_has_family": ref_has_family,  # la ref realmente tiene esa familia?
                        "n_dbd_ref": n_seg,
                        "len_ref": len(r_seq),
                        "len_target": len(t_seq),
                        "target_translation_id": t_tid,
                        "ref_seq_used": r_seq if include_seqs else None,
                        "target_seq": t_seq if include_seqs else None,
                    }
                    row.update(
                        met
                        if met
                        else {
                            "score": None,
                            "pct_identity": None,
                            "pct_positives": None,
                            "coverage_ref_pct": None,
                            "target_start": None,
                            "target_end": None,
                            "alignment_fasta_oneline": None,
                        }
                    )
                    rows.append(row)

    # ordenar y devolver
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(
        ["Gene name (HGNC)", "pfam_family", "target_isoform", "ref_isoform", "score"],
        ascending=[True, True, True, True, False],
        na_position="last",
    )


# === Ejemplo de uso ===
df_long = align_dbds_samegene_stream(
    df,
    mode="global",
    single_gaps=(10.0, 0.5),
    multi_gaps=(6.0, 0.2),
    min_ref_len=8,
    include_seqs=True,  # guarda ref_seq, target_seq
)
