TF_domains

A lightweight pipeline for constructing, verifying, and annotating DNA-binding domains (DBDs) across human transcription factor (TF) isoforms. Integrates data from APPRIS, Ensembl, UniProt, CisBP 2.0, and dbPTM.

Pipeline Overview
1. 01_database.py — Build the core APPRIS-based TF dataset

Initializes the main table containing:

Gene IDs, transcript IDs, translation IDs

APPRIS isoform annotations

Baseline metadata for all TFs

2. 02_functions.py — Complete missing identifiers and sequences

Adds:

Full amino-acid sequences

Isoform names

UniProt canonical and isoform IDs

Optional:
02_helper_uniprot_isoforms.py retrieves UniProt isoform IDs but is not required for downstream steps.

3. 03_obtain_coordinates.py — Extract DBD coordinates from CisBP

Extracts:

DBD coordinates

Domain type and family

DBD amino-acid sequences

Coordinates exceeding protein length are filtered out.

4. 04_alignment_dbd.py — Map DBDs onto isoforms lacking annotations

Aligns known DBDs from one isoform to isoforms of the same gene without coordinates.

Two alignment strategies are used:

Local alignment (Smith–Waterman) for short, conserved domains

Constrained global alignment for longer/multi-repeat domains

Outputs inferred coordinates, aligned sequences, and alignment scores.

5. 05_obtain_ptm.py — Integrate post-translational modifications

Adds PTM annotations from dbPTM:

Phosphorylation

Acetylation

Ubiquitination

Methylation

Mapped using translation IDs.
