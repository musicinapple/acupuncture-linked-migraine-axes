# Figure 1 Generation Prompt

To ensure full transparency and reproducibility of the visual elements in this study, we provide the specific prompt used to generate **Figure 1 (Study Workflow)**. 

Figure 1 was generated using an AI-based image generation tool based on the parameters and quantitative results of our study. 

## Prompt Text

```text
Task: Create a premium, portrait-oriented (3:4) workflow schematic for a bioinformatics SCI paper (Nature/Cell style). Flat vector only (no gradients, no shadows, no 3D). White background, generous margins, strict grid alignment.

Global style:
- Font: Arial/Helvetica only.
- Text hierarchy: Section headers 11pt bold; box titles 10pt bold; body 8.5–9pt; small numeric badges 7.5pt.
- Line/arrow: 1.0–1.2 pt, straight arrows, consistent spacing.
- Palette (flat): Input/Data = Blue #4C78A8; Analysis = Orange #F28E2B; Interpretation/Output = Green #59A14F; Neutral/Limitations = Gray #B0B0B0.
- Do NOT include “Figure 1”, “Fig.1”, or any figure label/title outside panel letters. No watermarks.

Layout (top→bottom, 3 stages):

STAGE 1 — KNOWLEDGE & DATA HARVESTING (blue header band)
- Top box: “Knowledge‑Driven Target Discovery” (icon: document + database).
- Left box: “PubMed Mining (Endogenous Mediators) + CTD (Migraine Targets)”
  - small note: “Intersection: 4 shared targets (CALCA, CALCB, TAC1, TNF)”
- Right stacked dataset cards (A/B/C, perfectly aligned):
  - A. “E‑MTAB‑13397 — Whole blood RNA‑seq (repeated measures)”
    - badge: “N=104 libraries; 52 subjects”; subtitle: “Interictal migraine vs control”
  - B. “PRJEB40032 — PBMC RNA‑seq (multi‑phase)”
    - badge: “N=42 libraries”
    - sub‑badges: “Cases: interictal(20), ictal(10); Controls: basal(12)”
    - small gray footnote: “Library-level counts; subject IDs reconstructed from ENA sample_alias (ictal: 10 libraries from 8 subjects)”
  - C. “GSE198274 — Serum exosomal miRNA”
    - badge: “N=30 samples”
    - sub‑badges: “10 controls; 10 case-before; 10 case-after (paired cases)”

STAGE 2 — ANALYSIS & MULTI‑OMICS INTEGRATION (orange header band)
- Center hub: “miRNA–mRNA Integration & Pathway Profiling”
- Left module: “Subject‑aware ML Feature Selection (LASSO / RF / SVM‑RFE)”
  - orange pills (must match exactly):
    - “Internal subject-held-out AUC = 0.531”
    - “External AUC = 0.483 (inverted-score AUC = 0.517)”
  - tiny gray footnote: “AUC<0.5 indicates opposite separation under fixed direction”
- Small box into hub: “miRNA baseline signal: 38 miRNAs (FDR<0.05) — case-before vs control”
- Small gray box into hub: “Post-treatment delta: 0 miRNAs (FDR<0.05) — case-after vs case-before”
- Right module: “Two‑sample Mendelian Randomization (Driver check)”
  - dashed gray arrow down: “Nominal signals; FDR>0.05 (no robust drivers)”

STAGE 3 — RESULTS & INTERPRETATION (green header band)
- Output box 1: “32‑Gene Migraine Signature (fixed intersection)”
- Output box 2: “Molecular Subtyping (k=3; cases only)”
  - badge: “General clustering: C1=7, C2=15, C3=8 (libraries)”
  - sub-badge (small): “Acupuncture-linked sensitivity: C1=9, C2=16, C3=5 (libraries)”
- Bottom conclusion box (green border; conservative wording only; no causal claim):
  - “Findings are consistent with state‑dependent molecular heterogeneity and downstream response pathways rather than primary inherited drivers.”

Legend (bottom-right):
- Blue = Input/Data; Orange = Analysis; Green = Interpretation/Output; Dashed gray = Limited/no robust genetic support

Hard constraints (must follow):
- Every number must match exactly: Internal AUC 0.531; External AUC 0.483; inverted 0.517; general subtypes 7/15/8; acupuncture sensitivity subtypes 9/16/5; miRNA baseline 38; post-treatment 0; dataset Ns as specified.
- Do NOT add any new numbers, claims, or dataset IDs.
```
