#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import logging
import math
import os
import re
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

LOG = logging.getLogger(__name__)


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _run(cmd: list[str]) -> None:
    LOG.info("Run: %s", " ".join(map(str, cmd)))
    subprocess.check_call(cmd)


def _ensure_kallisto() -> None:
    try:
        subprocess.check_output(["kallisto", "version"], text=True)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError("kallisto not found; install first (e.g., `brew install kallisto`).") from exc


def _parse_gencode_header_to_gene_symbol(header: str) -> Optional[str]:
    # GENCODE transcript fasta header typically:
    # >ENST...|ENSG...|...|...|TRANSCRIPT_NAME|GENE_SYMBOL|...
    if not header.startswith(">"):
        return None
    token = header[1:].split()[0]
    fields = token.split("|")
    # Prefer gene symbol field when present; in GENCODE this is typically the 6th field.
    if len(fields) >= 6:
        gene = fields[5].strip()
        return gene if gene else None
    if len(fields) >= 5:
        gene = fields[4].strip()
        return gene if gene else None
    return None


def _build_tx_to_gene_map(transcripts_fa_gz: Path) -> dict[str, str]:
    tx_to_gene: dict[str, str] = {}
    with gzip.open(transcripts_fa_gz, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line.strip()
            # kallisto uses the fasta record name (first whitespace-delimited token)
            tx = header[1:].split()[0].strip()
            gene = _parse_gencode_header_to_gene_symbol(header)
            if tx and gene:
                tx_to_gene[tx] = gene
    if not tx_to_gene:
        raise RuntimeError(f"Failed to build transcript->gene map from: {transcripts_fa_gz}")
    return tx_to_gene


def _read_abundance_tsv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"target_id", "tpm"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"abundance.tsv missing columns {sorted(missing)}: {path}")
    df = df[["target_id", "tpm"]].copy()
    df["target_id"] = df["target_id"].astype(str)
    df["tpm"] = pd.to_numeric(df["tpm"], errors="coerce")
    df = df.dropna(subset=["tpm"]).copy()
    return df


def _aggregate_tpm_to_gene(
    abundance: pd.DataFrame,
    tx_to_gene: dict[str, str],
) -> pd.DataFrame:
    abundance = abundance.copy()
    abundance["gene_symbol"] = abundance["target_id"].map(tx_to_gene)
    abundance = abundance.dropna(subset=["gene_symbol"]).copy()
    gene = abundance.groupby("gene_symbol", as_index=False)["tpm"].sum()
    return gene


@dataclass(frozen=True)
class Outputs:
    sample_metadata: Path
    expression_matrix: Path
    run_config: Path


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Quantify PRJEB40032 FASTQs with kallisto and build a gene-level validation matrix.",
    )
    p.add_argument("--manifest", required=True, type=Path, help="TSV (DATA_099) with sample_accession, run_accession, https_url, fastq_md5, ...")
    p.add_argument("--fastq-dir", required=True, type=Path, help="Directory containing downloaded *.fastq.gz files")
    p.add_argument("--kallisto-index", required=True, type=Path, help="kallisto index file (.idx)")
    p.add_argument("--transcripts-fasta", required=True, type=Path, help="GENCODE transcripts fasta.gz used for tx->gene mapping")
    p.add_argument("--outdir", required=True, type=Path, help="Output directory for per-sample kallisto results")
    p.add_argument("--out-expression-matrix", required=True, type=Path, help="Gene-level expression matrix TSV (log2(TPM+1))")
    p.add_argument("--out-sample-metadata", required=True, type=Path, help="Sample metadata TSV aligned to matrix columns")
    p.add_argument("--out-run-config", required=True, type=Path, help="JSON run config for traceability")

    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--fragment-length", type=float, default=200.0, help="Single-end fragment length (assumption if unknown)")
    p.add_argument("--fragment-sd", type=float, default=20.0, help="Single-end fragment length SD")
    p.add_argument("--resume", action="store_true", help="Skip kallisto quant if abundance.tsv already exists for a sample")
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(message)s",
    )
    _ensure_kallisto()

    manifest = pd.read_csv(args.manifest, sep="\t", dtype="string")
    required = {"sample_accession", "run_accession", "file_name", "fastq_md5", "group", "phase"}
    missing = required - set(manifest.columns)
    if missing:
        # Backwards-compatible: infer file_name from https_url if needed
        if "https_url" in manifest.columns and "file_name" in missing:
            manifest["file_name"] = manifest["https_url"].astype(str).map(lambda u: Path(u).name)
            missing = required - set(manifest.columns)
        if missing:
            raise ValueError(f"Manifest missing required columns {sorted(missing)}: {args.manifest}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    args.out_expression_matrix.parent.mkdir(parents=True, exist_ok=True)
    args.out_sample_metadata.parent.mkdir(parents=True, exist_ok=True)
    args.out_run_config.parent.mkdir(parents=True, exist_ok=True)

    LOG.info("Building transcript->gene mapping from: %s", args.transcripts_fasta)
    tx_to_gene = _build_tx_to_gene_map(args.transcripts_fasta)
    LOG.info("Transcript->gene mappings: %d", len(tx_to_gene))

    # Run kallisto quant per sample, then aggregate
    sample_rows = []
    gene_tables = []
    for row in manifest.to_dict(orient="records"):
        sample = str(row["sample_accession"]).strip()
        run = str(row["run_accession"]).strip()
        group = str(row["group"]).strip()
        phase = str(row["phase"]).strip()
        file_name = str(row["file_name"]).strip()
        md5 = str(row["fastq_md5"]).strip()

        fastq_path = args.fastq_dir / file_name
        if not fastq_path.exists():
            raise FileNotFoundError(f"FASTQ not found: {fastq_path}")

        sample_out = args.outdir / sample
        abundance_path = sample_out / "abundance.tsv"
        if args.resume and abundance_path.exists():
            LOG.info("Resume: using existing %s", abundance_path)
        else:
            sample_out.mkdir(parents=True, exist_ok=True)
            cmd = [
                "kallisto",
                "quant",
                "-i",
                str(args.kallisto_index),
                "-o",
                str(sample_out),
                "--threads",
                str(int(args.threads)),
                "--single",
                "-l",
                str(float(args.fragment_length)),
                "-s",
                str(float(args.fragment_sd)),
                str(fastq_path),
            ]
            _run(cmd)

        abundance = _read_abundance_tsv(abundance_path)
        gene = _aggregate_tpm_to_gene(abundance, tx_to_gene)
        gene = gene.rename(columns={"tpm": sample})
        gene_tables.append(gene)

        sample_rows.append(
            {
                "sample_id": sample,
                "run_accession": run,
                "group": group,
                "phase": phase,
                "fastq_file": str(fastq_path),
                "fastq_md5_expected": md5,
                "kallisto_outdir": str(sample_out),
                "abundance_tsv": str(abundance_path),
                "retrieved_at": _utc_now_iso(),
            }
        )

    meta = pd.DataFrame(sample_rows).sort_values(["group", "phase", "sample_id"]).reset_index(drop=True)

    # Merge gene tables on gene_symbol
    merged: Optional[pd.DataFrame] = None
    for gt in gene_tables:
        if merged is None:
            merged = gt
        else:
            merged = merged.merge(gt, on="gene_symbol", how="outer")
    if merged is None:
        raise RuntimeError("No gene tables produced.")

    merged = merged.fillna(0.0)
    # log2(TPM + 1)
    sample_cols = [c for c in merged.columns if c != "gene_symbol"]
    merged[sample_cols] = np.log2(merged[sample_cols].astype(float) + 1.0)
    merged = merged.sort_values("gene_symbol").reset_index(drop=True)

    merged.to_csv(args.out_expression_matrix, sep="\t", index=False)
    meta.to_csv(args.out_sample_metadata, sep="\t", index=False)

    cfg = {
        "created_at": _utc_now_iso(),
        "inputs": {
            "manifest": str(args.manifest),
            "fastq_dir": str(args.fastq_dir),
            "kallisto_index": str(args.kallisto_index),
            "transcripts_fasta": str(args.transcripts_fasta),
        },
        "params": {
            "threads": int(args.threads),
            "fragment_length": float(args.fragment_length),
            "fragment_sd": float(args.fragment_sd),
            "resume": bool(args.resume),
        },
        "outputs": {
            "expression_matrix": str(args.out_expression_matrix),
            "sample_metadata": str(args.out_sample_metadata),
            "kallisto_outdir": str(args.outdir),
        },
        "notes": [
            "Single-end kallisto quant requires fragment length/SD; defaults are assumptions unless provided by the study.",
            "Expression matrix values are log2(TPM+1) aggregated from transcript-level TPM to gene_symbol using GENCODE transcript headers.",
        ],
    }
    args.out_run_config.write_text(json.dumps(cfg, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")

    LOG.info("Wrote expression matrix: %s", args.out_expression_matrix)
    LOG.info("Wrote sample metadata: %s", args.out_sample_metadata)
    LOG.info("Done. genes=%d samples=%d", int(len(merged)), int(len(sample_cols)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
