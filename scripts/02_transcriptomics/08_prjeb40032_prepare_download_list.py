#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import logging
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import pandas as pd

LOG = logging.getLogger(__name__)


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _read_tsv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype="string")
    if df.empty:
        raise ValueError(f"Empty runlist TSV: {path}")
    return df


def _write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def _write_text_lines(lines: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for line in lines:
            f.write(line.rstrip("\n") + "\n")


def _fetch_json(url: str, timeout_s: int = 60, retries: int = 4, backoff_s: float = 1.2) -> dict:
    last_exc: Optional[Exception] = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout_s) as r:
                payload = r.read().decode("utf-8", errors="replace")
            return json.loads(payload)
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_exc = exc
            sleep_s = backoff_s * (2**attempt)
            LOG.warning("Fetch failed (attempt %d/%d) url=%s err=%s; sleep=%.1fs", attempt + 1, retries, url, exc, sleep_s)
            time.sleep(sleep_s)
    raise RuntimeError(f"Failed to fetch JSON after {retries} attempts: {url} ({last_exc})")


def _biosample_label(sample_accession: str) -> tuple[str, str]:
    url = f"https://www.ebi.ac.uk/biosamples/samples/{sample_accession}"
    obj = _fetch_json(url)
    ch = obj.get("characteristics", {}) or {}
    title = ((ch.get("title") or [{}])[0].get("text") or "").strip()
    desc = ((ch.get("description") or [{}])[0].get("text") or "").strip()
    return title, desc


def _normalize_label(title: str, desc: str) -> str:
    label = (title or desc or "").strip()
    if not label:
        return ""
    return " ".join(label.split())


def _label_to_group_phase(label: str) -> tuple[str, str]:
    lab = label.lower()
    if "healthy control" in lab or lab == "control":
        return "control", "basal"
    if "migraine" in lab and "interictal" in lab:
        return "case", "interictal"
    if "migraine" in lab and "ictal" in lab:
        return "case", "ictal"
    # fallback: if it says migraine but no phase, keep as case/unknown
    if "migraine" in lab:
        return "case", "unknown"
    return "unknown", "unknown"


@dataclass(frozen=True)
class Outputs:
    labeled_runlist: Path
    sample_metadata: Path
    download_tsv: Path
    download_urls: Path


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare a labeled runlist + download manifest for ENA BioProject PRJEB40032.",
    )
    p.add_argument(
        "--runlist",
        required=True,
        type=Path,
        help="ENA read_run TSV (run_accession, sample_accession, fastq_ftp, fastq_md5, ...).",
    )
    p.add_argument(
        "--out-labeled-runlist",
        required=True,
        type=Path,
        help="Output TSV with biosample labels and derived group/phase.",
    )
    p.add_argument(
        "--out-sample-metadata",
        required=True,
        type=Path,
        help="Output TSV with one row per sample_accession.",
    )
    p.add_argument(
        "--out-download-tsv",
        required=True,
        type=Path,
        help="Output TSV listing filtered runs to download (includes https_url).",
    )
    p.add_argument(
        "--out-download-urls",
        required=True,
        type=Path,
        help="Output text file with one https URL per line (filtered).",
    )
    p.add_argument(
        "--include-phases",
        default="interictal",
        help="Comma-separated phases to include for download (e.g., interictal,ictal). Controls are always included.",
    )
    p.add_argument("--sleep-between-requests", type=float, default=0.25)
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(message)s",
    )

    include_phases = {p.strip().lower() for p in (args.include_phases or "").split(",") if p.strip()}
    if not include_phases:
        include_phases = {"interictal"}

    df = _read_tsv(args.runlist)
    required = {"run_accession", "sample_accession", "fastq_ftp", "fastq_md5"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Runlist missing required columns {sorted(missing)}: {args.runlist}")

    df = df.copy()
    df["sample_accession"] = df["sample_accession"].astype(str).str.strip()
    df["run_accession"] = df["run_accession"].astype(str).str.strip()
    df["fastq_ftp"] = df["fastq_ftp"].astype(str).str.strip()
    df["fastq_md5"] = df["fastq_md5"].astype(str).str.strip()

    samples = sorted({s for s in df["sample_accession"].tolist() if s and s != "nan"})
    LOG.info("Unique samples: %d", len(samples))

    rows_meta = []
    label_cache: dict[str, tuple[str, str, str, str]] = {}
    for i, sample in enumerate(samples, start=1):
        title, desc = _biosample_label(sample)
        label = _normalize_label(title, desc)
        group, phase = _label_to_group_phase(label)
        label_cache[sample] = (title, desc, label, phase)
        rows_meta.append(
            {
                "sample_accession": sample,
                "label": label,
                "group": group,
                "phase": phase,
                "title": title,
                "description": desc,
                "retrieved_at": _utc_now_iso(),
            }
        )
        if args.sleep_between_requests > 0:
            time.sleep(args.sleep_between_requests)
        if i % 10 == 0:
            LOG.info("Fetched biosamples: %d/%d", i, len(samples))

    meta = pd.DataFrame(rows_meta)
    if meta.empty:
        raise RuntimeError("Failed to build biosample metadata.")

    df["label"] = df["sample_accession"].map(lambda s: label_cache.get(s, ("", "", "", ""))[2])
    df["phase"] = df["sample_accession"].map(lambda s: label_cache.get(s, ("", "", "", ""))[3])
    df["group"] = df["label"].map(lambda lab: _label_to_group_phase(lab)[0] if isinstance(lab, str) else "unknown")
    df["https_url"] = df["fastq_ftp"].map(lambda p: f"https://{p}" if p and not p.startswith("http") else p)
    df["retrieved_at"] = _utc_now_iso()

    # Download subset: controls + specified phases for cases
    want_phase = {p for p in include_phases}
    subset = df[(df["group"] == "control") | ((df["group"] == "case") & (df["phase"].str.lower().isin(want_phase)))].copy()

    LOG.info(
        "Subset rows=%d (controls=%d; cases_included=%d; phases=%s)",
        len(subset),
        int((subset["group"] == "control").sum()),
        int((subset["group"] == "case").sum()),
        ",".join(sorted(want_phase)),
    )

    _write_tsv(df, args.out_labeled_runlist)
    _write_tsv(meta, args.out_sample_metadata)
    _write_tsv(subset, args.out_download_tsv)
    _write_text_lines(subset["https_url"].dropna().astype(str).tolist(), args.out_download_urls)

    LOG.info("Wrote labeled runlist: %s", args.out_labeled_runlist)
    LOG.info("Wrote sample metadata: %s", args.out_sample_metadata)
    LOG.info("Wrote download TSV: %s", args.out_download_tsv)
    LOG.info("Wrote download URLs: %s", args.out_download_urls)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

