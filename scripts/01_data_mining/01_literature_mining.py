#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import dataclasses
import datetime as dt
import logging
import os
import re
import sys
import textwrap
from typing import Iterable


LOGGER = logging.getLogger("literature_mining")


DEFAULT_QUERY = (
    '("acupuncture"[Title/Abstract] OR "electroacupuncture"[Title/Abstract]) '
    'AND ("migraine"[Title/Abstract] OR "headache"[Title/Abstract]) '
    'AND ("mechanism"[Title/Abstract] OR "pathway"[Title/Abstract] OR "neurotransmitter"[Title/Abstract] '
    'OR "neuropeptide"[Title/Abstract] OR "cytokine"[Title/Abstract]) '
    'AND ("2020"[Date - Publication] : "3000"[Date - Publication])'
)


MOLECULE_PATTERNS: list[tuple[str, str, re.Pattern[str]]] = [
    ("CGRP", "neuropeptide", re.compile(r"\bCGRP\b|\bcalcitonin gene-related peptide\b", re.I)),
    ("PACAP", "neuropeptide", re.compile(r"\bPACAP\b|\bpituitary adenylate cyclase-activating polypeptide\b", re.I)),
    ("Substance P", "neuropeptide", re.compile(r"\bsubstance\s*p\b|\bSP\b", re.I)),
    ("NPY", "neuropeptide", re.compile(r"\bNPY\b|\bneuropeptide\s*y\b", re.I)),
    ("5-HT", "neurotransmitter", re.compile(r"\b5-HT\b|\bserotonin\b", re.I)),
    ("Dopamine", "neurotransmitter", re.compile(r"\bdopamine\b|\bDA\b", re.I)),
    ("Norepinephrine", "neurotransmitter", re.compile(r"\bnorepinephrine\b|\bnoradrenaline\b|\bNE\b", re.I)),
    ("GABA", "neurotransmitter", re.compile(r"\bGABA\b|\bgamma-aminobutyric acid\b", re.I)),
    ("Glutamate", "neurotransmitter", re.compile(r"\bglutamate\b|\bGlu\b", re.I)),
    ("ET-1", "peptide", re.compile(r"\bET-?1\b|\bendothelin-?1\b", re.I)),
    ("Nitric oxide", "small_molecule", re.compile(r"(?i:\bnitric oxide\b)|\bNO\b")),
    ("TRPV1", "receptor", re.compile(r"\bTRPV1\b|\btransient receptor potential vanilloid(?:\s+subfamily)?\s+member\s+1\b", re.I)),
    ("BDNF", "protein", re.compile(r"\bBDNF\b|\bbrain-derived neurotrophic factor\b", re.I)),
    ("IL-17", "cytokine", re.compile(r"\bIL-?17(?:A|F)?\b|\binterleukin-?17(?:A|F)?\b", re.I)),
    ("TNF", "cytokine", re.compile(r"\bTNF-?alpha\b|\bTNF\b|\btumor necrosis factor\b", re.I)),
    ("IL-1β", "cytokine", re.compile(r"\bIL-?1\s*β\b|\bIL-?1beta\b|\binterleukin-?1\b", re.I)),
    ("IL-6", "cytokine", re.compile(r"\bIL-?6\b|\binterleukin-?6\b", re.I)),
    ("IL-10", "cytokine", re.compile(r"\bIL-?10\b|\binterleukin-?10\b", re.I)),
]


DIRECTION_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("increased", re.compile(r"\b(increase|increased|elevated|elevation|higher|upregulated|up-regulated)\b", re.I)),
    ("decreased", re.compile(r"\b(decrease|decreased|reduced|reduction|lower|downregulated|down-regulated)\b", re.I)),
]


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def split_sentences(text: str) -> list[str]:
    if not text:
        return []
    normalized = re.sub(r"\s+", " ", text.strip())
    if not normalized:
        return []
    parts = re.split(r"(?<=[.!?])\s+", normalized)
    return [part.strip() for part in parts if part.strip()]


def infer_direction(sentence: str) -> str:
    matched: list[str] = []
    for label, pattern in DIRECTION_PATTERNS:
        if pattern.search(sentence):
            matched.append(label)
    if not matched:
        return "unknown"
    if len(set(matched)) > 1:
        return "mixed"
    return matched[0]


@dataclasses.dataclass(frozen=True)
class PubMedRecord:
    pmid: str
    title: str
    journal: str
    year: str
    abstract: str


def _require_biopython() -> None:
    try:
        import Bio  # noqa: F401
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(
            "Biopython is required. Install with: pip3 install biopython"
        ) from exc


def pubmed_search(*, query: str, retmax: int, email: str, api_key: str | None) -> list[str]:
    _require_biopython()
    from Bio import Entrez  # type: ignore

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    LOGGER.info("Searching PubMed")
    LOGGER.info("Query: %s", query)
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax, usehistory="y")
    result = Entrez.read(handle)
    handle.close()

    id_list = result.get("IdList", [])
    if not id_list:
        LOGGER.warning("No PubMed results for query.")
        return []

    LOGGER.info("PubMed hits: %s", result.get("Count", "?"))
    return list(id_list)


def pubmed_fetch(*, pmids: list[str], email: str, api_key: str | None) -> list[PubMedRecord]:
    _require_biopython()
    from Bio import Entrez  # type: ignore

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    records: list[PubMedRecord] = []
    batch_size = 100
    for start in range(0, len(pmids), batch_size):
        batch = pmids[start : start + batch_size]
        LOGGER.info("Fetching PubMed records: %s-%s / %s", start + 1, min(start + batch_size, len(pmids)), len(pmids))
        handle = Entrez.efetch(db="pubmed", id=",".join(batch), rettype="medline", retmode="text")
        medline_text = handle.read()
        handle.close()

        records.extend(parse_medline(medline_text))
    return records


def parse_medline(medline_text: str) -> list[PubMedRecord]:
    _require_biopython()
    from Bio import Medline  # type: ignore
    import io

    parsed: list[PubMedRecord] = []
    for record in Medline.parse(io.StringIO(medline_text)):
        pmid = str(record.get("PMID", "")).strip()
        title = str(record.get("TI", "")).strip()
        journal = str(record.get("JT", "")).strip()
        abstract = str(record.get("AB", "")).strip()

        year = ""
        dp = str(record.get("DP", "")).strip()
        if dp:
            match = re.search(r"\b(19|20)\d{2}\b", dp)
            if match:
                year = match.group(0)

        if not pmid:
            continue
        parsed.append(PubMedRecord(pmid=pmid, title=title, journal=journal, year=year, abstract=abstract))
    return parsed


@dataclasses.dataclass(frozen=True)
class MoleculeHit:
    pmid: str
    year: str
    title: str
    journal: str
    molecule: str
    molecule_type: str
    direction: str
    evidence_sentence: str


def extract_hits(records: Iterable[PubMedRecord]) -> list[MoleculeHit]:
    hits: list[MoleculeHit] = []

    for record in records:
        if not record.abstract:
            continue

        for sentence in split_sentences(record.abstract):
            for molecule_name, molecule_type, pattern in MOLECULE_PATTERNS:
                if not pattern.search(sentence):
                    continue
                direction = infer_direction(sentence)
                hits.append(
                    MoleculeHit(
                        pmid=record.pmid,
                        year=record.year,
                        title=record.title,
                        journal=record.journal,
                        molecule=molecule_name,
                        molecule_type=molecule_type,
                        direction=direction,
                        evidence_sentence=sentence,
                    )
                )

    return hits


def write_csv(
    *,
    output_path: str,
    hits: list[MoleculeHit],
    query: str,
    retrieved_at: str,
) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fieldnames = [
        "pmid",
        "year",
        "title",
        "journal",
        "molecule",
        "molecule_type",
        "direction",
        "evidence_sentence",
        "query",
        "retrieved_at",
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as file_handle:
        writer = csv.DictWriter(file_handle, fieldnames=fieldnames)
        writer.writeheader()
        for hit in hits:
            writer.writerow(
                {
                    "pmid": hit.pmid,
                    "year": hit.year,
                    "title": hit.title,
                    "journal": hit.journal,
                    "molecule": hit.molecule,
                    "molecule_type": hit.molecule_type,
                    "direction": hit.direction,
                    "evidence_sentence": hit.evidence_sentence,
                    "query": query,
                    "retrieved_at": retrieved_at,
                }
            )

    LOGGER.info("Wrote %s rows to %s", len(hits), output_path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="PubMed literature mining for acupuncture–migraine mechanistic mediators (2020+).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """
            Notes:
              - Provide --email (or set ENTREZ_EMAIL) for NCBI Entrez.
              - Use --api-key (or set ENTREZ_API_KEY) for higher rate limits.
            """
        ).strip(),
    )
    parser.add_argument("--query", default=DEFAULT_QUERY, help="PubMed query (default: a reproducible 2020+ query).")
    parser.add_argument("--retmax", type=int, default=200, help="Maximum number of PubMed records to fetch.")
    parser.add_argument(
        "--output",
        default="data/processed/literature_compounds.csv",
        help="Output CSV path (default: data/processed/literature_compounds.csv)",
    )
    parser.add_argument("--email", default=os.environ.get("ENTREZ_EMAIL", ""), help="Entrez email (or ENV ENTREZ_EMAIL).")
    parser.add_argument(
        "--api-key", default=os.environ.get("ENTREZ_API_KEY", ""), help="Entrez API key (or ENV ENTREZ_API_KEY)."
    )
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging.")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    configure_logging(args.verbose)

    email = (args.email or "").strip()
    api_key = (args.api_key or "").strip() or None
    if not email:
        LOGGER.error("Missing Entrez email. Provide --email or set ENTREZ_EMAIL.")
        return 2

    retrieved_at = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    pmids = pubmed_search(query=args.query, retmax=args.retmax, email=email, api_key=api_key)
    if not pmids:
        LOGGER.warning("No PMIDs to fetch; output will be created with only a header.")
        write_csv(output_path=args.output, hits=[], query=args.query, retrieved_at=retrieved_at)
        return 0

    records = pubmed_fetch(pmids=pmids, email=email, api_key=api_key)
    LOGGER.info("Fetched %s records", len(records))

    hits = extract_hits(records)
    write_csv(output_path=args.output, hits=hits, query=args.query, retrieved_at=retrieved_at)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
