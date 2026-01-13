#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import logging
import os
import re
import sys
from dataclasses import dataclass, asdict
from typing import Iterable


LOGGER = logging.getLogger("find_human_migraine_datasets")


GSE_RE = re.compile(r"\bGSE\d+\b")
ARRAYEXP_RE = re.compile(r"\bE-MTAB-\d+\b|\bE-GEOD-\d+\b", re.I)
BIOPROJECT_RE = re.compile(r"\bPRJNA\d+\b|\bPRJEB\d+\b|\bPRJDB\d+\b", re.I)
SRA_RE = re.compile(r"\bSRP\d+\b|\bERP\d+\b|\bDRP\d+\b", re.I)


DEFAULT_PUBMED_QUERY = (
    '("migraine"[Title/Abstract] OR "chronic migraine"[Title/Abstract] OR "vestibular migraine"[Title/Abstract]) '
    'AND ("transcriptome"[Title/Abstract] OR "gene expression"[Title/Abstract] OR "RNA-seq"[Title/Abstract] OR microarray[Title/Abstract]) '
    'AND ("2020"[Date - Publication] : "3000"[Date - Publication]) '
    'AND Humans[MeSH Terms]'
)

DEFAULT_GDS_QUERY = (
    '(migraine[All Fields] OR "chronic migraine"[All Fields] OR "vestibular migraine"[All Fields]) '
    'AND "Homo sapiens"[Organism]'
)


@dataclass(frozen=True)
class Candidate:
    source: str
    pmid: str
    year: str
    title: str
    journal: str
    geo_accession: str
    geo_entry_type: str
    geo_gds_type: str
    geo_taxon: str
    geo_n_samples: str
    geo_pubmed_ids: str
    geo_ftp_link: str
    geo_supp_files: str
    accessions: str
    gse_accessions: str
    arrayexpress_accessions: str
    bioproject_accessions: str
    sra_accessions: str
    note: str


def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def require_biopython() -> None:
    try:
        import Bio  # noqa: F401
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("Biopython is required. Install with: pip3 install biopython") from exc


def pubmed_search(query: str, email: str, api_key: str | None, retmax: int) -> list[str]:
    require_biopython()
    from Bio import Entrez  # type: ignore

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    LOGGER.info("PubMed search query: %s", query)
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
    result = Entrez.read(handle)
    handle.close()
    ids = list(result.get("IdList", []))
    LOGGER.info("PubMed hits: %s (returned %s)", result.get("Count", "?"), len(ids))
    return ids


def pubmed_fetch_medline(pmids: list[str], email: str, api_key: str | None) -> str:
    require_biopython()
    from Bio import Entrez  # type: ignore

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    batch_size = 200
    chunks: list[str] = []
    for start in range(0, len(pmids), batch_size):
        batch = pmids[start : start + batch_size]
        handle = Entrez.efetch(db="pubmed", id=",".join(batch), rettype="medline", retmode="text")
        chunks.append(handle.read())
        handle.close()
    return "\n".join(chunks)

def gds_search(query: str, email: str, api_key: str | None, retmax: int) -> list[str]:
    require_biopython()
    from Bio import Entrez  # type: ignore

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    LOGGER.info("GDS search query: %s", query)
    handle = Entrez.esearch(db="gds", term=query, retmax=retmax)
    result = Entrez.read(handle)
    handle.close()
    ids = list(result.get("IdList", []))
    LOGGER.info("GDS hits: %s (returned %s)", result.get("Count", "?"), len(ids))
    return ids


def gds_fetch_summaries(gds_ids: list[str], email: str, api_key: str | None) -> list[dict]:
    require_biopython()
    from Bio import Entrez  # type: ignore

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    batch_size = 200
    docs: list[dict] = []
    for start in range(0, len(gds_ids), batch_size):
        batch = gds_ids[start : start + batch_size]
        handle = Entrez.esummary(db="gds", id=",".join(batch))
        res = Entrez.read(handle)
        handle.close()
        docs.extend(list(res))
    return docs


def _as_str(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        parts: list[str] = []
        for item in value:
            atom = _as_str(item)
            if atom:
                parts.append(atom)
        return "|".join(parts)
    if isinstance(value, (str, bytes)):
        return str(value).strip()
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, (int, float)):
        if isinstance(value, float) and not value.is_integer():
            return str(value)
        return str(int(value))
    try:
        return str(int(value))  # Biopython IntegerElement and friends
    except Exception:
        return str(value).strip()


def parse_medline_records(medline_text: str) -> Iterable[Candidate]:
    require_biopython()
    from Bio import Medline  # type: ignore
    import io

    for rec in Medline.parse(io.StringIO(medline_text)):
        pmid = str(rec.get("PMID", "")).strip()
        title = str(rec.get("TI", "")).strip()
        journal = str(rec.get("JT", "")).strip()
        abstract = str(rec.get("AB", "")).strip()
        dp = str(rec.get("DP", "")).strip()
        year = ""
        match = re.search(r"\b(19|20)\d{2}\b", dp)
        if match:
            year = match.group(0)

        text = f"{title} {abstract}"
        gses = sorted(set(GSE_RE.findall(text)))
        aexp = sorted(set(ARRAYEXP_RE.findall(text)))
        biop = sorted(set(BIOPROJECT_RE.findall(text)))
        sra = sorted(set(SRA_RE.findall(text)))

        accessions = sorted(set(gses + aexp + biop + sra))
        note = ""
        if not accessions:
            note = "No accession in title/abstract; may require full text/supplement search."

        if not pmid:
            continue
        yield Candidate(
            source="pubmed",
            pmid=pmid,
            year=year,
            title=title,
            journal=journal,
            geo_accession="",
            geo_entry_type="",
            geo_gds_type="",
            geo_taxon="",
            geo_n_samples="",
            geo_pubmed_ids="",
            geo_ftp_link="",
            geo_supp_files="",
            accessions="|".join(accessions),
            gse_accessions="|".join(gses),
            arrayexpress_accessions="|".join(aexp),
            bioproject_accessions="|".join(biop),
            sra_accessions="|".join(sra),
            note=note,
        )

def parse_gds_records(gds_docs: list[dict]) -> list[Candidate]:
    def is_allowed_gds_type(gds_type: str) -> bool:
        if not gds_type:
            return False
        allowed_prefixes = ("expression profiling", "non-coding rna profiling")
        parts = [p.strip().lower() for p in gds_type.split(";") if p.strip()]
        return any(any(part.startswith(prefix) for prefix in allowed_prefixes) for part in parts) or (
            any(prefix in gds_type.lower() for prefix in allowed_prefixes)
        )

    rows: list[Candidate] = []
    for doc in gds_docs:
        entry_type = _as_str(doc.get("entryType", ""))
        accession = _as_str(doc.get("Accession", ""))
        taxon = _as_str(doc.get("taxon", ""))
        gds_type = _as_str(doc.get("gdsType", ""))

        if entry_type != "GSE":
            continue
        if taxon and "Homo sapiens" not in taxon:
            continue
        if gds_type and not is_allowed_gds_type(gds_type):
            continue

        n_samples = _as_str(doc.get("n_samples", ""))
        title = _as_str(doc.get("title", "")) or _as_str(doc.get("SeriesTitle", ""))
        series_title = _as_str(doc.get("SeriesTitle", "")) or title
        pubmed_ids = _as_str(doc.get("PubMedIds", ""))
        ftp_link = _as_str(doc.get("FTPLink", ""))
        supp_files = _as_str(doc.get("suppFile", ""))

        note = ""
        if not supp_files and not ftp_link:
            note = "No supp/FTP links in GDS summary; may need manual GEO page check."
        if "migraine" not in (series_title or "").lower() and "migraine" not in (title or "").lower():
            note = (note + " " if note else "") + "Title does not contain 'migraine'; verify relevance."

        rows.append(
            Candidate(
                source="gds",
                pmid="",
                year="",
                title=series_title,
                journal="",
                geo_accession=accession,
                geo_entry_type=entry_type,
                geo_gds_type=gds_type,
                geo_taxon=taxon,
                geo_n_samples=n_samples,
                geo_pubmed_ids=pubmed_ids,
                geo_ftp_link=ftp_link,
                geo_supp_files=supp_files,
                accessions=accession,
                gse_accessions=accession if accession.startswith("GSE") else "",
                arrayexpress_accessions="",
                bioproject_accessions="",
                sra_accessions="",
                note=note,
            )
        )
    return rows


def write_csv(path: str, rows: list[Candidate], query: str, retrieved_at: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fieldnames = (
        list(asdict(rows[0]).keys())
        if rows
        else list(
            asdict(
                Candidate(
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                )
            ).keys()
        )
    )
    # Add provenance columns
    fieldnames += ["query", "retrieved_at"]

    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            d = asdict(row)
            d["query"] = query
            d["retrieved_at"] = retrieved_at
            w.writerow(d)


def main() -> int:
    parser = argparse.ArgumentParser(description="Find candidate human migraine transcriptome studies and extract accessions from PubMed title/abstract.")
    parser.add_argument("--source", choices=["pubmed", "gds", "both"], default="both", help="Candidate source to query.")
    parser.add_argument("--pubmed-query", default=DEFAULT_PUBMED_QUERY, help="PubMed query.")
    parser.add_argument("--gds-query", default=DEFAULT_GDS_QUERY, help="GEO/GDS (db=gds) query.")
    parser.add_argument("--retmax", type=int, default=200, help="Max PubMed records to fetch.")
    parser.add_argument("--email", default=os.environ.get("ENTREZ_EMAIL", ""), help="Entrez email (or ENV ENTREZ_EMAIL).")
    parser.add_argument("--api-key", default=os.environ.get("ENTREZ_API_KEY", ""), help="Entrez API key (optional).")
    parser.add_argument("--out", default="data/processed/human_migraine_dataset_candidates.csv", help="Output CSV path.")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging.")
    args = parser.parse_args()

    configure_logging(args.verbose)
    email = (args.email or "").strip()
    api_key = (args.api_key or "").strip() or None
    if not email:
        LOGGER.error("Missing Entrez email. Set ENTREZ_EMAIL or pass --email.")
        return 2

    retrieved_at = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    rows: list[Candidate] = []
    queries: list[str] = []

    if args.source in {"pubmed", "both"}:
        pmids = pubmed_search(args.pubmed_query, email=email, api_key=api_key, retmax=args.retmax)
        if pmids:
            medline = pubmed_fetch_medline(pmids, email=email, api_key=api_key)
            rows.extend(list(parse_medline_records(medline)))
        queries.append(f"pubmed:{args.pubmed_query}")

    if args.source in {"gds", "both"}:
        gds_ids = gds_search(args.gds_query, email=email, api_key=api_key, retmax=args.retmax)
        if gds_ids:
            docs = gds_fetch_summaries(gds_ids, email=email, api_key=api_key)
            rows.extend(parse_gds_records(docs))
        queries.append(f"gds:{args.gds_query}")

    query_str = " || ".join(queries)
    write_csv(args.out, rows, query=query_str, retrieved_at=retrieved_at)
    LOGGER.info("Wrote %s rows to %s", len(rows), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
