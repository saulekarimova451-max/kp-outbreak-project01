"""
Fetch run metadata for a BioProject from ENA and summarise it.

Usage (from the project folder):
    python scripts/fetch_metadata.py PRJNA1251496

Writes:
    data/meta/<ACC>_ena_runs.tsv   full run table from ENA (raw, unedited)
    data/meta/<ACC>_summary.tsv    one row per run: platform, coverage, size, match key
    data/meta/provenance.tsv       when and from which URL the table was fetched
Prints:
    runs per platform, total download size, coverage per run,
    and which Illumina and Nanopore runs look like the same isolate.

Needs only the Python standard library (no pip install).
"""
import csv
import io
import re
import sys
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

GENOME_SIZE = 5_500_000  # K. pneumoniae genome, about 5.5 Mb

FIELDS = [
    "run_accession", "sample_accession", "secondary_sample_accession",
    "sample_alias", "sample_title", "experiment_title",
    "instrument_platform", "instrument_model", "library_layout",
    "read_count", "base_count", "fastq_bytes", "fastq_md5", "fastq_ftp",
]

# Words that describe the platform, not the isolate. Removed when matching names.
PLATFORM_WORDS = r"(?<![a-z])(final|trimmed|raw|nanopore|illumina|hiseq|novaseq|miseq|gridion|minion|promethion|ont|ill|long|short|lr|sr)(?![a-z])"


def ena_url(accession):
    params = urllib.parse.urlencode({
        "accession": accession,
        "result": "read_run",
        "fields": ",".join(FIELDS),
        "format": "tsv",
        "limit": 0,  # 0 = all rows
    })
    return "https://www.ebi.ac.uk/ena/portal/api/filereport?" + params


def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read().decode("utf-8")


def total_bytes(value):
    # Paired runs list two files separated by ';'
    return sum(int(x) for x in value.split(";") if x.strip().isdigit())


def match_key(name):
    """Turn a sample name into a key shared by both platforms, e.g. 'KP12_ONT' -> 'kp12'."""
    key = name.lower()
    key = re.sub(PLATFORM_WORDS, "", key)
    return re.sub(r"[^a-z0-9]", "", key)


def summarise(rows):
    out = []
    for r in rows:
        bases = int(r["base_count"] or 0)
        name = r["sample_alias"] or r["sample_title"] or r["sample_accession"]
        out.append({
            "run": r["run_accession"],
            "biosample": r["sample_accession"],
            "sample_name": name,
            "match_key": match_key(name),
            "platform": r["instrument_platform"],
            "model": r["instrument_model"],
            "layout": r["library_layout"],
            "reads": int(r["read_count"] or 0),
            "bases": bases,
            "coverage_x": round(bases / GENOME_SIZE),
            "size_gb": round(total_bytes(r["fastq_bytes"]) / 1e9, 2),
        })
    return out


def report(summary):
    print(f"\nRuns found: {len(summary)}")
    by_platform = defaultdict(list)
    for s in summary:
        by_platform[s["platform"]].append(s)
    for platform, runs in by_platform.items():
        size = sum(s["size_gb"] for s in runs)
        print(f"  {platform:<17} {len(runs):>3} runs   {size:7.2f} GB")
    print(f"  {'TOTAL download':<17}       {sum(s['size_gb'] for s in summary):11.2f} GB")

    print("\nPer run:")
    print(f"  {'run':<13} {'platform':<16} {'layout':<7} {'cov':>5}  {'GB':>5}  sample name")
    for s in sorted(summary, key=lambda s: (s["match_key"], s["platform"])):
        print(f"  {s['run']:<13} {s['platform']:<16} {s['layout']:<7} {s['coverage_x']:>4}x  "
              f"{s['size_gb']:>5}  {s['sample_name']}")

    groups = defaultdict(set)
    for s in summary:
        groups[s["match_key"]].add(s["platform"])
    both = [k for k, p in groups.items() if len(p) > 1]
    print(f"\nIsolates with BOTH Illumina and Nanopore (by name): {len(both)}")
    for k in sorted(both):
        runs = [f"{s['platform'][:3]}={s['run']}" for s in summary if s["match_key"] == k]
        print(f"  {k}: {', '.join(runs)}")
    if not both:
        print("  None matched automatically. Check the 'sample name' column above by eye.")


def main():
    accession = sys.argv[1] if len(sys.argv) > 1 else "PRJNA1251496"
    outdir = Path("data/meta")
    outdir.mkdir(parents=True, exist_ok=True)

    url = ena_url(accession)
    print(f"Fetching {accession} from ENA ...")
    text = fetch(url)
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows:
        sys.exit(f"No runs returned for {accession}. Check the accession.")

    raw_path = outdir / f"{accession}_ena_runs.tsv"
    raw_path.write_text(text, encoding="utf-8")

    summary = summarise(rows)
    sum_path = outdir / f"{accession}_summary.tsv"
    with sum_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    prov = outdir / "provenance.tsv"
    new = not prov.exists()
    with prov.open("a", encoding="utf-8") as f:
        if new:
            f.write("file\tsource_url\tfetched_utc\n")
        f.write(f"{raw_path}\t{url}\t{datetime.now(timezone.utc).isoformat(timespec='seconds')}\n")

    report(summary)
    print(f"\nSaved: {raw_path}\n       {sum_path}\n       {prov}")


if __name__ == "__main__":
    main()
