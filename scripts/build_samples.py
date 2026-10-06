"""
Build the master sample table: one row per isolate, linking its Illumina and Nanopore
runs and BioSamples, plus metadata pulled from NCBI BioSample. Also logs every
identifier or metadata problem found (Task 2: ID reconciliation).

Usage (from the project folder, after fetch_metadata.py):
    python scripts/build_samples.py PRJNA1251496

Reads:
    data/meta/<ACC>_summary.tsv          from fetch_metadata.py
Writes:
    data/meta/biosample/<SAMN>.xml       raw NCBI BioSample records (provenance)
    config/samples.tsv                   master table used by the pipeline (committed)
    report/tables/id_issues.tsv          problems found and how they were resolved (committed)
"""
import csv
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
KEEP = ["strain", "isolate", "collection_date", "isolation_source", "host",
        "host_disease", "geo_loc_name", "collected_by", "sample_name"]


def fetch_biosample(acc, cache_dir):
    """Download one BioSample XML (cached) and return its attributes as a dict."""
    path = cache_dir / f"{acc}.xml"
    if not path.exists():
        url = EUTILS + "?" + urllib.parse.urlencode({"db": "biosample", "id": acc})
        with urllib.request.urlopen(url, timeout=60) as r:
            path.write_bytes(r.read())
        time.sleep(0.4)  # NCBI allows max 3 requests per second without an API key
    root = ET.parse(path).getroot()
    attrs = {"title": root.findtext(".//Description/Title", default="")}
    for a in root.iter("Attribute"):
        key = a.get("harmonized_name") or a.get("attribute_name")
        attrs[key] = (a.text or "").strip()
    return attrs


def main():
    acc = sys.argv[1] if len(sys.argv) > 1 else "PRJNA1251496"
    summary = Path(f"data/meta/{acc}_summary.tsv")
    cache = Path("data/meta/biosample")
    cache.mkdir(parents=True, exist_ok=True)
    Path("report/tables").mkdir(parents=True, exist_ok=True)

    with summary.open(encoding="utf-8") as f:
        runs = list(csv.DictReader(f, delimiter="\t"))

    by_isolate = defaultdict(dict)
    for r in runs:
        tech = "ont" if r["platform"] == "OXFORD_NANOPORE" else "illumina"
        by_isolate[r["match_key"]][tech] = r

    issues = []
    rows = []
    for key in sorted(by_isolate):
        pair = by_isolate[key]
        ill, ont = pair.get("illumina"), pair.get("ont")
        if not ill or not ont:
            issues.append(["missing_platform", key, f"only {list(pair)} present", "isolate kept for available platform only"])
        meta = {}
        for tech, run in (("illumina", ill), ("ont", ont)):
            if run:
                meta[tech] = fetch_biosample(run["biosample"], cache)

        # Name used in all outputs: the shared part of the sample names
        any_run = ill or ont
        isolate = any_run["sample_name"].rsplit("_", 1)[0]

        if ill and ont:
            if ill["biosample"] != ont["biosample"]:
                issues.append(["split_biosample", isolate,
                               f"Illumina {ill['biosample']} vs Nanopore {ont['biosample']}",
                               "linked by shared sample name; both accessions kept in samples.tsv"])
            if ill["sample_name"] != ont["sample_name"]:
                issues.append(["name_suffix", isolate,
                               f"'{ill['sample_name']}' vs '{ont['sample_name']}'",
                               "processing-stage suffixes (_trimmed/_final) stripped for matching"])
            for field in KEEP:
                a, b = meta["illumina"].get(field, ""), meta["ont"].get(field, "")
                if a and b and a != b:
                    issues.append(["attribute_conflict", isolate, f"{field}: Illumina '{a}' vs Nanopore '{b}'",
                                   "Illumina BioSample value used; flagged in limitations"])

        merged = {}
        for field in KEEP:
            value = (meta.get("illumina", {}).get(field) or meta.get("ont", {}).get(field) or "")
            if value.lower() in {"missing", "not collected", "not applicable", "na", "unknown", "restricted access"}:
                value = ""
            if not value and field in ("collection_date", "isolation_source", "host"):
                issues.append(["missing_metadata", isolate, f"{field} is empty or 'missing'", "left blank; noted in report"])
            # 'isolate' is also a BioSample attribute; keep it under another name
            merged["biosample_" + field if field == "isolate" else field] = value

        rows.append({
            "isolate": isolate,
            "illumina_run": ill["run"] if ill else "",
            "illumina_biosample": ill["biosample"] if ill else "",
            "illumina_model": ill["model"] if ill else "",
            "illumina_coverage_x": ill["coverage_x"] if ill else "",
            "ont_run": ont["run"] if ont else "",
            "ont_biosample": ont["biosample"] if ont else "",
            "ont_model": ont["model"] if ont else "",
            "ont_coverage_x": ont["coverage_x"] if ont else "",
            **merged,
        })

    # Isolates that share a numeric prefix (e.g. 0327436-KPC and 0327436-NDM)
    prefixes = defaultdict(list)
    for r in rows:
        prefixes[r["isolate"].split("-")[0]].append(r["isolate"])
    for p, names in prefixes.items():
        if len(names) > 1:
            issues.append(["shared_sample_prefix", p, f"{', '.join(names)} share one sample number",
                           "treated as separate isolates (likely one patient/specimen, two strains); checked by SNP distance"])

    out = Path("config/samples.tsv")
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    iss = Path("report/tables/id_issues.tsv")
    with iss.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["issue_type", "record", "detail", "how_resolved"])
        w.writerows(issues)

    print(f"{len(rows)} isolates written to {out}")
    print(f"{len(issues)} issues written to {iss}\n")
    cols = ["isolate", "illumina_run", "ont_run", "collection_date", "isolation_source", "host"]
    print("\t".join(cols))
    for r in rows:
        print("\t".join(str(r[c]) for c in cols))
    print("\nIssue counts:")
    counts = defaultdict(int)
    for i in issues:
        counts[i[0]] += 1
    for k, v in counts.items():
        print(f"  {k:<22} {v}")


if __name__ == "__main__":
    main()
