"""
Download all FASTQ files listed in an ENA run table, with resume and md5 check.

Usage (from the project folder):
    python scripts/fetch_reads.py data/meta/PRJNA1251496_ena_runs.tsv data/raw

- Uses the ENA 'fastq_ftp' links over HTTPS (faster than fasterq-dump).
- 'wget -c' resumes a partly downloaded file if the connection drops.
- Every file is checked against ENA's 'fastq_md5'. Files that already pass are skipped,
  so the script is safe to run again.
- A file at full size but with the wrong md5 was corrupted in transit: it is deleted and
  downloaded again (up to 3 attempts). A smaller file is resumed.
- Writes data/meta/download_log.tsv (file, expected md5, status, time).
"""
import csv
import hashlib
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

MAX_ATTEMPTS = 3  # tries per file before giving up


def md5sum(path, chunk=8 * 1024 * 1024):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python scripts/fetch_reads.py <ena_runs.tsv> <output_dir>")
    table, outdir = Path(sys.argv[1]), Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = Path("data/meta/download_log.tsv")
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with table.open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    jobs = []
    for r in rows:
        urls = [u for u in r["fastq_ftp"].split(";") if u]
        md5s = r["fastq_md5"].split(";")
        sizes = [int(s) for s in r["fastq_bytes"].split(";") if s]
        jobs += list(zip(urls, md5s, sizes))
    print(f"{len(jobs)} files to check/download into {outdir}/")

    failed = 0
    with log_path.open("a", encoding="utf-8") as log:
        for n, (url, expected, size) in enumerate(jobs, 1):
            dest = outdir / Path(url).name
            status = "MD5_FAIL"
            for attempt in range(1, MAX_ATTEMPTS + 1):
                if dest.exists() and md5sum(dest) == expected:
                    status = "ok" if attempt > 1 else "ok_existing"
                    break
                # Full size but wrong checksum = corrupted in transit. Resuming cannot fix it,
                # so delete and download again. A smaller file is a partial download: resume it.
                if dest.exists() and dest.stat().st_size >= size:
                    print(f"[{n}/{len(jobs)}] corrupted, deleting {dest.name}", flush=True)
                    dest.unlink()
                print(f"[{n}/{len(jobs)}] downloading {dest.name} (attempt {attempt})", flush=True)
                subprocess.run(["wget", "-q", "-c", "-P", str(outdir), "https://" + url], check=False)
            else:
                if dest.exists() and md5sum(dest) == expected:
                    status = "ok"
            if status == "MD5_FAIL":
                failed += 1
            print(f"[{n}/{len(jobs)}] {status:<12} {dest.name}", flush=True)
            log.write(f"{dest}\t{expected}\t{status}\t{datetime.now(timezone.utc).isoformat(timespec='seconds')}\n")

    print(f"\nDone. {len(jobs) - failed} OK, {failed} failed.")
    if failed:
        print("Run the same command again: it resumes and retries only the failed files.")


if __name__ == "__main__":
    main()
