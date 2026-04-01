#!/usr/bin/env python3
"""Fix GPL570 and GPL96 platform annotation downloads."""
import requests
import gzip
import shutil
import os

pdir = "/Users/ambuj.mehta/code/Git_Repos/Bio_Claude_testing/platforms"
os.makedirs(pdir, exist_ok=True)

# Correct URL pattern for short GPL IDs
fixes = {
    "GPL570": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
    "GPL96": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz",
}

for gpl, url in fixes.items():
    print(f"Downloading {gpl}...")
    dest_gz = os.path.join(pdir, f"{gpl}.annot.gz")
    dest = os.path.join(pdir, f"{gpl}.annot")

    if os.path.exists(dest):
        print(f"  Already exists: {dest}")
        continue

    try:
        r = requests.get(url, timeout=120)
        r.raise_for_status()
        with open(dest_gz, "wb") as f:
            f.write(r.content)
        with gzip.open(dest_gz, "rb") as fi, open(dest, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        size_mb = os.path.getsize(dest) / (1024 * 1024)
        print(f"  OK: {size_mb:.1f} MB")
    except Exception as e:
        print(f"  annot failed: {e}")
        # Try SOFT file instead
        soft_url = url.replace("/annot/", "/soft/").replace(".annot.gz", "_family.soft.gz")
        print(f"  Trying SOFT: {soft_url}")
        try:
            r = requests.get(soft_url, timeout=120)
            r.raise_for_status()
            soft_gz = os.path.join(pdir, f"{gpl}_family.soft.gz")
            soft_f = os.path.join(pdir, f"{gpl}_family.soft")
            with open(soft_gz, "wb") as f:
                f.write(r.content)
            with gzip.open(soft_gz, "rb") as fi, open(soft_f, "wb") as fo:
                shutil.copyfileobj(fi, fo)
            size_mb = os.path.getsize(soft_f) / (1024 * 1024)
            print(f"  SOFT OK: {size_mb:.1f} MB")
        except Exception as e2:
            print(f"  SOFT also failed: {e2}")

print("\nPlatform files:", os.listdir(pdir))
