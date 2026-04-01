#!/usr/bin/env python3
"""Download GPL platform annotation files for probe-to-gene mapping."""
import os
import gzip
import shutil
import requests
import time

PLATFORMS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "platforms")
os.makedirs(PLATFORMS_DIR, exist_ok=True)

PLATFORMS = {
    "GPL570": "Affymetrix Human Genome U133 Plus 2.0",
    "GPL1261": "Affymetrix Mouse Genome 430 2.0",
    "GPL8321": "Affymetrix Mouse Gene 1.0 ST",
    "GPL96": "Affymetrix Human Genome U133A",
}

for gpl_id, desc in PLATFORMS.items():
    print(f"\nDownloading {gpl_id} ({desc})...")
    # GEO SOFT annotation file
    prefix = gpl_id.replace("GPL", "")
    nnn = prefix[:-3] + "nnn" if len(prefix) > 3 else prefix + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl_id[:3]}{nnn}/{gpl_id}/annot/{gpl_id}.annot.gz"
    dest_gz = os.path.join(PLATFORMS_DIR, f"{gpl_id}.annot.gz")
    dest_txt = os.path.join(PLATFORMS_DIR, f"{gpl_id}.annot")

    if os.path.exists(dest_txt):
        print(f"  Already exists: {dest_txt}")
        continue

    try:
        r = requests.get(url, timeout=120, stream=True)
        r.raise_for_status()
        with open(dest_gz, "wb") as f:
            for chunk in r.iter_content(8192):
                f.write(chunk)
        size_mb = os.path.getsize(dest_gz) / (1024 * 1024)
        print(f"  Downloaded: {size_mb:.1f} MB")

        # Decompress
        with gzip.open(dest_gz, "rb") as f_in, open(dest_txt, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"  Decompressed: {dest_txt}")
    except Exception as e:
        print(f"  Error: {e}")
        # Try alternative: SOFT file
        alt_url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl_id[:3]}{nnn}/{gpl_id}/soft/{gpl_id}_family.soft.gz"
        print(f"  Trying SOFT file: {alt_url}")
        try:
            r = requests.get(alt_url, timeout=120, stream=True)
            r.raise_for_status()
            soft_gz = os.path.join(PLATFORMS_DIR, f"{gpl_id}_family.soft.gz")
            with open(soft_gz, "wb") as f:
                for chunk in r.iter_content(8192):
                    f.write(chunk)
            with gzip.open(soft_gz, "rb") as f_in, open(os.path.join(PLATFORMS_DIR, f"{gpl_id}_family.soft"), "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            print(f"  SOFT file downloaded and decompressed")
        except Exception as e2:
            print(f"  SOFT download also failed: {e2}")

    time.sleep(1)

print("\nDone!")
print("Files:", os.listdir(PLATFORMS_DIR))
