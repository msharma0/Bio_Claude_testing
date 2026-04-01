#!/usr/bin/env python3
"""Download E-MEXP-2280 from EBI BioStudies with correct paths."""
import requests
import os
import json

base = "https://ftp.ebi.ac.uk/biostudies/fire/E-MEXP-/280/E-MEXP-2280"
out_dir = "/Users/ambuj.mehta/code/Git_Repos/Bio_Claude_testing/neurodeg_datasets/E-MEXP-2280"
os.makedirs(out_dir, exist_ok=True)

# Get file list from BioStudies API
resp = requests.get("https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MEXP-2280", timeout=60)
study = resp.json()


def find_files(obj, files=None):
    if files is None:
        files = []
    if isinstance(obj, dict):
        if "path" in obj and obj.get("type") == "file":
            files.append(obj["path"])
        for v in obj.values():
            find_files(v, files)
    elif isinstance(obj, list):
        for item in obj:
            find_files(item, files)
    return files


files = find_files(study)
print(f"Found {len(files)} files in study metadata")

# Save full file list
with open(os.path.join(out_dir, "file_list.json"), "w") as f:
    json.dump(files, f, indent=2)

# Download key files (sdrf, idf, processed data)
key_extensions = [".sdrf.txt", ".idf.txt", ".processed", "-processed"]
downloaded = 0
for fpath in files:
    fname = os.path.basename(fpath)
    is_key = any(ext in fname.lower() for ext in key_extensions)
    if not is_key:
        continue

    dest = os.path.join(out_dir, fname)
    print(f"Downloading: {fname}")

    # Try multiple URL patterns
    urls_to_try = [
        f"{base}/Files/{fpath}",
        f"{base}/{fpath}",
        f"{base}/Files/{fname}",
    ]

    success = False
    for url in urls_to_try:
        try:
            r = requests.get(url, timeout=120)
            r.raise_for_status()
            with open(dest, "wb") as f:
                f.write(r.content)
            size_mb = len(r.content) / (1024 * 1024)
            print(f"  Saved ({size_mb:.2f} MB)")
            downloaded += 1
            success = True
            break
        except Exception:
            continue

    if not success:
        print(f"  Failed all URL patterns for {fname}")

# If no processed data found, try downloading the whole study zip
if downloaded == 0:
    print("\nNo individual files found. Trying full study download...")
    zip_url = f"{base}/E-MEXP-2280.zip"
    zip_dest = os.path.join(out_dir, "E-MEXP-2280.zip")
    try:
        r = requests.get(zip_url, timeout=180)
        r.raise_for_status()
        with open(zip_dest, "wb") as f:
            f.write(r.content)
        size_mb = len(r.content) / (1024 * 1024)
        print(f"  Saved zip ({size_mb:.2f} MB)")
        downloaded += 1
    except Exception as e:
        print(f"  Zip download failed: {e}")

print(f"\nDownloaded {downloaded} files for E-MEXP-2280")
print(f"Files in: {out_dir}")
print(f"\nAll files in study: {files[:20]}")
