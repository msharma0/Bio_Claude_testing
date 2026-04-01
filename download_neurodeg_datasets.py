#!/usr/bin/env python3
"""
Download all 11 neurodegenerative disease datasets from BioGPS/GEO/ArrayExpress.

E-GEOD-* datasets are downloaded from NCBI GEO (as GSE* series).
E-MEXP-* datasets are downloaded from ArrayExpress (EBI).
"""

import os
import sys
import requests
import gzip
import shutil
import time

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "neurodeg_datasets")

# All 11 datasets from BioGPS neurodegenerative disease tag
DATASETS = [
    # (BioGPS ID, Description)
    ("E-GEOD-21942", "Expression data from PBMCs in multiple sclerosis patients and controls"),
    ("E-GEOD-8762",  "Transcription profiling of human lymphocytes from Huntingtons disease patients"),
    ("E-GEOD-26600", "Cycad Genotoxin MAM modulates pathways in cancer and neurodegeneration"),
    ("E-GEOD-63012", "Expression data from brain-specific IkBa conditional knockout mice"),
    ("E-GEOD-52118", "Gene expression in motor pools with differential vulnerability in ALS"),
    ("E-GEOD-23182", "Systemic inflammation modulates Fc receptor expression on microglia"),
    ("E-GEOD-47516", "Gene expression in cerebellum of Cstb knockout mouse"),
    ("E-GEOD-5786",  "Transcription profiling of striata from PGC-1a knockout mice"),
    ("E-GEOD-31458", "Expression data from naive and MPTP-exposed cholinergic transgenic mice"),
    ("E-GEOD-35642", "Transcriptome analysis of chronic in vitro model of Parkinsonism"),
    ("E-MEXP-2280",  "Transcription profiling of patients with four neurodegenerative disorders"),
]


def download_file(url, dest_path, description=""):
    """Download a file with progress indication."""
    print(f"  Downloading: {url}")
    try:
        resp = requests.get(url, stream=True, timeout=120)
        resp.raise_for_status()
        total = int(resp.headers.get("content-length", 0))
        downloaded = 0
        with open(dest_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded += len(chunk)
        size_mb = os.path.getsize(dest_path) / (1024 * 1024)
        print(f"  Saved: {os.path.basename(dest_path)} ({size_mb:.2f} MB)")
        return True
    except Exception as e:
        print(f"  ERROR downloading {url}: {e}")
        return False


def download_geo_dataset(geo_id, dataset_dir):
    """Download a GEO dataset (series matrix + supplementary files)."""
    os.makedirs(dataset_dir, exist_ok=True)

    # GEO FTP base: series matrix file
    # URL pattern: https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSExxxx/matrix/
    prefix = geo_id[:-3] + "nnn" if len(geo_id) > 3 else geo_id + "nnn"

    # Series matrix (main expression data)
    matrix_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{geo_id}/matrix/{geo_id}_series_matrix.txt.gz"
    matrix_path = os.path.join(dataset_dir, f"{geo_id}_series_matrix.txt.gz")

    success = download_file(matrix_url, matrix_path, "series matrix")

    if success:
        # Decompress the .gz file
        txt_path = matrix_path.replace(".gz", "")
        try:
            with gzip.open(matrix_path, "rb") as f_in, open(txt_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            print(f"  Decompressed: {os.path.basename(txt_path)}")
        except Exception as e:
            print(f"  Warning: Could not decompress {matrix_path}: {e}")

    # Also try to download the SOFT file (full metadata)
    soft_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{geo_id}/soft/{geo_id}_family.soft.gz"
    soft_path = os.path.join(dataset_dir, f"{geo_id}_family.soft.gz")
    download_file(soft_url, soft_path, "SOFT file")

    return success


def download_arrayexpress_dataset(ae_id, dataset_dir):
    """Download an ArrayExpress dataset from EBI BioStudies."""
    os.makedirs(dataset_dir, exist_ok=True)

    # Try the BioStudies API to list files
    api_url = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{ae_id}"
    print(f"  Fetching metadata from: {api_url}")

    try:
        resp = requests.get(api_url, timeout=60)
        resp.raise_for_status()
        study = resp.json()

        # Save metadata JSON
        import json
        meta_path = os.path.join(dataset_dir, f"{ae_id}_metadata.json")
        with open(meta_path, "w") as f:
            json.dump(study, f, indent=2)
        print(f"  Saved metadata: {os.path.basename(meta_path)}")

        # Try to find and download data files from the files section
        files_url = f"https://www.ebi.ac.uk/biostudies/files/{ae_id}/"
        # Download the processed data file(s) via ArrayExpress files API
        file_list_url = f"https://www.ebi.ac.uk/arrayexpress/files/{ae_id}/"
        
        # Try the newer BioStudies file listing
        ftp_base = f"https://ftp.ebi.ac.uk/biostudies/fire/E-MEXP-/{ae_id.replace('E-MEXP-', '')}/{ae_id}/"
        
        # Attempt direct download of common file patterns
        common_files = [
            f"{ae_id}.processed.1.zip",
            f"{ae_id}-processed-data-matrix.txt",
            f"{ae_id}.sdrf.txt",
            f"{ae_id}.idf.txt",
        ]

        any_success = False
        for fname in common_files:
            file_url = f"https://ftp.ebi.ac.uk/biostudies/nfs/E-MEXP-/{ae_id.replace('E-MEXP-', '')}/{ae_id}/Files/{fname}"
            file_path = os.path.join(dataset_dir, fname)
            if download_file(file_url, file_path, fname):
                any_success = True
            time.sleep(0.5)

        return any_success

    except Exception as e:
        print(f"  ERROR fetching ArrayExpress metadata: {e}")
        return False


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}\n")
    print(f"Downloading {len(DATASETS)} neurodegenerative disease datasets...\n")
    print("=" * 70)

    results = []

    for i, (dataset_id, description) in enumerate(DATASETS, 1):
        print(f"\n[{i}/{len(DATASETS)}] {dataset_id}")
        print(f"  {description}")
        print("-" * 50)

        dataset_dir = os.path.join(OUTPUT_DIR, dataset_id)

        if dataset_id.startswith("E-GEOD-"):
            # Convert to GEO ID: E-GEOD-21942 -> GSE21942
            geo_num = dataset_id.replace("E-GEOD-", "")
            geo_id = f"GSE{geo_num}"
            success = download_geo_dataset(geo_id, dataset_dir)
        elif dataset_id.startswith("E-MEXP-"):
            success = download_arrayexpress_dataset(dataset_id, dataset_dir)
        else:
            print(f"  Unknown dataset type: {dataset_id}")
            success = False

        results.append((dataset_id, description, success))
        time.sleep(1)  # Be polite to servers

    # Summary
    print("\n" + "=" * 70)
    print("DOWNLOAD SUMMARY")
    print("=" * 70)
    succeeded = sum(1 for _, _, s in results if s)
    failed = sum(1 for _, _, s in results if not s)

    for dataset_id, description, success in results:
        status = "OK" if success else "FAILED"
        print(f"  [{status}] {dataset_id} - {description}")

    print(f"\nTotal: {succeeded} succeeded, {failed} failed out of {len(results)}")
    print(f"Files saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
