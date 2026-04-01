#!/usr/bin/env python3
"""Debug platform annotation loading."""
import os

BASE = "/Users/ambuj.mehta/code/Git_Repos/Bio_Claude_testing/platforms"
annot = os.path.join(BASE, "GPL570.annot")

with open(annot, "r", errors="replace") as f:
    header_found = False
    gene_col = None
    count = 0
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if not header_found:
            for i, col in enumerate(parts):
                if col.strip().lower() == "gene symbol":
                    gene_col = i
                    break
            header_found = True
            print(f"Gene symbol column index: {gene_col}")
            continue
        if gene_col is not None and len(parts) > gene_col:
            sym = parts[gene_col].strip()
            targets = ["LATS1", "LATS2", "YAP1", "MTOR", "AKT1", "PTEN", "BDNF", "TP53", "SIRT1"]
            if sym and sym in targets:
                print(f"  Probe: {parts[0].strip()} -> Gene: {sym}")
                count += 1
            if count > 20:
                break

# Now test the load function
import sys
sys.path.insert(0, "/Users/ambuj.mehta/code/Git_Repos/Bio_Claude_testing")
from kibra_analysis import load_platform_annotation

p = load_platform_annotation("GPL570")
print(f"\nTotal probes with annotations: {len(p)}")
print(f"Sample entries:")
for k, v in list(p.items())[:5]:
    print(f"  {k}: {v} (type={type(v)})")

# Search for pathway genes
pathway_genes = ["LATS1", "LATS2", "YAP1", "WWTR1", "STK3", "STK4", "AKT1", "MTOR", "PTEN", "BDNF", "TP53", "SIRT1"]
found = {}
for probe, genes in p.items():
    for g in genes:
        if g in pathway_genes:
            found[g] = found.get(g, []) + [probe]

print(f"\nPathway genes found in GPL570:")
for g, probes in sorted(found.items()):
    print(f"  {g}: {probes}")

if not found:
    print("  NONE FOUND - checking why...")
    # Check what the actual gene symbols look like
    all_genes = set()
    for genes in p.values():
        all_genes.update(genes)
    # Check partial matches
    for target in pathway_genes:
        matches = [g for g in all_genes if target.lower() in g.lower()]
        if matches:
            print(f"  Partial match for {target}: {matches[:5]}")
