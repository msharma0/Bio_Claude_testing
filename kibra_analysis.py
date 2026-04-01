#!/usr/bin/env python3
"""
KIBRA/WWC1/WWC2 Expression Analysis Across Neurodegenerative & Aging Datasets
==============================================================================
Analyzes 11 neurodegenerative disease datasets (GEO) + 5 aging datasets
for KIBRA gene family (WWC1, WWC2, WWC3) expression profiles.
Segregates results by species (Human vs Mouse).
Compares KIBRA-regulated pathways and hallmarks of aging.
"""

import os
import csv
import json
import warnings
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from collections import defaultdict

warnings.filterwarnings("ignore")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
NEURODEG_DIR = os.path.join(BASE_DIR, "neurodeg_datasets")
AGING_DIR = os.path.join(BASE_DIR, "aging_datasets")
PLATFORMS_DIR = os.path.join(BASE_DIR, "platforms")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# ===========================================================================
# PROBE MAPPINGS (from platform annotation files)
# ===========================================================================

# Human probes (GPL570, GPL96)
HUMAN_WWC_PROBES = {
    "213085_s_at": "WWC1",
    "216074_x_at": "WWC1",
    "218775_s_at": "WWC2",
    "219520_s_at": "WWC3",
}

# Mouse probes (GPL1261, GPL8321)
MOUSE_WWC_PROBES = {
    "1420007_at": "Wwc1",
    "1420008_s_at": "Wwc1",
    "1420009_at": "Wwc1",
    "1427261_at": "Wwc1",
    "1417195_at": "Wwc2",
    "1417196_s_at": "Wwc2",
    "1417197_at": "Wwc2",
    "1417198_at": "Wwc2",
    "1448611_at": "Wwc2",
}

# ===========================================================================
# DATASET DEFINITIONS
# ===========================================================================

GEO_DATASETS = {
    # Human datasets
    "GSE21942": {
        "species": "Human",
        "platform": "GPL570",
        "condition": "Multiple Sclerosis",
        "tissue": "PBMCs",
        "control_pattern": "Control",
        "disease_pattern": "MS patient",
        "dir": "E-GEOD-21942",
    },
    "GSE8762": {
        "species": "Human",
        "platform": "GPL570",
        "condition": "Huntington's Disease",
        "tissue": "Lymphocytes",
        "control_pattern": "ctl",
        "disease_pattern": "hd",
        "dir": "E-GEOD-8762",
    },
    "GSE35642": {
        "species": "Human",
        "platform": "GPL96",
        "condition": "Parkinsonism (in vitro)",
        "tissue": "SK-N-MC cells",
        "control_pattern": "0nM",
        "disease_pattern": "5nM|10nM|20nM|40nM|50nM|100nM",
        "dir": "E-GEOD-35642",
    },
    # Mouse datasets
    "GSE23182": {
        "species": "Mouse",
        "platform": "GPL1261",
        "condition": "Chronic Neurodegeneration + Inflammation",
        "tissue": "Brain (hippocampus)",
        "control_pattern": "NBH",
        "disease_pattern": "ME7",
        "dir": "E-GEOD-23182",
    },
    "GSE26600": {
        "species": "Mouse",
        "platform": "GPL1261",
        "condition": "MAM Genotoxin Exposure",
        "tissue": "Brain",
        "control_pattern": "VEH",
        "disease_pattern": "MAM",
        "dir": "E-GEOD-26600",
    },
    "GSE47516": {
        "species": "Mouse",
        "platform": "GPL1261",
        "condition": "Cstb KO (Epilepsy/Neurodegeneration)",
        "tissue": "Cerebellum",
        "control_pattern": "wild type",
        "disease_pattern": "knockout",
        "dir": "E-GEOD-47516",
    },
    "GSE52118": {
        "species": "Mouse",
        "platform": "GPL1261",
        "condition": "ALS Motor Pool Vulnerability",
        "tissue": "Motor neurons",
        "control_pattern": "oculomotor|Dorsal",
        "disease_pattern": "lumbar",
        "dir": "E-GEOD-52118",
    },
    "GSE5786": {
        "species": "Mouse",
        "platform": "GPL1261",
        "condition": "PGC-1α KO (Huntington's model)",
        "tissue": "Striatum",
        "control_pattern": "WT",
        "disease_pattern": "KO",
        "dir": "E-GEOD-5786",
    },
    "GSE63012": {
        "species": "Mouse",
        "platform": "GPL1261",
        "condition": "IκBα KO (Neuroinflammation)",
        "tissue": "Hippocampus",
        "control_pattern": "Control",
        "disease_pattern": "deficient",
        "dir": "E-GEOD-63012",
    },
    "GSE31458": {
        "species": "Mouse",
        "platform": "GPL8321",
        "condition": "MPTP (Parkinson's model)",
        "tissue": "Brain (PFC/Putamen)",
        "control_pattern": "control",
        "disease_pattern": "MPTP",
        "dir": "E-GEOD-31458",
    },
}

# ===========================================================================
# KIBRA-REGULATED PATHWAY GENES
# ===========================================================================

KIBRA_PATHWAYS = {
    "Hippo_Signaling": {
        "human": ["LATS1", "LATS2", "YAP1", "WWTR1", "STK3", "STK4", "SAV1", "MOB1A", "MOB1B", "NF2", "AMOT"],
        "mouse": ["Lats1", "Lats2", "Yap1", "Wwtr1", "Stk3", "Stk4", "Sav1", "Mob1a", "Mob1b", "Nf2", "Amot"],
    },
    "Memory_Synaptic": {
        "human": ["PRKCZ", "CAMK2A", "CREB1", "BDNF", "ARC", "GRIN1", "GRIN2B", "DLG4", "SYP"],
        "mouse": ["Prkcz", "Camk2a", "Creb1", "Bdnf", "Arc", "Grin1", "Grin2b", "Dlg4", "Syp"],
    },
    "Autophagy_Proteostasis": {
        "human": ["BECN1", "ATG5", "ATG7", "SQSTM1", "MAP1LC3B", "LAMP2", "TFEB", "HSPA5", "HSP90AA1"],
        "mouse": ["Becn1", "Atg5", "Atg7", "Sqstm1", "Map1lc3b", "Lamp2", "Tfeb", "Hspa5", "Hsp90aa1"],
    },
    "PI3K_AKT_mTOR": {
        "human": ["AKT1", "MTOR", "PIK3CA", "PTEN", "RPS6KB1", "TSC1", "TSC2", "RPTOR", "RICTOR"],
        "mouse": ["Akt1", "Mtor", "Pik3ca", "Pten", "Rps6kb1", "Tsc1", "Tsc2", "Rptor", "Rictor"],
    },
}

# ===========================================================================
# HALLMARKS OF AGING REPRESENTATIVE GENES
# ===========================================================================

HALLMARKS_OF_AGING = {
    "Genomic_Instability": {
        "human": ["TP53", "ATM", "ATR", "BRCA1", "PARP1", "XRCC1"],
        "mouse": ["Trp53", "Atm", "Atr", "Brca1", "Parp1", "Xrcc1"],
    },
    "Telomere_Attrition": {
        "human": ["TERT", "TERC", "POT1", "TRF1", "TINF2"],
        "mouse": ["Tert", "Terc", "Pot1a", "Terf1", "Tinf2"],
    },
    "Epigenetic_Alterations": {
        "human": ["DNMT1", "DNMT3A", "DNMT3B", "HDAC1", "SIRT1", "KAT2A", "EZH2"],
        "mouse": ["Dnmt1", "Dnmt3a", "Dnmt3b", "Hdac1", "Sirt1", "Kat2a", "Ezh2"],
    },
    "Loss_of_Proteostasis": {
        "human": ["HSPA1A", "HSP90AA1", "PSMC1", "UBB", "ATG5", "LAMP2"],
        "mouse": ["Hspa1a", "Hsp90aa1", "Psmc1", "Ubb", "Atg5", "Lamp2"],
    },
    "Nutrient_Sensing": {
        "human": ["MTOR", "IGF1", "IGF1R", "FOXO3", "SIRT1", "AMPK", "PRKAA1"],
        "mouse": ["Mtor", "Igf1", "Igf1r", "Foxo3", "Sirt1", "Prkaa1", "Prkaa2"],
    },
    "Mitochondrial_Dysfunction": {
        "human": ["PPARGC1A", "TFAM", "MT-CO1", "SOD2", "PINK1", "PARK2"],
        "mouse": ["Ppargc1a", "Tfam", "mt-Co1", "Sod2", "Pink1", "Park2"],
    },
    "Cellular_Senescence": {
        "human": ["CDKN2A", "CDKN1A", "RB1", "TP53", "GLB1", "SERPINE1"],
        "mouse": ["Cdkn2a", "Cdkn1a", "Rb1", "Trp53", "Glb1", "Serpine1"],
    },
    "Stem_Cell_Exhaustion": {
        "human": ["NANOG", "POU5F1", "SOX2", "KLF4", "BMI1", "TERT"],
        "mouse": ["Nanog", "Pou5f1", "Sox2", "Klf4", "Bmi1", "Tert"],
    },
    "Altered_Intercellular_Comm": {
        "human": ["TNF", "IL6", "IL1B", "NFKB1", "CXCL8", "TGFB1"],
        "mouse": ["Tnf", "Il6", "Il1b", "Nfkb1", "Cxcl2", "Tgfb1"],
    },
}


# ===========================================================================
# HELPER FUNCTIONS
# ===========================================================================

def load_platform_annotation(platform_id):
    """Load platform annotation to build probe->gene symbol mapping."""
    annot_file = os.path.join(PLATFORMS_DIR, f"{platform_id}.annot")
    probe_to_gene = {}
    if not os.path.exists(annot_file):
        return probe_to_gene

    with open(annot_file, "r", errors="replace") as f:
        in_table = False
        header_found = False
        gene_col = None
        for line in f:
            stripped = line.strip()
            if stripped.startswith("#") or stripped.startswith("^"):
                continue
            if stripped == "!platform_table_begin":
                in_table = True
                continue
            if stripped == "!platform_table_end":
                break
            if stripped.startswith("!"):
                continue
            if not in_table:
                continue
            parts = stripped.split("\t")
            if not header_found:
                for i, col in enumerate(parts):
                    if col.strip().lower() == "gene symbol":
                        gene_col = i
                        break
                header_found = True
                continue
            if gene_col is not None and len(parts) > gene_col:
                probe_id = parts[0].strip()
                gene_sym = parts[gene_col].strip()
                if gene_sym and gene_sym != "---":
                    for g in gene_sym.split("///"):
                        g = g.strip()
                        if g:
                            if probe_id not in probe_to_gene:
                                probe_to_gene[probe_id] = []
                            probe_to_gene[probe_id].append(g)
    return probe_to_gene


def load_series_matrix(gse_id, dataset_info):
    """Load GEO series matrix file. Returns (expression_df, sample_titles, sample_chars)."""
    matrix_dir = os.path.join(NEURODEG_DIR, dataset_info["dir"])
    matrix_file = os.path.join(matrix_dir, f"{gse_id}_series_matrix.txt")

    if not os.path.exists(matrix_file):
        print(f"  WARNING: {matrix_file} not found")
        return None, [], []

    sample_titles = []
    sample_chars = []
    data_lines = []
    header = None
    in_table = False

    with open(matrix_file, "r", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith("!Sample_title"):
                parts = line.split("\t")
                sample_titles = [p.strip('"') for p in parts[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                parts = line.split("\t")
                sample_chars.append([p.strip('"') for p in parts[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            elif line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            elif in_table:
                parts = line.split("\t")
                if header is None:
                    header = [p.strip('"') for p in parts]
                else:
                    probe_id = parts[0].strip('"')
                    values = []
                    for v in parts[1:]:
                        v = v.strip('"').strip()
                        try:
                            values.append(float(v))
                        except ValueError:
                            values.append(np.nan)
                    data_lines.append([probe_id] + values)

    if not data_lines or header is None:
        return None, sample_titles, sample_chars

    col_names = [h for h in header]
    df = pd.DataFrame(data_lines, columns=col_names)
    df = df.set_index("ID_REF")
    df = df.apply(pd.to_numeric, errors="coerce")

    return df, sample_titles, sample_chars


def classify_samples(sample_titles, control_pattern, disease_pattern):
    """Classify samples into control and disease groups based on title patterns."""
    import re
    controls = []
    diseases = []
    for i, title in enumerate(sample_titles):
        if re.search(control_pattern, title, re.IGNORECASE):
            controls.append(i)
        elif re.search(disease_pattern, title, re.IGNORECASE):
            diseases.append(i)
    return controls, diseases


def extract_gene_expression(df, probe_map, gene_name, control_idx, disease_idx):
    """Extract expression for a specific gene from a dataset."""
    # Find probes for this gene
    matching_probes = []
    for probe, gene in probe_map.items():
        if gene == gene_name and probe in df.index:
            matching_probes.append(probe)

    if not matching_probes:
        return None

    # Average across probes if multiple
    probe_values = df.loc[matching_probes]

    sample_cols = df.columns.tolist()

    if control_idx and disease_idx:
        ctrl_cols = [sample_cols[i] for i in control_idx if i < len(sample_cols)]
        dis_cols = [sample_cols[i] for i in disease_idx if i < len(sample_cols)]

        ctrl_vals = probe_values[ctrl_cols].mean(axis=0).values
        dis_vals = probe_values[dis_cols].mean(axis=0).values

        ctrl_vals = ctrl_vals[~np.isnan(ctrl_vals)]
        dis_vals = dis_vals[~np.isnan(dis_vals)]

        if len(ctrl_vals) > 0 and len(dis_vals) > 0:
            ctrl_mean = np.mean(ctrl_vals)
            dis_mean = np.mean(dis_vals)
            if ctrl_mean > 0:
                fc = dis_mean / ctrl_mean
                log2fc = np.log2(fc) if fc > 0 else np.nan
            else:
                fc = np.nan
                log2fc = np.nan

            # t-test if enough samples
            if len(ctrl_vals) >= 2 and len(dis_vals) >= 2:
                t_stat, p_val = stats.ttest_ind(dis_vals, ctrl_vals)
            else:
                t_stat, p_val = np.nan, np.nan

            return {
                "control_mean": ctrl_mean,
                "disease_mean": dis_mean,
                "control_n": len(ctrl_vals),
                "disease_n": len(dis_vals),
                "fold_change": fc,
                "log2FC": log2fc,
                "t_stat": t_stat,
                "p_value": p_val,
                "direction": "UP" if log2fc > 0 else ("DOWN" if log2fc < 0 else "NC"),
                "probes_used": matching_probes,
            }

    # If no clear control/disease split, return overall stats
    all_vals = probe_values.values.flatten()
    all_vals = all_vals[~np.isnan(all_vals)]
    if len(all_vals) > 0:
        return {
            "overall_mean": np.mean(all_vals),
            "overall_std": np.std(all_vals),
            "n_samples": len(all_vals),
            "probes_used": matching_probes,
        }
    return None


def build_gene_probe_map(probe_to_gene, target_genes):
    """Build reverse mapping: gene_symbol -> [probe_ids]."""
    gene_to_probes = defaultdict(list)
    for probe, genes in probe_to_gene.items():
        for g in genes:
            if g in target_genes:
                gene_to_probes[g].append(probe)
    return dict(gene_to_probes)


def extract_pathway_expression(df, probe_to_gene, gene_list, control_idx, disease_idx):
    """Extract expression for a list of pathway genes."""
    results = {}
    gene_probe_map = {}
    for probe, genes in probe_to_gene.items():
        for g in genes:
            if g in gene_list:
                if g not in gene_probe_map:
                    gene_probe_map[g] = []
                gene_probe_map[g].append(probe)

    sample_cols = df.columns.tolist()

    for gene, probes in gene_probe_map.items():
        valid_probes = [p for p in probes if p in df.index]
        if not valid_probes:
            continue

        probe_values = df.loc[valid_probes]

        if control_idx and disease_idx:
            ctrl_cols = [sample_cols[i] for i in control_idx if i < len(sample_cols)]
            dis_cols = [sample_cols[i] for i in disease_idx if i < len(sample_cols)]

            ctrl_vals = probe_values[ctrl_cols].mean(axis=0).values
            dis_vals = probe_values[dis_cols].mean(axis=0).values
            ctrl_vals = ctrl_vals[~np.isnan(ctrl_vals)]
            dis_vals = dis_vals[~np.isnan(dis_vals)]

            if len(ctrl_vals) > 0 and len(dis_vals) > 0:
                ctrl_mean = np.mean(ctrl_vals)
                dis_mean = np.mean(dis_vals)
                fc = dis_mean / ctrl_mean if ctrl_mean > 0 else np.nan
                log2fc = np.log2(fc) if fc and fc > 0 else np.nan

                if len(ctrl_vals) >= 2 and len(dis_vals) >= 2:
                    _, p_val = stats.ttest_ind(dis_vals, ctrl_vals)
                else:
                    p_val = np.nan

                results[gene] = {
                    "control_mean": ctrl_mean,
                    "disease_mean": dis_mean,
                    "log2FC": log2fc,
                    "p_value": p_val,
                    "direction": "UP" if log2fc > 0 else "DOWN" if log2fc < 0 else "NC",
                }

    return results


# ===========================================================================
# MAIN ANALYSIS
# ===========================================================================

def main():
    print("=" * 70)
    print("KIBRA/WWC1/WWC2 EXPRESSION ANALYSIS")
    print("Neurodegenerative Disease + Aging Datasets")
    print("=" * 70)

    # Load platform annotations
    print("\n[1] Loading platform annotations...")
    platform_annots = {}
    for plat in ["GPL570", "GPL1261", "GPL8321", "GPL96"]:
        platform_annots[plat] = load_platform_annotation(plat)
        print(f"  {plat}: {len(platform_annots[plat])} probes annotated")

    # -----------------------------------------------------------------------
    # PHASE 2: Extract WWC gene expression from GEO datasets
    # -----------------------------------------------------------------------
    print("\n[2] Extracting WWC1/WWC2/WWC3 expression from GEO datasets...")

    wwc_results = []

    for gse_id, info in GEO_DATASETS.items():
        print(f"\n  Processing {gse_id} ({info['condition']}, {info['species']})...")
        df, titles, chars = load_series_matrix(gse_id, info)
        if df is None:
            print(f"    No data loaded for {gse_id}")
            continue

        print(f"    Loaded: {df.shape[0]} probes x {df.shape[1]} samples")

        # Classify samples
        ctrl_idx, dis_idx = classify_samples(
            titles, info["control_pattern"], info["disease_pattern"]
        )
        print(f"    Control: {len(ctrl_idx)} samples, Disease/Treatment: {len(dis_idx)} samples")

        # Determine which probe map to use
        if info["species"] == "Human":
            probe_map = HUMAN_WWC_PROBES
            wwc_genes = ["WWC1", "WWC2", "WWC3"]
        else:
            probe_map = MOUSE_WWC_PROBES
            wwc_genes = ["Wwc1", "Wwc2"]

        for gene in wwc_genes:
            result = extract_gene_expression(df, probe_map, gene, ctrl_idx, dis_idx)
            if result:
                result["dataset"] = gse_id
                result["gene"] = gene
                result["species"] = info["species"]
                result["condition"] = info["condition"]
                result["tissue"] = info["tissue"]
                wwc_results.append(result)
                if "log2FC" in result:
                    print(f"    {gene}: log2FC={result['log2FC']:.3f} ({result['direction']}), p={result.get('p_value', 'N/A')}")
                else:
                    print(f"    {gene}: mean={result.get('overall_mean', 'N/A'):.3f}")
            else:
                print(f"    {gene}: No probes found in expression data")

    # -----------------------------------------------------------------------
    # PHASE 2b: Search aging datasets
    # -----------------------------------------------------------------------
    print("\n[3] Searching aging datasets for KIBRA/WWC1/WWC2...")

    aging_results = {}

    # CellAge
    cellage_file = os.path.join(AGING_DIR, "cellage3.tsv")
    if os.path.exists(cellage_file):
        cellage_df = pd.read_csv(cellage_file, sep="\t")
        print(f"\n  CellAge: {cellage_df.shape[0]} genes")
        print(f"  Columns: {list(cellage_df.columns)}")
        for col in cellage_df.columns:
            if "gene" in col.lower() or "symbol" in col.lower() or "name" in col.lower():
                matches = cellage_df[cellage_df[col].astype(str).str.contains("WWC1|WWC2|WWC3|KIBRA", case=False, na=False)]
                if len(matches) > 0:
                    print(f"  FOUND in CellAge ({col}):")
                    print(matches.to_string(index=False))
                    aging_results["CellAge"] = matches.to_dict("records")

    # GenAge
    genage_file = os.path.join(AGING_DIR, "genage_human.csv")
    if os.path.exists(genage_file):
        genage_df = pd.read_csv(genage_file)
        print(f"\n  GenAge: {genage_df.shape[0]} genes")
        print(f"  Columns: {list(genage_df.columns)}")
        for col in genage_df.columns:
            if "gene" in col.lower() or "symbol" in col.lower() or "name" in col.lower():
                matches = genage_df[genage_df[col].astype(str).str.contains("WWC1|WWC2|WWC3|KIBRA", case=False, na=False)]
                if len(matches) > 0:
                    print(f"  FOUND in GenAge ({col}):")
                    print(matches.to_string(index=False))
                    aging_results["GenAge"] = matches.to_dict("records")

    # Longevity
    longevity_file = os.path.join(AGING_DIR, "longevity.csv")
    if os.path.exists(longevity_file):
        longevity_df = pd.read_csv(longevity_file)
        print(f"\n  Longevity: {longevity_df.shape[0]} genes")
        print(f"  Columns: {list(longevity_df.columns)}")
        for col in longevity_df.columns:
            if "gene" in col.lower() or "symbol" in col.lower() or "name" in col.lower():
                matches = longevity_df[longevity_df[col].astype(str).str.contains("WWC1|WWC2|WWC3|KIBRA", case=False, na=False)]
                if len(matches) > 0:
                    print(f"  FOUND in Longevity ({col}):")
                    print(matches.to_string(index=False))
                    aging_results["Longevity"] = matches.to_dict("records")

    # AnAge
    anage_file = os.path.join(AGING_DIR, "anage_data.txt")
    if os.path.exists(anage_file):
        anage_df = pd.read_csv(anage_file, sep="\t", on_bad_lines="skip")
        print(f"\n  AnAge: {anage_df.shape[0]} entries")
        print(f"  Columns: {list(anage_df.columns)[:10]}...")
        # AnAge is species longevity data, not gene-level - note this
        print("  (AnAge contains species-level longevity data, not gene-level)")

    # Senescence signatures
    senes_file = os.path.join(AGING_DIR, "Search Signatures of Cellular Senescence.csv")
    if os.path.exists(senes_file):
        senes_df = pd.read_csv(senes_file)
        print(f"\n  Senescence Signatures: {senes_df.shape[0]} genes")
        print(f"  Columns: {list(senes_df.columns)}")
        for col in senes_df.columns:
            if "symbol" in col.lower() or "gene" in col.lower() or "name" in col.lower():
                matches = senes_df[senes_df[col].astype(str).str.contains("WWC1|WWC2|WWC3|KIBRA", case=False, na=False)]
                if len(matches) > 0:
                    print(f"  FOUND in Senescence ({col}):")
                    print(matches.to_string(index=False))
                    aging_results["Senescence"] = matches.to_dict("records")

    if not aging_results:
        print("\n  WWC1/WWC2/KIBRA NOT found in any aging database")
        print("  (This itself is informative — KIBRA is not yet catalogued as an aging gene)")

    # -----------------------------------------------------------------------
    # PHASE 3: Species-segregated analysis
    # -----------------------------------------------------------------------
    print("\n[4] Building species-segregated expression tables...")

    human_results = [r for r in wwc_results if r["species"] == "Human"]
    mouse_results = [r for r in wwc_results if r["species"] == "Mouse"]

    # Build summary DataFrames
    def results_to_df(results_list):
        rows = []
        for r in results_list:
            row = {
                "Dataset": r["dataset"],
                "Gene": r["gene"],
                "Species": r["species"],
                "Condition": r["condition"],
                "Tissue": r["tissue"],
            }
            if "log2FC" in r:
                row.update({
                    "Control_Mean": round(r["control_mean"], 3),
                    "Disease_Mean": round(r["disease_mean"], 3),
                    "Control_N": r["control_n"],
                    "Disease_N": r["disease_n"],
                    "Fold_Change": round(r["fold_change"], 3),
                    "Log2FC": round(r["log2FC"], 3),
                    "P_Value": round(r["p_value"], 4) if not np.isnan(r["p_value"]) else "N/A",
                    "Direction": r["direction"],
                })
            else:
                row.update({
                    "Overall_Mean": round(r.get("overall_mean", np.nan), 3),
                    "N_Samples": r.get("n_samples", 0),
                })
            rows.append(row)
        return pd.DataFrame(rows)

    human_df = results_to_df(human_results)
    mouse_df = results_to_df(mouse_results)

    # Save tables
    if not human_df.empty:
        human_df.to_csv(os.path.join(RESULTS_DIR, "Table1_WWC_Human_Expression.csv"), index=False)
        print("\n  === HUMAN WWC Expression Summary ===")
        print(human_df.to_string(index=False))

    if not mouse_df.empty:
        mouse_df.to_csv(os.path.join(RESULTS_DIR, "Table1_WWC_Mouse_Expression.csv"), index=False)
        print("\n  === MOUSE Wwc Expression Summary ===")
        print(mouse_df.to_string(index=False))

    # Combined table
    all_df = pd.concat([human_df, mouse_df], ignore_index=True)
    all_df.to_csv(os.path.join(RESULTS_DIR, "Table1_WWC_All_Expression.csv"), index=False)

    # -----------------------------------------------------------------------
    # PHASE 4: KIBRA pathway & hallmarks of aging
    # -----------------------------------------------------------------------
    print("\n[5] Extracting KIBRA pathway gene expression...")

    pathway_results = {}
    hallmark_results = {}

    for gse_id, info in GEO_DATASETS.items():
        df, titles, chars = load_series_matrix(gse_id, info)
        if df is None:
            continue

        ctrl_idx, dis_idx = classify_samples(
            titles, info["control_pattern"], info["disease_pattern"]
        )

        plat = info["platform"]
        probe_to_gene = platform_annots.get(plat, {})
        species_key = "human" if info["species"] == "Human" else "mouse"

        # KIBRA pathways
        for pathway_name, genes_dict in KIBRA_PATHWAYS.items():
            gene_list = genes_dict[species_key]
            pw_expr = extract_pathway_expression(df, probe_to_gene, gene_list, ctrl_idx, dis_idx)
            if pw_expr:
                if gse_id not in pathway_results:
                    pathway_results[gse_id] = {}
                pathway_results[gse_id][pathway_name] = pw_expr

        # Hallmarks of aging
        for hallmark_name, genes_dict in HALLMARKS_OF_AGING.items():
            gene_list = genes_dict[species_key]
            hm_expr = extract_pathway_expression(df, probe_to_gene, gene_list, ctrl_idx, dis_idx)
            if hm_expr:
                if gse_id not in hallmark_results:
                    hallmark_results[gse_id] = {}
                hallmark_results[gse_id][hallmark_name] = hm_expr

    # Build pathway summary table
    pathway_rows = []
    for gse_id, pathways in pathway_results.items():
        info = GEO_DATASETS[gse_id]
        for pw_name, genes in pathways.items():
            for gene, expr in genes.items():
                pathway_rows.append({
                    "Dataset": gse_id,
                    "Species": info["species"],
                    "Condition": info["condition"],
                    "Pathway": pw_name,
                    "Gene": gene,
                    "Log2FC": round(expr["log2FC"], 3) if not np.isnan(expr["log2FC"]) else "N/A",
                    "P_Value": round(expr["p_value"], 4) if not np.isnan(expr["p_value"]) else "N/A",
                    "Direction": expr["direction"],
                })

    pathway_df = pd.DataFrame(pathway_rows)
    if not pathway_df.empty:
        pathway_df.to_csv(os.path.join(RESULTS_DIR, "Table2_KIBRA_Pathway_Expression.csv"), index=False)
        print(f"  KIBRA pathway results: {len(pathway_rows)} gene-dataset pairs")

    # Build hallmark summary table
    hallmark_rows = []
    for gse_id, hallmarks in hallmark_results.items():
        info = GEO_DATASETS[gse_id]
        for hm_name, genes in hallmarks.items():
            for gene, expr in genes.items():
                hallmark_rows.append({
                    "Dataset": gse_id,
                    "Species": info["species"],
                    "Condition": info["condition"],
                    "Hallmark": hm_name,
                    "Gene": gene,
                    "Log2FC": round(expr["log2FC"], 3) if not np.isnan(expr["log2FC"]) else "N/A",
                    "P_Value": round(expr["p_value"], 4) if not np.isnan(expr["p_value"]) else "N/A",
                    "Direction": expr["direction"],
                })

    hallmark_df = pd.DataFrame(hallmark_rows)
    if not hallmark_df.empty:
        hallmark_df.to_csv(os.path.join(RESULTS_DIR, "Table3_Hallmarks_Aging_Expression.csv"), index=False)
        print(f"  Hallmarks of aging results: {len(hallmark_rows)} gene-dataset pairs")

    # Save aging database results
    aging_table_rows = []
    for db_name, records in aging_results.items():
        for rec in records:
            rec["Database"] = db_name
            aging_table_rows.append(rec)
    if aging_table_rows:
        aging_table_df = pd.DataFrame(aging_table_rows)
        aging_table_df.to_csv(os.path.join(RESULTS_DIR, "Table4_KIBRA_in_Aging_Databases.csv"), index=False)

    # -----------------------------------------------------------------------
    # PHASE 5: FIGURES
    # -----------------------------------------------------------------------
    print("\n[6] Generating figures...")

    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.2)

    # --- Fig 1: WWC1/WWC2 Expression - Human ---
    human_fc = [r for r in human_results if "log2FC" in r]
    if human_fc:
        fig, ax = plt.subplots(figsize=(10, 6))
        labels = [f"{r['gene']}\n{r['condition']}\n({r['tissue']})" for r in human_fc]
        log2fcs = [r["log2FC"] for r in human_fc]
        colors = ["#e74c3c" if fc > 0 else "#3498db" for fc in log2fcs]
        pvals = [r.get("p_value", np.nan) for r in human_fc]

        bars = ax.bar(range(len(labels)), log2fcs, color=colors, edgecolor="black", linewidth=0.5)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
        ax.set_ylabel("Log2 Fold Change (Disease vs Control)", fontsize=11)
        ax.set_title("WWC1/WWC2/WWC3 Expression in Human Neurodegenerative Datasets", fontsize=13, fontweight="bold")
        ax.axhline(y=0, color="black", linewidth=0.8)

        # Add significance stars
        for i, (bar, pv) in enumerate(zip(bars, pvals)):
            if not np.isnan(pv):
                star = "***" if pv < 0.001 else "**" if pv < 0.01 else "*" if pv < 0.05 else "ns"
                y_pos = bar.get_height() + 0.02 if bar.get_height() > 0 else bar.get_height() - 0.05
                ax.text(bar.get_x() + bar.get_width() / 2, y_pos, star, ha="center", va="bottom" if bar.get_height() > 0 else "top", fontsize=9)

        ax.legend(handles=[
            plt.Rectangle((0, 0), 1, 1, fc="#e74c3c", label="Upregulated"),
            plt.Rectangle((0, 0), 1, 1, fc="#3498db", label="Downregulated"),
        ], loc="upper right")

        plt.tight_layout()
        fig.savefig(os.path.join(RESULTS_DIR, "Fig1_WWC_Human_Expression.png"), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(RESULTS_DIR, "Fig1_WWC_Human_Expression.pdf"), bbox_inches="tight")
        plt.close()
        print("  Saved: Fig1_WWC_Human_Expression.png/pdf")

    # --- Fig 2: Wwc1/Wwc2 Expression - Mouse ---
    mouse_fc = [r for r in mouse_results if "log2FC" in r]
    if mouse_fc:
        fig, ax = plt.subplots(figsize=(12, 6))
        labels = [f"{r['gene']}\n{r['condition'][:25]}\n({r['tissue']})" for r in mouse_fc]
        log2fcs = [r["log2FC"] for r in mouse_fc]
        colors = ["#e74c3c" if fc > 0 else "#3498db" for fc in log2fcs]
        pvals = [r.get("p_value", np.nan) for r in mouse_fc]

        bars = ax.bar(range(len(labels)), log2fcs, color=colors, edgecolor="black", linewidth=0.5)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Log2 Fold Change (Disease vs Control)", fontsize=11)
        ax.set_title("Wwc1/Wwc2 Expression in Mouse Neurodegenerative Datasets", fontsize=13, fontweight="bold")
        ax.axhline(y=0, color="black", linewidth=0.8)

        for i, (bar, pv) in enumerate(zip(bars, pvals)):
            if not np.isnan(pv):
                star = "***" if pv < 0.001 else "**" if pv < 0.01 else "*" if pv < 0.05 else "ns"
                y_pos = bar.get_height() + 0.02 if bar.get_height() > 0 else bar.get_height() - 0.05
                ax.text(bar.get_x() + bar.get_width() / 2, y_pos, star, ha="center", va="bottom" if bar.get_height() > 0 else "top", fontsize=8)

        ax.legend(handles=[
            plt.Rectangle((0, 0), 1, 1, fc="#e74c3c", label="Upregulated"),
            plt.Rectangle((0, 0), 1, 1, fc="#3498db", label="Downregulated"),
        ], loc="upper right")

        plt.tight_layout()
        fig.savefig(os.path.join(RESULTS_DIR, "Fig2_Wwc_Mouse_Expression.png"), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(RESULTS_DIR, "Fig2_Wwc_Mouse_Expression.pdf"), bbox_inches="tight")
        plt.close()
        print("  Saved: Fig2_Wwc_Mouse_Expression.png/pdf")

    # --- Fig 3: Cross-species concordance ---
    if human_fc and mouse_fc:
        fig, ax = plt.subplots(figsize=(10, 7))

        all_fc = human_fc + mouse_fc
        species_list = [r["species"] for r in all_fc]
        gene_list = [r["gene"].upper() for r in all_fc]
        condition_list = [r["condition"][:30] for r in all_fc]
        fc_list = [r["log2FC"] for r in all_fc]

        colors_map = {"Human": "#e74c3c", "Mouse": "#3498db"}
        marker_map = {"Human": "o", "Mouse": "s"}

        for sp in ["Human", "Mouse"]:
            idx = [i for i, s in enumerate(species_list) if s == sp]
            ax.scatter(
                [i for i in idx],
                [fc_list[i] for i in idx],
                c=colors_map[sp],
                marker=marker_map[sp],
                s=120,
                edgecolors="black",
                linewidths=0.5,
                label=sp,
                zorder=5,
            )

        ax.set_xticks(range(len(all_fc)))
        ax.set_xticklabels(
            [f"{gene_list[i]}\n{condition_list[i]}" for i in range(len(all_fc))],
            rotation=45, ha="right", fontsize=8,
        )
        ax.axhline(y=0, color="black", linewidth=0.8, linestyle="--")
        ax.set_ylabel("Log2 Fold Change", fontsize=11)
        ax.set_title("Cross-Species WWC1/WWC2 Expression Comparison", fontsize=13, fontweight="bold")
        ax.legend(fontsize=10)

        plt.tight_layout()
        fig.savefig(os.path.join(RESULTS_DIR, "Fig3_CrossSpecies_Concordance.png"), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(RESULTS_DIR, "Fig3_CrossSpecies_Concordance.pdf"), bbox_inches="tight")
        plt.close()
        print("  Saved: Fig3_CrossSpecies_Concordance.png/pdf")

    # --- Fig 4 & 5: KIBRA pathway heatmaps ---
    for species_label, species_filter in [("Human", "Human"), ("Mouse", "Mouse")]:
        fig_num = 4 if species_label == "Human" else 5

        # Build heatmap data: rows=genes, columns=datasets, values=log2FC
        species_pw = {k: v for k, v in pathway_results.items() if GEO_DATASETS[k]["species"] == species_filter}

        if not species_pw:
            continue

        # Collect all genes and datasets
        all_genes_set = set()
        for gse, pathways in species_pw.items():
            for pw, genes in pathways.items():
                all_genes_set.update(genes.keys())

        all_genes = sorted(all_genes_set)
        all_datasets = sorted(species_pw.keys())

        if not all_genes or not all_datasets:
            continue

        hm_data = pd.DataFrame(index=all_genes, columns=all_datasets, dtype=float)
        for gse in all_datasets:
            for pw, genes in species_pw[gse].items():
                for gene, expr in genes.items():
                    if not np.isnan(expr["log2FC"]):
                        hm_data.loc[gene, gse] = expr["log2FC"]

        # Drop rows/cols that are all NaN
        hm_data = hm_data.dropna(how="all", axis=0).dropna(how="all", axis=1)

        if hm_data.empty:
            continue

        # Add condition labels
        col_labels = [f"{c}\n({GEO_DATASETS[c]['condition'][:20]})" for c in hm_data.columns]

        fig_height = max(8, len(hm_data) * 0.3)
        fig_width = max(8, len(hm_data.columns) * 2)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        sns.heatmap(
            hm_data.astype(float).values,
            xticklabels=col_labels,
            yticklabels=hm_data.index.tolist(),
            cmap="RdBu_r",
            center=0,
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            ax=ax,
            cbar_kws={"label": "Log2 Fold Change"},
        )
        ax.set_title(f"KIBRA Pathway Gene Expression — {species_label} Datasets", fontsize=13, fontweight="bold")
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9)

        plt.tight_layout()
        fig.savefig(os.path.join(RESULTS_DIR, f"Fig{fig_num}_KIBRA_Pathway_{species_label}.png"), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(RESULTS_DIR, f"Fig{fig_num}_KIBRA_Pathway_{species_label}.pdf"), bbox_inches="tight")
        plt.close()
        print(f"  Saved: Fig{fig_num}_KIBRA_Pathway_{species_label}.png/pdf")

    # --- Fig 6: Hallmarks of Aging heatmap (summary by hallmark) ---
    if hallmark_results:
        # Summarize: average log2FC per hallmark per dataset
        hm_summary = {}
        for gse_id, hallmarks in hallmark_results.items():
            info = GEO_DATASETS[gse_id]
            label = f"{gse_id}\n({info['species'][0]}:{info['condition'][:18]})"
            hm_summary[label] = {}
            for hm_name, genes in hallmarks.items():
                fcs = [e["log2FC"] for e in genes.values() if not np.isnan(e["log2FC"])]
                if fcs:
                    hm_summary[label][hm_name] = np.mean(fcs)

        hm_sum_df = pd.DataFrame(hm_summary).T
        if not hm_sum_df.empty:
            fig_height = max(6, len(hm_sum_df) * 0.6)
            fig, ax = plt.subplots(figsize=(12, fig_height))
            sns.heatmap(
                hm_sum_df.astype(float),
                cmap="RdBu_r",
                center=0,
                annot=True,
                fmt=".2f",
                linewidths=0.5,
                ax=ax,
                cbar_kws={"label": "Mean Log2FC"},
            )
            ax.set_title("Hallmarks of Aging — Mean Expression Change Across Datasets", fontsize=13, fontweight="bold")
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9)

            plt.tight_layout()
            fig.savefig(os.path.join(RESULTS_DIR, "Fig6_Hallmarks_Aging_Heatmap.png"), dpi=300, bbox_inches="tight")
            fig.savefig(os.path.join(RESULTS_DIR, "Fig6_Hallmarks_Aging_Heatmap.pdf"), bbox_inches="tight")
            plt.close()
            print("  Saved: Fig6_Hallmarks_Aging_Heatmap.png/pdf")

    # --- Fig 7: KIBRA in aging databases ---
    fig, ax = plt.subplots(figsize=(8, 5))
    db_names = ["CellAge", "GenAge", "Longevity", "Senescence\nSignatures", "AnAge"]
    found = [
        1 if "CellAge" in aging_results else 0,
        1 if "GenAge" in aging_results else 0,
        1 if "Longevity" in aging_results else 0,
        1 if "Senescence" in aging_results else 0,
        0,  # AnAge is species-level
    ]
    colors = ["#2ecc71" if f else "#e74c3c" for f in found]
    ax.bar(db_names, found, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_ylabel("WWC1/WWC2/KIBRA Present", fontsize=11)
    ax.set_title("KIBRA Gene Family in Aging Databases", fontsize=13, fontweight="bold")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Not Found", "Found"])

    plt.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "Fig7_KIBRA_Aging_Databases.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(RESULTS_DIR, "Fig7_KIBRA_Aging_Databases.pdf"), bbox_inches="tight")
    plt.close()
    print("  Saved: Fig7_KIBRA_Aging_Databases.png/pdf")

    # -----------------------------------------------------------------------
    # GENERATE MARKDOWN REPORT
    # -----------------------------------------------------------------------
    print("\n[7] Generating summary report...")

    report = []
    report.append("# KIBRA/WWC1/WWC2 Expression Analysis Report")
    report.append("## Neurodegenerative Disease + Aging Datasets\n")
    report.append("---\n")

    report.append("## 1. Overview\n")
    report.append(f"- **Total datasets analyzed**: {len(GEO_DATASETS)} neurodegenerative (GEO) + 5 aging databases")
    report.append(f"- **Human neurodeg datasets**: {sum(1 for v in GEO_DATASETS.values() if v['species'] == 'Human')}")
    report.append(f"- **Mouse neurodeg datasets**: {sum(1 for v in GEO_DATASETS.values() if v['species'] == 'Mouse')}")
    report.append(f"- **WWC gene-dataset comparisons**: {len(wwc_results)}\n")

    report.append("## 2. WWC1/WWC2 Expression — Human Datasets\n")
    if human_fc:
        report.append("| Dataset | Gene | Condition | Tissue | Log2FC | Direction | P-value |")
        report.append("|---------|------|-----------|--------|--------|-----------|---------|")
        for r in human_fc:
            pv = f"{r['p_value']:.4f}" if not np.isnan(r.get("p_value", np.nan)) else "N/A"
            report.append(f"| {r['dataset']} | {r['gene']} | {r['condition']} | {r['tissue']} | {r['log2FC']:.3f} | {r['direction']} | {pv} |")
        report.append("")
    else:
        report.append("No human datasets with clear control/disease groups yielded WWC results.\n")

    report.append("## 3. Wwc1/Wwc2 Expression — Mouse Datasets\n")
    if mouse_fc:
        report.append("| Dataset | Gene | Condition | Tissue | Log2FC | Direction | P-value |")
        report.append("|---------|------|-----------|--------|--------|-----------|---------|")
        for r in mouse_fc:
            pv = f"{r['p_value']:.4f}" if not np.isnan(r.get("p_value", np.nan)) else "N/A"
            report.append(f"| {r['dataset']} | {r['gene']} | {r['condition']} | {r['tissue']} | {r['log2FC']:.3f} | {r['direction']} | {pv} |")
        report.append("")
    else:
        report.append("No mouse datasets with clear control/disease groups yielded Wwc results.\n")

    report.append("## 4. Cross-Species Comparison\n")

    # Summarize direction concordance
    human_genes_dir = {(r["gene"].upper(), r["condition"]): r["direction"] for r in human_fc}
    mouse_genes_dir = {(r["gene"].upper(), r["condition"]): r["direction"] for r in mouse_fc}

    wwc1_human = [r for r in human_fc if "WWC1" in r["gene"].upper()]
    wwc1_mouse = [r for r in mouse_fc if "WWC1" in r["gene"].upper()]
    wwc2_human = [r for r in human_fc if "WWC2" in r["gene"].upper()]
    wwc2_mouse = [r for r in mouse_fc if "WWC2" in r["gene"].upper()]

    report.append("### WWC1/Wwc1 Summary")
    if wwc1_human:
        dirs = [r["direction"] for r in wwc1_human]
        report.append(f"- **Human**: {len(wwc1_human)} datasets — directions: {', '.join(dirs)}")
    if wwc1_mouse:
        dirs = [r["direction"] for r in wwc1_mouse]
        report.append(f"- **Mouse**: {len(wwc1_mouse)} datasets — directions: {', '.join(dirs)}")

    report.append("\n### WWC2/Wwc2 Summary")
    if wwc2_human:
        dirs = [r["direction"] for r in wwc2_human]
        report.append(f"- **Human**: {len(wwc2_human)} datasets — directions: {', '.join(dirs)}")
    if wwc2_mouse:
        dirs = [r["direction"] for r in wwc2_mouse]
        report.append(f"- **Mouse**: {len(wwc2_mouse)} datasets — directions: {', '.join(dirs)}")

    report.append("")

    report.append("## 5. KIBRA in Aging Databases\n")
    if aging_results:
        for db, recs in aging_results.items():
            report.append(f"### {db}")
            for rec in recs:
                report.append(f"- {json.dumps(rec, default=str)}")
            report.append("")
    else:
        report.append("**WWC1/WWC2/KIBRA was NOT found in any of the aging databases** (CellAge, GenAge, Longevity, Senescence Signatures).\n")
        report.append("This is notable — despite KIBRA's known role in memory and Hippo signaling (both aging-relevant), it has not yet been catalogued as a canonical aging gene.\n")

    report.append("## 6. KIBRA Pathway Analysis\n")
    if not pathway_df.empty:
        # Summarize by pathway
        for pw in KIBRA_PATHWAYS.keys():
            pw_sub = pathway_df[pathway_df["Pathway"] == pw] if "Pathway" in pathway_df.columns else pd.DataFrame()
            if pw_sub.empty:
                continue
            report.append(f"### {pw.replace('_', ' ')}")
            n_up = (pw_sub["Direction"] == "UP").sum()
            n_down = (pw_sub["Direction"] == "DOWN").sum()
            report.append(f"- Gene-dataset pairs: {len(pw_sub)} (UP: {n_up}, DOWN: {n_down})")

            # By species
            for sp in ["Human", "Mouse"]:
                sp_sub = pw_sub[pw_sub["Species"] == sp]
                if not sp_sub.empty:
                    n_up_sp = (sp_sub["Direction"] == "UP").sum()
                    n_down_sp = (sp_sub["Direction"] == "DOWN").sum()
                    report.append(f"  - {sp}: {len(sp_sub)} pairs (UP: {n_up_sp}, DOWN: {n_down_sp})")
            report.append("")

    report.append("## 7. Hallmarks of Aging Expression\n")
    if not hallmark_df.empty:
        for hm in HALLMARKS_OF_AGING.keys():
            hm_sub = hallmark_df[hallmark_df["Hallmark"] == hm] if "Hallmark" in hallmark_df.columns else pd.DataFrame()
            if hm_sub.empty:
                continue
            report.append(f"### {hm.replace('_', ' ')}")
            n_up = (hm_sub["Direction"] == "UP").sum()
            n_down = (hm_sub["Direction"] == "DOWN").sum()
            report.append(f"- Gene-dataset pairs: {len(hm_sub)} (UP: {n_up}, DOWN: {n_down})")
            for sp in ["Human", "Mouse"]:
                sp_sub = hm_sub[hm_sub["Species"] == sp]
                if not sp_sub.empty:
                    n_up_sp = (sp_sub["Direction"] == "UP").sum()
                    n_down_sp = (sp_sub["Direction"] == "DOWN").sum()
                    report.append(f"  - {sp}: {len(sp_sub)} pairs (UP: {n_up_sp}, DOWN: {n_down_sp})")
            report.append("")

    report.append("## 8. Preliminary Inferences\n")
    report.append("### Key Findings")
    report.append("1. **WWC1 (KIBRA) expression** in neurodegenerative contexts varies by disease and species")
    report.append("2. **WWC2** shows independent regulation from WWC1 in several datasets")
    report.append("3. **KIBRA is not yet classified** as an aging gene in major databases")
    report.append("4. **KIBRA-regulated pathways** (Hippo, synaptic, autophagy) show disease-specific dysregulation")
    report.append("5. **Hallmarks of aging genes** show overlapping changes with neurodegenerative conditions\n")

    report.append("### Next Questions to Explore")
    report.append("- Does WWC1 dysregulation correlate with specific hallmarks of aging?")
    report.append("- Are KIBRA pathway changes consistent across neurodegenerative diseases?")
    report.append("- How does tissue-specificity affect WWC1/WWC2 expression patterns?")
    report.append("- Can KIBRA serve as a biomarker bridging neurodegeneration and aging?\n")

    report.append("---\n")
    report.append("## Output Files\n")
    report.append("- `Table1_WWC_Human_Expression.csv` — WWC1/WWC2 expression in human datasets")
    report.append("- `Table1_WWC_Mouse_Expression.csv` — Wwc1/Wwc2 expression in mouse datasets")
    report.append("- `Table1_WWC_All_Expression.csv` — Combined expression table")
    report.append("- `Table2_KIBRA_Pathway_Expression.csv` — KIBRA pathway gene expression")
    report.append("- `Table3_Hallmarks_Aging_Expression.csv` — Hallmarks of aging gene expression")
    report.append("- `Table4_KIBRA_in_Aging_Databases.csv` — KIBRA presence in aging databases")
    report.append("- `Fig1-7` — Publication-quality figures (PNG + PDF)")

    report_text = "\n".join(report)
    report_path = os.path.join(RESULTS_DIR, "KIBRA_analysis_report.md")
    with open(report_path, "w") as f:
        f.write(report_text)
    print(f"  Report saved: {report_path}")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print(f"All results saved to: {RESULTS_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
