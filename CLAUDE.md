# KIBRA/WWC Gene Family — Neurodegenerative & Aging Analysis

## Project Goal

Explore the expression profiles of KIBRA gene homologs (WWC1/WWC2/WWC3) across neurodegenerative disease and aging datasets, segregated by species (Human vs Mouse). Compare KIBRA-regulated pathway activity and hallmarks of aging expression across conditions to assess whether KIBRA could bridge neurodegeneration and aging biology.

---

## Workspace Layout

```
Bio_Claude_testing/
├── CLAUDE.md                          # This file — project context
├── kibra_analysis.py                  # Main analysis script (all phases)
├── improve_fig6.py                    # Improved hallmarks of aging heatmap
├── download_neurodeg_datasets.py      # Download 11 neurodeg GEO datasets
├── fix_emexp.py                       # Fix E-MEXP-2280 download
├── download_platforms.py              # Download GPL platform annotations
├── fix_platforms.py                   # Fix GPL570/GPL96 annotation download
├── neurodeg_datasets/                 # 10 GEO series matrix datasets
│   ├── E-GEOD-21942/                  # Human — Multiple Sclerosis (GPL570)
│   ├── E-GEOD-8762/                   # Human — Huntington's Disease (GPL570)
│   ├── E-GEOD-35642/                  # Human — Parkinsonism in vitro (GPL96)
│   ├── E-GEOD-23182/                  # Mouse — Chronic Neurodegeneration (GPL1261)
│   ├── E-GEOD-26600/                  # Mouse — MAM Genotoxin (GPL1261)
│   ├── E-GEOD-47516/                  # Mouse — Cstb KO / Epilepsy (GPL1261)
│   ├── E-GEOD-52118/                  # Mouse — ALS Motor Pools (GPL1261)
│   ├── E-GEOD-5786/                   # Mouse — PGC-1α KO / Huntington's (GPL1261)
│   ├── E-GEOD-63012/                  # Mouse — IκBα KO / Neuroinflammation (GPL1261)
│   ├── E-GEOD-31458/                  # Mouse — MPTP Parkinson's (GPL8321)
│   └── E-MEXP-2280/                   # ArrayExpress — raw CEL files only (no matrix)
├── aging_datasets/                    # 5 aging databases
│   ├── cellage3.tsv                   # CellAge — cellular senescence genes
│   ├── genage_human.csv               # GenAge — human aging genes
│   ├── longevity.csv                  # LongevityMap — longevity variants
│   ├── anage_data.txt                 # AnAge — species longevity (not gene-level)
│   └── Search Signatures of Cellular Senescence.csv  # SenMayo/senescence signatures
├── platforms/                         # GEO platform annotations (.annot files)
│   ├── GPL570.annot                   # Human — Affymetrix HG-U133 Plus 2.0
│   ├── GPL96.annot                    # Human — Affymetrix HG-U133A
│   ├── GPL1261.annot                  # Mouse — Affymetrix Mouse 430 2.0
│   └── GPL8321.annot                  # Mouse — Affymetrix Mouse Gene 1.0 ST
└── results/                           # All output
    ├── KIBRA_analysis_report.md       # Full results report (tables + inferences)
    ├── Table1_WWC_Human_Expression.csv
    ├── Table1_WWC_Mouse_Expression.csv
    ├── Table1_WWC_All_Expression.csv
    ├── Table2_KIBRA_Pathway_Expression.csv
    ├── Table3_Hallmarks_Aging_Expression.csv
    ├── Fig1_WWC_Human_Expression.png/pdf
    ├── Fig2_Wwc_Mouse_Expression.png/pdf
    ├── Fig3_CrossSpecies_Concordance.png/pdf
    ├── Fig4_KIBRA_Pathway_Human.png/pdf
    ├── Fig5_KIBRA_Pathway_Mouse.png/pdf
    ├── Fig6_Hallmarks_Aging_Heatmap.png/pdf
    └── Fig7_KIBRA_Aging_Databases.png/pdf
```

---

## Gene Targets

| Gene   | Species | Also Known As |
|--------|---------|---------------|
| WWC1   | Human   | KIBRA         |
| WWC2   | Human   |               |
| WWC3   | Human   |               |
| Wwc1   | Mouse   | Kibra         |
| Wwc2   | Mouse   |               |

---

## Probe-to-Gene Mappings

### Human Probes (GPL570 / GPL96)

| Probe ID      | Gene |
|---------------|------|
| 213085_s_at   | WWC1 |
| 216074_x_at   | WWC1 |
| 218775_s_at   | WWC2 (also maps to CLDN22) |
| 219520_s_at   | WWC3 |

### Mouse Probes (GPL1261 / GPL8321)

| Probe ID      | Gene |
|---------------|------|
| 1420007_at    | Wwc1 |
| 1420008_s_at  | Wwc1 |
| 1420009_at    | Wwc1 |
| 1427261_at    | Wwc1 |
| 1417195_at    | Wwc2 |
| 1417196_s_at  | Wwc2 |
| 1417197_at    | Wwc2 |
| 1417198_at    | Wwc2 |
| 1448611_at    | Wwc2 |

---

## Dataset Details

### Human Neurodegenerative Datasets (3)

| GEO ID   | Condition              | Tissue          | Platform | Control Pattern | Disease Pattern |
|----------|------------------------|-----------------|----------|-----------------|-----------------|
| GSE21942 | Multiple Sclerosis     | PBMCs           | GPL570   | "Control"       | "MS patient"    |
| GSE8762  | Huntington's Disease   | Lymphocytes     | GPL570   | "ctl"           | "hd"            |
| GSE35642 | Parkinsonism (in vitro)| SK-N-MC cells   | GPL96    | "0nM"           | "5nM\|10nM\|..." |

### Mouse Neurodegenerative Datasets (7)

| GEO ID   | Condition                          | Tissue              | Platform | Control Pattern | Disease Pattern |
|----------|------------------------------------|----------------------|----------|-----------------|-----------------|
| GSE23182 | Chronic Neurodegeneration + LPS    | Brain (hippocampus) | GPL1261  | "NBH"           | "ME7"           |
| GSE26600 | MAM Genotoxin Exposure             | Brain                | GPL1261  | "VEH"           | "MAM"           |
| GSE47516 | Cstb KO (Epilepsy/Neurodegeneration)| Cerebellum          | GPL1261  | "wild type"     | "knockout"      |
| GSE52118 | ALS Motor Pool Vulnerability       | Motor neurons        | GPL1261  | "oculomotor\|Dorsal" | "lumbar"   |
| GSE5786  | PGC-1α KO (Huntington's model)    | Striatum             | GPL1261  | "WT"            | "KO"            |
| GSE63012 | IκBα KO (Neuroinflammation)        | Hippocampus          | GPL1261  | "Control"       | "deficient"     |
| GSE31458 | MPTP (Parkinson's model)           | Brain (PFC/Putamen)  | GPL8321  | "control"       | "MPTP"          |

### Aging Databases (5)

| File | Database | Content |
|------|----------|---------|
| cellage3.tsv | CellAge | 949 cellular senescence genes |
| genage_human.csv | GenAge | 307 human aging-associated genes |
| longevity.csv | LongevityMap | 550 longevity variant associations |
| anage_data.txt | AnAge | 4645 species-level longevity records (not gene-level) |
| Search Signatures of Cellular Senescence.csv | SenMayo/Senescence Signatures | 1259 senescence marker genes |

### E-MEXP-2280 (Not Used in Main Analysis)

Downloaded raw CEL files from ArrayExpress/BioStudies. No series matrix available — would require CEL file processing with R/affy or oligo packages for future use.

---

## KIBRA-Regulated Pathways Analyzed

1. **Hippo Signaling**: LATS1/2, YAP1, WWTR1 (TAZ), STK3/4 (MST1/2), SAV1, MOB1A/B, NF2, AMOT
2. **Memory / Synaptic**: PRKCZ (PKMζ), CAMK2A, CREB1, BDNF, ARC, GRIN1, GRIN2B, DLG4 (PSD-95), SYP
3. **Autophagy / Proteostasis**: BECN1, ATG5, ATG7, SQSTM1 (p62), MAP1LC3B, LAMP2, TFEB, HSPA5, HSP90AA1
4. **PI3K-AKT-mTOR**: AKT1, MTOR, PIK3CA, PTEN, RPS6KB1, TSC1/2, RPTOR, RICTOR

## Hallmarks of Aging Analyzed (9)

Genomic Instability, Telomere Attrition, Epigenetic Alterations, Loss of Proteostasis, Deregulated Nutrient Sensing, Mitochondrial Dysfunction, Cellular Senescence, Stem Cell Exhaustion, Altered Intercellular Communication

---

## Key Findings Summary

### WWC1/Wwc1 (KIBRA)
- **Predominant trend**: downregulation in neurodegenerative disease models
- **Strongest signal**: Mouse ALS motor neurons — Wwc1 log2FC = −0.749, **p = 0.003**
- Human datasets show minimal WWC1 changes (blood-based tissue, not brain)
- WWC1 slightly upregulated in MS and HD blood cells (non-significant)

### WWC2/Wwc2
- Regulated **independently** from WWC1 across datasets
- **Significantly upregulated** in mouse chronic neurodegeneration hippocampus (p = 0.005)
- **Significantly downregulated** in mouse ALS motor neurons (p = 0.023) and PGC-1α KO striatum (p = 0.039)
- Human WWC2 upregulated in Parkinsonism model (log2FC = 0.37, trending)

### WWC3
- Only present on human platforms
- **Significantly upregulated** in Parkinsonism in vitro model (p = 0.006)

### KIBRA in Aging Databases
- **NOT found** in any of the 5 aging databases (CellAge, GenAge, LongevityMap, AnAge, Senescence Signatures)
- This is itself a notable finding — KIBRA is not yet catalogued as an aging gene despite its roles in memory consolidation and Hippo signaling

### Hallmarks of Aging Interpretation (Fig 6)
- Each cell = average Log2FC of representative hallmark genes (disease vs control)
- Red = upregulated in disease, Blue = downregulated
- **Cstb KO (epilepsy, mouse)**: most extreme — genomic instability massively UP (+6.52), nutrient sensing UP (+3.82)
- **ALS motor neurons (mouse)**: broad downregulation — stem cell exhaustion (−0.61), nutrient sensing (−0.63)
- **Chronic neurodegeneration (mouse)**: cellular senescence (+0.73), stem cell exhaustion (+0.73), altered communication (+0.86) — classic aging signatures
- **Human datasets**: smaller effects (blood tissue, not brain)
- **Multiple Sclerosis**: altered intercellular communication UP (+0.36) — consistent with inflammation

---

## Methodology Notes

- **Expression extraction**: GEO series matrix files parsed directly; expression values averaged across multiple probes per gene
- **Fold change**: Disease mean / Control mean → log2 transformed
- **Statistics**: Welch's two-sample t-test (scipy.stats.ttest_ind) when n ≥ 2 per group
- **Platform annotation parsing**: GEO .annot files have `!platform_table_begin` / `!platform_table_end` delimiters; `Gene symbol` column used; multi-gene annotations split on `///`
- **Sample classification**: Regex pattern matching on `!Sample_title` fields
- **Pathway genes**: Literature-curated representative gene sets (not exhaustive)
- **Hallmark scores**: Mean log2FC across all representative genes per hallmark per dataset

### Limitations
- Human datasets use **peripheral blood** (PBMCs/lymphocytes), not brain tissue — expression changes are attenuated
- E-MEXP-2280 (ArrayExpress) lacks a series matrix — raw CEL files need R/affy processing
- Small sample sizes in some datasets (n=3 per group)
- No multiple testing correction applied (preliminary exploration)
- Pathway gene sets are representative, not comprehensive

---

## Dependencies

```
python3 (3.9+)
numpy, pandas, scipy, matplotlib, seaborn
requests (for downloads)
```

Install with: `pip3 install --user numpy pandas scipy matplotlib seaborn requests`

(Note: Salesforce machines may need `--cert /private/etc/ssl/cert.pem` due to TLS proxy config)

---

## Next Steps / Future Questions

1. Does WWC1 dysregulation correlate with specific hallmarks of aging?
2. Are KIBRA pathway changes consistent across neurodegenerative diseases?
3. How does tissue-specificity affect WWC1/WWC2 expression patterns?
4. Can KIBRA serve as a biomarker bridging neurodegeneration and aging?
5. Process E-MEXP-2280 CEL files with R/affy for additional Parkinson's data
6. Add more brain-tissue human datasets (e.g., Allen Brain Atlas, GTEx aging)
7. Perform GSEA or pathway enrichment correlating with WWC1 expression
8. Investigate KIBRA protein interactome overlap with aging pathways
9. Look at KIBRA expression in single-cell RNA-seq neurodegenerative datasets
