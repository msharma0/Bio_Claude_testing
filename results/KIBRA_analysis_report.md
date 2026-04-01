# KIBRA/WWC1/WWC2 Expression Analysis Report
## Neurodegenerative Disease + Aging Datasets

---

## 1. Overview

- **Total datasets analyzed**: 10 neurodegenerative (GEO) + 5 aging databases
- **Human neurodeg datasets**: 3
- **Mouse neurodeg datasets**: 7
- **WWC gene-dataset comparisons**: 23

## 2. WWC1/WWC2 Expression — Human Datasets

| Dataset | Gene | Condition | Tissue | Log2FC | Direction | P-value |
|---------|------|-----------|--------|--------|-----------|---------|
| GSE21942 | WWC1 | Multiple Sclerosis | PBMCs | 0.073 | UP | 0.1189 |
| GSE21942 | WWC2 | Multiple Sclerosis | PBMCs | -0.046 | DOWN | 0.6776 |
| GSE21942 | WWC3 | Multiple Sclerosis | PBMCs | -0.016 | DOWN | 0.8168 |
| GSE8762 | WWC1 | Huntington's Disease | Lymphocytes | 0.013 | UP | 0.6652 |
| GSE8762 | WWC2 | Huntington's Disease | Lymphocytes | -0.027 | DOWN | 0.4270 |
| GSE8762 | WWC3 | Huntington's Disease | Lymphocytes | -0.022 | DOWN | 0.5432 |
| GSE35642 | WWC1 | Parkinsonism (in vitro) | SK-N-MC cells | -0.031 | DOWN | 0.7492 |
| GSE35642 | WWC2 | Parkinsonism (in vitro) | SK-N-MC cells | 0.370 | UP | 0.1028 |
| GSE35642 | WWC3 | Parkinsonism (in vitro) | SK-N-MC cells | 0.307 | UP | 0.0057 |

## 3. Wwc1/Wwc2 Expression — Mouse Datasets

| Dataset | Gene | Condition | Tissue | Log2FC | Direction | P-value |
|---------|------|-----------|--------|--------|-----------|---------|
| GSE23182 | Wwc1 | Chronic Neurodegeneration + Inflammation | Brain (hippocampus) | -0.146 | DOWN | 0.1292 |
| GSE23182 | Wwc2 | Chronic Neurodegeneration + Inflammation | Brain (hippocampus) | 0.234 | UP | 0.0051 |
| GSE26600 | Wwc1 | MAM Genotoxin Exposure | Brain | -0.001 | DOWN | 0.8479 |
| GSE26600 | Wwc2 | MAM Genotoxin Exposure | Brain | 0.002 | UP | 0.7004 |
| GSE47516 | Wwc1 | Cstb KO (Epilepsy/Neurodegeneration) | Cerebellum | 3.149 | UP | 0.4162 |
| GSE47516 | Wwc2 | Cstb KO (Epilepsy/Neurodegeneration) | Cerebellum | nan | NC | 0.0477 |
| GSE52118 | Wwc1 | ALS Motor Pool Vulnerability | Motor neurons | -0.749 | DOWN | 0.0026 |
| GSE52118 | Wwc2 | ALS Motor Pool Vulnerability | Motor neurons | -0.647 | DOWN | 0.0230 |
| GSE5786 | Wwc1 | PGC-1α KO (Huntington's model) | Striatum | -0.032 | DOWN | 0.8126 |
| GSE5786 | Wwc2 | PGC-1α KO (Huntington's model) | Striatum | -0.127 | DOWN | 0.0394 |
| GSE63012 | Wwc1 | IκBα KO (Neuroinflammation) | Hippocampus | 0.004 | UP | 0.5360 |
| GSE63012 | Wwc2 | IκBα KO (Neuroinflammation) | Hippocampus | 0.002 | UP | 0.8221 |
| GSE31458 | Wwc1 | MPTP (Parkinson's model) | Brain (PFC/Putamen) | -0.012 | DOWN | 0.4052 |
| GSE31458 | Wwc2 | MPTP (Parkinson's model) | Brain (PFC/Putamen) | -0.021 | DOWN | 0.0872 |

## 4. Cross-Species Comparison

### WWC1/Wwc1 Summary
- **Human**: 3 datasets — directions: UP, UP, DOWN
- **Mouse**: 7 datasets — directions: DOWN, DOWN, UP, DOWN, DOWN, UP, DOWN

### WWC2/Wwc2 Summary
- **Human**: 3 datasets — directions: DOWN, DOWN, UP
- **Mouse**: 7 datasets — directions: UP, UP, NC, DOWN, DOWN, UP, DOWN

## 5. KIBRA in Aging Databases

**WWC1/WWC2/KIBRA was NOT found in any of the aging databases** (CellAge, GenAge, Longevity, Senescence Signatures).

This is notable — despite KIBRA's known role in memory and Hippo signaling (both aging-relevant), it has not yet been catalogued as a canonical aging gene.

## 6. KIBRA Pathway Analysis

### Hippo Signaling
- Gene-dataset pairs: 108 (UP: 55, DOWN: 43)
  - Human: 31 pairs (UP: 19, DOWN: 12)
  - Mouse: 77 pairs (UP: 36, DOWN: 31)

### Memory Synaptic
- Gene-dataset pairs: 90 (UP: 35, DOWN: 47)
  - Human: 27 pairs (UP: 13, DOWN: 14)
  - Mouse: 63 pairs (UP: 22, DOWN: 33)

### Autophagy Proteostasis
- Gene-dataset pairs: 90 (UP: 41, DOWN: 41)
  - Human: 27 pairs (UP: 8, DOWN: 19)
  - Mouse: 63 pairs (UP: 33, DOWN: 22)

### PI3K AKT mTOR
- Gene-dataset pairs: 87 (UP: 35, DOWN: 46)
  - Human: 25 pairs (UP: 13, DOWN: 12)
  - Mouse: 62 pairs (UP: 22, DOWN: 34)

## 7. Hallmarks of Aging Expression

### Genomic Instability
- Gene-dataset pairs: 60 (UP: 26, DOWN: 29)
  - Human: 18 pairs (UP: 8, DOWN: 10)
  - Mouse: 42 pairs (UP: 18, DOWN: 19)

### Telomere Attrition
- Gene-dataset pairs: 37 (UP: 16, DOWN: 20)
  - Human: 9 pairs (UP: 4, DOWN: 5)
  - Mouse: 28 pairs (UP: 12, DOWN: 15)

### Epigenetic Alterations
- Gene-dataset pairs: 70 (UP: 40, DOWN: 23)
  - Human: 21 pairs (UP: 13, DOWN: 8)
  - Mouse: 49 pairs (UP: 27, DOWN: 15)

### Loss of Proteostasis
- Gene-dataset pairs: 60 (UP: 24, DOWN: 32)
  - Human: 18 pairs (UP: 7, DOWN: 11)
  - Mouse: 42 pairs (UP: 17, DOWN: 21)

### Nutrient Sensing
- Gene-dataset pairs: 65 (UP: 28, DOWN: 32)
  - Human: 18 pairs (UP: 9, DOWN: 9)
  - Mouse: 47 pairs (UP: 19, DOWN: 23)

### Mitochondrial Dysfunction
- Gene-dataset pairs: 50 (UP: 21, DOWN: 26)
  - Human: 15 pairs (UP: 9, DOWN: 6)
  - Mouse: 35 pairs (UP: 12, DOWN: 20)

### Cellular Senescence
- Gene-dataset pairs: 60 (UP: 29, DOWN: 25)
  - Human: 18 pairs (UP: 8, DOWN: 10)
  - Mouse: 42 pairs (UP: 21, DOWN: 15)

### Stem Cell Exhaustion
- Gene-dataset pairs: 60 (UP: 28, DOWN: 29)
  - Human: 18 pairs (UP: 11, DOWN: 7)
  - Mouse: 42 pairs (UP: 17, DOWN: 22)

### Altered Intercellular Comm
- Gene-dataset pairs: 60 (UP: 37, DOWN: 17)
  - Human: 18 pairs (UP: 13, DOWN: 5)
  - Mouse: 42 pairs (UP: 24, DOWN: 12)

## 8. Preliminary Inferences

### Key Findings
1. **WWC1 (KIBRA) expression** in neurodegenerative contexts varies by disease and species
2. **WWC2** shows independent regulation from WWC1 in several datasets
3. **KIBRA is not yet classified** as an aging gene in major databases
4. **KIBRA-regulated pathways** (Hippo, synaptic, autophagy) show disease-specific dysregulation
5. **Hallmarks of aging genes** show overlapping changes with neurodegenerative conditions

### Next Questions to Explore
- Does WWC1 dysregulation correlate with specific hallmarks of aging?
- Are KIBRA pathway changes consistent across neurodegenerative diseases?
- How does tissue-specificity affect WWC1/WWC2 expression patterns?
- Can KIBRA serve as a biomarker bridging neurodegeneration and aging?

---

## Output Files

- `Table1_WWC_Human_Expression.csv` — WWC1/WWC2 expression in human datasets
- `Table1_WWC_Mouse_Expression.csv` — Wwc1/Wwc2 expression in mouse datasets
- `Table1_WWC_All_Expression.csv` — Combined expression table
- `Table2_KIBRA_Pathway_Expression.csv` — KIBRA pathway gene expression
- `Table3_Hallmarks_Aging_Expression.csv` — Hallmarks of aging gene expression
- `Table4_KIBRA_in_Aging_Databases.csv` — KIBRA presence in aging databases
- `Fig1-7` — Publication-quality figures (PNG + PDF)