# Federated GWAS — Data Description

This document describes the synthetic federated GWAS dataset used in the thesis experiments.  
The dataset consists of **500 individuals** spread across **6 simulated research sites** (clients). This data is inspired from the structure of **toy_dataset** from `sPLINK` paper

---

## 1. Overall Design

| Property | Value |
|---|---|
| Total individuals | 500 |
| Number of clients | 6 |
| SNP panel | 100 autosomal SNPs (assumed after QC phase: MAF > 0.05, no missing, LD-pruned) |
| Covariates | Sex (1 = Male, 2 = Female), Age (years, continuous), Smoking_Status (0/1, 0 = former, 1 = current) |
| Phenotype | Binary case/control (0 = control, 1 = case; PLINK FAM: 1 = control, 2 = case) |


### Covariate distributions (population-level targets)

| Covariate | Distribution |
|---|---|
| Sex | ~50 % Male / 50 % Female |
| Age | Gaussian (μ = 62.5, σ = 7.5), clipped to [45, 80] years, rounded to 1 decimal place |
| Smoking_Status | 0 = former (50 %), 1 = current (50 %) |
| Phenotype | ~50 % cases / 50 % controls |

### SNP panel

All 100 SNPs are drawn from a shared panel spanning chromosomes 1–22.  
Minor-allele frequencies (MAFs) are generated per client by perturbing a shared base MAF vector (Uniform[0.10, 0.45]) with site-specific noise (Uniform[−0.08, +0.08]), clamped to [0.06, 0.49].  
Genotypes follow Hardy–Weinberg proportions (dosage 0 / 1 / 2).

---

## 2. File Formats

Each client directory contains five files:

| File | Format | Description |
|---|---|---|
| `client{N}.bim` | PLINK BIM | SNPs metadata |
| `client{N}.fam` | PLINK FAM | Samples metadata - FID, IID, father, mother, sex (1/2), phenotype (1=ctrl, 2=case) |
| `client{N}.cov` | COV | FID, IID, Sex (1/2), Age (float), Smoking_Status (0/1, 0 = former, 1 = current) |
| `client{N}.bed` | PLINK BED binary (SNP-major) | - |
| `client{N}_geno.txt` | TEXT FILE | Genotype matrix: rows = SNPs (100), cols = individuals; values 0/1/2 |

`pool.csv` contains the permuted pool data from all clients  with a header row:

```
IID(Identifier across all clients)  FID(Identifier within a client)  Client  <SNP1> … <SNP100>  Sex  Age  Smoking_Status  Phenotype
```

---

## 3. Client Profiles

The table below summarises each client's sample composition.  
**Sex** is coded 1 = Male, 2 = Female.  
**Smoking** columns count individuals with status 0 (former) and 1 (current).  
**IID range** gives the first and last synthetic individual identifier assigned to that site.  
**MAF** are the frequency of the less common allele at a SNP, computed from each client’s genotype matrix across all individuals..

### 3.1 Summary Table

| Client | n | IID Range | Male | Female | Cases | Controls | Age Mean ± SD | Age Range | Smoke 0 | Smoke 1 | MAF Mean | MAF Range |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| client1 | 80 | SYN0001 – SYN0080 | 41 | 39 | 43 | 37 | 62.7 ± 8.1 | 45.0 – 80.0 | 40 | 40 | 0.260 | 0.050 – 0.556 |
| client2 | 85 | SYN0081 – SYN0165 | 38 | 47 | 40 | 45 | 61.4 ± 6.8 | 45.8 – 75.2 | 38 | 47 | 0.263 | 0.041 – 0.529 |
| client3 | 90 | SYN0166 – SYN0255 | 51 | 39 | 33 | 57 | 63.9 ± 8.2 | 46.0 – 79.8 | 44 | 46 | 0.263 | 0.033 – 0.511 |
| client4 | 85 | SYN0256 – SYN0340 | 45 | 40 | 38 | 47 | 62.8 ± 7.2 | 45.0 – 80.0 | 36 | 49 | 0.265 | 0.024 – 0.547 |
| client5 | 80 | SYN0341 – SYN0420 | 46 | 34 | 46 | 34 | 62.1 ± 7.3 | 46.7 – 80.0 | 41 | 39 | 0.267 | 0.044 – 0.556 |
| client6 | 80 | SYN0421 – SYN0500 | 42 | 38 | 42 | 38 | 63.1 ± 8.0 | 45.0 – 79.0 | 47 | 33 | 0.256 | 0.031 – 0.538 |

### 3.2 Per-Client Details

#### client1 — 80 individuals (SYN0001 – SYN0080)

- **Sex:** 41 Male, 39 Female
- **Phenotype:** 43 cases, 37 controls
- **Age:** mean = 62.7, SD = 8.1, range [45.0, 80.0], Q1/median/Q3 = 57.8 / 62.4 / 67.9
- **Smoking:** 40 former, 40 current
- **MAF (100 SNPs):** mean = 0.260, range [0.050, 0.556]

---

#### client2 — 85 individuals (SYN0081 – SYN0165)

- **Sex:** 38 Male, 47 Female
- **Phenotype:** 40 cases, 45 controls
- **Age:** mean = 61.4, SD = 6.8, range [45.8, 75.2], Q1/median/Q3 = 56.0 / 60.9 / 65.8
- **Smoking:** 38 former, 47 current
- **MAF (100 SNPs):** mean = 0.263, range [0.041, 0.529]

---

#### client3 — 90 individuals (SYN0166 – SYN0255)

- **Sex:** 51 Male, 39 Female
- **Phenotype:** 33 cases, 57 controls
- **Age:** mean = 63.9, SD = 8.2, range [46.0, 79.8], Q1/median/Q3 = 56.7 / 64.7 / 70.5
- **Smoking:** 44 former, 46 current
- **MAF (100 SNPs):** mean = 0.263, range [0.033, 0.511]

---

#### client4 — 85 individuals (SYN0256 – SYN0340)

- **Sex:** 45 Male, 40 Female
- **Phenotype:** 38 cases, 47 controls
- **Age:** mean = 62.8, SD = 7.2, range [45.0, 80.0], Q1/median/Q3 = 58.4 / 62.4 / 67.5
- **Smoking:** 36 former, 49 current
- **MAF (100 SNPs):** mean = 0.265, range [0.024, 0.547]

---

#### client5 — 80 individuals (SYN0341 – SYN0420)

- **Sex:** 46 Male, 34 Female
- **Phenotype:** 46 cases, 34 controls
- **Age:** mean = 62.1, SD = 7.3, range [46.7, 80.0], Q1/median/Q3 = 56.5 / 62.2 / 67.5
- **Smoking:** 41 former, 39 current
- **MAF (100 SNPs):** mean = 0.267, range [0.044, 0.556]

---

#### client6 — 80 individuals (SYN0421 – SYN0500)

- **Sex:** 42 Male, 38 Female
- **Phenotype:** 42 cases, 38 controls
- **Age:** mean = 63.1, SD = 8.0, range [45.0, 79.0], Q1/median/Q3 = 58.0 / 64.3 / 68.6
- **Smoking:** 47 former, 33 current
- **MAF (100 SNPs):** mean = 0.256, range [0.031, 0.538]

---

## 4. Sample eigenvectors across pool data

The sample eigenvectors for the pool data are stored in `data/pool/`:

| File | Description |
|---|---|
| `G_sample_eigenvectors.csv` | PCA sample eigenvector matrix G of shape (500 × 20); columns IID, PC1 … PC20 |
| `H_SNP_eigenvectors.csv` | PCA feature eigenvector matrix G of shape (100 × 20); columns IID, PC1 … PC20 |
| `eigenvalues.csv` | PCA eigenvalues (G), with one value per principal component (PC1–PC20), representing the variance explained by each principal component. |

---
