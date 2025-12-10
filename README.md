# MDS DTA Mutational Burden Study

This repository contains the code, SQL queries, and Jupyter notebooks used for the computational analysis of **DTA mutational burden (DNMT3A, TET2, ASXL1)** as a prognostic biomarker in **myelodysplastic syndromes (MDS)**.  
The project evaluates the impact of cumulative mutational burden on **overall survival (OS)** and **leukemia-free survival (LFS)** using both traditional survival methods and machine-learning models.

---

##  Background

Mutations in **DNMT3A**, **TET2**, and **ASXL1** (collectively known as **DTA**) are highly prevalent in:

- Age-related clonal hematopoiesis (ARCH / CHIP)  
- Myelodysplastic syndromes (MDS)  
- Myeloid neoplasms more broadly

While individual DTA mutations have been associated with CHIP biology and clonal evolution, the **prognostic impact of cumulative DTA mutational burden** in treatment-naïve MDS has not been fully explored.  
This project evaluates whether *the number and type* of DTA variants contribute additional prognostic information beyond established clinical tools.

---


##  Methods Summary

### **Data Processing**
- SQL-based extraction of clinical, cytogenetic, and genomic features  
- Harmonization of mutation calls across DTA loci  
- Creation of cumulative DTA burden metrics:
  - `n_dta` (total DTA mutations)
  - gene-specific burden (e.g., `n_asxl1`, `n_dnmt3a`)
  - truncating vs. non-truncating subclassification

### **Statistical Modeling**
- Kaplan–Meier survival analysis  
- Cox proportional hazards models  
- Multiple comparison adjustments (Holm–Bonferroni)  
- Bootstrap confidence intervals for C-index and IBS

### **Machine Learning Survival Models**
- **CoxNet** (elastic net penalized Cox regression)  
- **Random Survival Forest (RSF)**  
- **Gradient Boosting Survival Tree (GBST)**  
- SHAP explainability for feature attribution  

### **Model Evaluation**
- C-index with 95% bootstrap confidence intervals  
- Integrated Brier Score (IBS) for calibration  
- SHAP-based variable importance for interpretability  

---

### **Clone the repository**
```bash
git clone git@github.com:pvhernandezlab/MDS_DTA_mutational_burden_study.git
cd MDS_DTA_mutational_burden_study
```

##  Contact

For questions, collaboration inquiries, or discussion about the analysis, please contact:

**Patricia Hernandez, MD**
Molecular Genetic Pathology & Clinical Informatics  
Email: pvhernandezlab@users.noreply.github.com
GitHub: https://github.com/pvhernandezlab


