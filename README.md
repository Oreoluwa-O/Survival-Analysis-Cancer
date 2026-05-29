# Survival Analysis of Cancer Patient Data

## Overview
Cancer remains one of the leading causes of death in the United States, with over 600,000 deaths annually. 
This project applies survival analysis techniques to publicly available cancer patient data from The Cancer 
Genome Atlas (TCGA) to identify factors that significantly affect patient survival outcomes.

The findings from this analysis contribute to the broader goal of strengthening medical research 
infrastructure in the United States through rigorous biostatistical methods.

## Public Health Significance
- Survival analysis is a core method in oncology clinical trials and drug approval processes in the US
- Understanding predictors of survival informs treatment decisions for millions of American patients
- This project demonstrates reproducible biostatistical workflows applicable to real clinical research

## Methods
- **Kaplan-Meier estimator** — visualizes survival probability over time across patient groups
- **Log-rank test** — statistically compares survival curves between groups (e.g. by cancer stage)
- **Cox Proportional Hazards model** — identifies which clinical variables significantly predict survival

## Data Source
Data is sourced from the [TCGA Research Network](https://www.cancer.gov/tcga), a publicly available 
landmark cancer genomics program funded by the US National Cancer Institute (NCI).

## Requirements
## How to Run
```bash
pip install -r requirements.txt
python src/survival_analysis.py
```

## Results
Key outputs include:
- Kaplan-Meier survival curves by cancer stage
- Hazard ratios from Cox regression with confidence intervals
- Statistical significance of clinical predictors on patient survival

## Author
**Oreoluwa Oyetubo** — M.S. Statistics  
Biostatistics & Health Analytics Researcher  
[GitHub Profile](https://github.com/Oreoluwa-O)
