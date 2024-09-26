# Protein-Peptide Docking Using ESMFold with Recycling

This repository contains code, data, and scripts related to our research on **protein-peptide docking using the ESMFold language model**. ESMFold is primarily designed for protein structure prediction, but here we explore its potential for docking peptides to proteins using a modified approach that includes **recycling** and **masking** strategies. These modifications are aimed at improving sampling diversity and docking accuracy.

---

## Contents

- **`Protein_peptide_docking_using_ESMFold_with_recycling.ipynb`**  
  This Jupyter notebook is an edited version of the ESMFold model, available via Colab. It enables protein-peptide docking by introducing a **recycling** mechanism for repeated model predictions and **masking** to increase model diversity. Detailed explanations of these methods can be found in our manuscript.
  
- **`weighted_avg_plddt.py`**  
  This Python script calculates the **weighted average pLDDT score**, a metric that evaluates the confidence of the model's predictions. Higher pLDDT values indicate greater confidence in predicted structures. The script is used to score the results of docking simulations.
  
- **Zipped Data Files**  
  Several zipped files contain the results of different protein-peptide docking simulations. These files are generated using ESMFold under different settings, which are outlined in the table below:
  
  | File Name                | Description |
  |--------------------------|-------------|
  | `default.zip`             | Default docking results without any masking or additional recycling. Represents baseline predictions using the standard ESMFold method with a 30-residue poly-glycine linker. |
  | `masking-c-200.zip`       | Docking results with **random masking** applied to the peptide sequence and a 200-residue poly-glycine linker at the C-terminal end of the peptide. |
  | `masking-c-30.zip`        | Docking results with random masking and a **30-residue poly-glycine linker** at the C-terminal end, which was found to be the optimal configuration for this study. |
  | `masking-n+c-100.zip`     | Docking results with random masking and **dual 100-residue poly-glycine linkers** at both the N- and C-terminal ends of the peptide. |
  | `masking-n-200.zip`       | Docking results with random masking and a **200-residue poly-glycine linker** at the N-terminal end of the peptide. |
  | `recycles-200.zip`        | Docking results using **adaptive recycling**, with up to 12 recycles and a **200-residue linker** to increase peptide-receptor contact. |
  | `recycles-30.zip`         | Docking results using adaptive recycling and a **30-residue poly-glycine linker**. |

---
