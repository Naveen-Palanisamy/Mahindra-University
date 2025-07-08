# 📊 Response Surface Methodology (RSM) in Bioprocess Optimization

This project uses R's `rsm` package to build and analyze a Response Surface Model (RSM) to optimize yield in a bioprocess system.

## 🔬 What is RSM?

Response Surface Methodology (RSM) is a set of statistical techniques for modeling and optimizing processes influenced by several variables. It is especially powerful in **bioprocess engineering**, where experiments can be complex and costly.

### ✅ Why use RSM in Bioprocessing?

- Efficient optimization of process parameters (e.g., pH, temperature, nutrients)
- Detection of interaction effects between factors
- Reduction in number of required experiments
- Prediction of optimal yield conditions using polynomial models

### 🧪 Project Description

This R script:
- Loads experimental data from Excel
- Fits a second-order RSM model: `Yield ~ SO(x1, x2, x3, x4)`
- Visualizes yield response using contour and 3D plots
- Predicts optimal yield based on specific conditions

> Designed for experimental design and optimization tasks in biological or biochemical research.

## 📁 Files

- `rsm_model_analysis.R` – Main R script for model fitting and visualization
- `README.md` – Project overview and background

---

🧬 Made for bioprocess engineers, scientists, and R users interested in process optimization.
