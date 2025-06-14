# Chapter 4: A Novel Geometric INGARCH Model

## Abstract
This paper introduces an integer-valued generalized autoregressive conditional heteroskedasticity (INGARCH) model based on the novel geometric distribution and discusses some of its properties. The parameter estimation problem of the models are studied by conditional maximum likelihood and Bayesian approach using Hamiltonian Monte Carlo (HMC) algorithm. The results of the simulation studies and real data analysis affirm the good performance of the estimators and the model.

## Code Contents

This folder contains the following scripts and files:



### **1. Simulation Studies**
- `NoGeAR1_estimation.R`: Implements parameter estimation using Two-Stage Conditional Least Squares.
- `NoGeAR1_saddlepoint.R`: Computes parameter estimates using the Saddlepoint Approximation method.
- `NoGeAR1_simulation_study.R`: Conducts extensive simulation studies comparing estimator performance.

### **2. Hepatitis Data Analysis**
- `NoGeAR1_data_analysis.R`: Applies the NoGeAR(1) model to the number of downloads of a TeX editor.


### **3. Transactions Data Analysis**
- `Downloads.txt`: Contains the  dataset used for real data analysis.
- `README.md`: This file, providing an overview of the available codes.

