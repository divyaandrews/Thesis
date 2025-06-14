# Chapter 4: A Novel Geometric INGARCH Model

## Abstract
This paper introduces an integer-valued generalized autoregressive conditional heteroskedasticity (INGARCH) model based on the novel geometric distribution and discusses some of its properties. The parameter estimation problem of the models are studied by conditional maximum likelihood and Bayesian approach using Hamiltonian Monte Carlo (HMC) algorithm. The results of the simulation studies and real data analysis affirm the good performance of the estimators and the model.

## Code Contents

This folder contains the following scripts and files:



### **1. Simulation Studies**
- `Simulation_CMLE.R`: Implements parameter estimation using Conditional Maximum Likelihood approach.
- `MATLAB_Bayesian`: Consists of MATLAB files that aid in computing parameter estimates using the Bayesian method.

### **2. Hepatitis Data Analysis**
- `hep_cml.R`: Applies the NoGe-INGARCH model to the weekly counts of Hepatitis - B cases and estimates CMLEs.
- `hep_hmc.R`: Applies the NoGe-INGARCH model to the Hepatitis data and Bayesian estimates using HMC algorithm.
- `hep_noge.stan`: Stan file required for implementation of HMC algorithm.
- `hepb.csv`: Contains the  dataset used for real data analysis.

### **3. Transactions Data Analysis**
- `trans_cmle.R`: Applies the NoGe-INGARCH model to the Transactions counts and estimates CMLEs.
- `trans_hmc.R`: Applies the NoGe-INGARCH model to the Transactions counts and Bayesian estimates using HMC algorithm.
- `trans_noge.stan`: Stan file required for implementation of HMC algorithm.
- `EricssonB_Jul2.txt`: Contains the  dataset used for real data analysis.
