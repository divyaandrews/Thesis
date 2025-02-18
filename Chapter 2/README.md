# Chapter 2: A Novel Geometric AR(1) Model

## Abstract
This paper proposes a first-order geometric autoregressive model based on a new thinning operator. The parameters of the model are estimated by the method of Two-Stage Conditional Least Squares and Saddlepoint Approximation. Simulation studies demonstrate that the saddlepoint approximation estimator performs well, even in the case of small sample sizes. To illustrate the significance of the model, a real data analysis of the number of downloads of a TeX editor is presented.

## Code Contents

This folder contains the following scripts and files:

### **1. Model Implementation**
- `NoGeAR.R`: Simulates data from the proposed Geometric AR(1) model.


### **2. Simulation Studies**
- `NoGeAR1_estimation.R`: Implements parameter estimation using Two-Stage Conditional Least Squares.
- `NoGeAR1_saddlepoint.R`: Computes parameter estimates using the Saddlepoint Approximation method.
- `NoGeAR1_simulation_study.R`: Conducts extensive simulation studies comparing estimator performance.

### **3. Real Data Analysis**
- `NoGeAR1_real_data_analysis.R`: Applies the Geometric AR(1) model to the number of downloads of a TeX editor.


### **5. Supporting Files**
- `downloads.txt`: Contains the  dataset used for real data analysis.
- `README.md`: This file, providing an overview of the available codes.

For details on how to run the scripts, please refer to the comments within each file or reach out for further assistance.
