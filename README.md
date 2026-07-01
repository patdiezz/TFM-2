# TFM 2026-Patricia 
R code for fitting and evaluating country-specific distance-based classification models for poor self-perceived health using SHARE Wave 9 data. The workflow includes Gower-based distance computation, dbglm modelling, sub-bagging, model evaluation, ROC analysis and sensitivity analysis for Spain, Germany, Finland and Latvia.

# Country-specific distance-based models for self-perceived health

This repository contains the R code used to fit, evaluate and interpret country-specific distance-based classification models for poor self-perceived health using SHARE Wave 9 data. The analysis was developed as part of a Master Thesis focused on self-perceived health in later life in four European countries: Spain, Germany, Finland and Latvia.

## Repository contents

The code includes the main steps of the modelling workflow:

- Data preprocessing and preparation of the country-specific datasets.
- Stratified train-test split by self-perceived health status.
- Construction of Gower-based distance matrices for mixed-type data.
- Implementation of distance-based generalized linear models using the `dbglm` framework.
- Sub-bagging strategy for robust G-Gower models.
- Comparison between robust G-Gower, classical Gower and balanced robust G-Gower approaches.
- Evaluation of model performance using accuracy, sensitivity, specificity, PPV, NPV, ROC curves and AUC.
- Measurement of computational time for the different distance-based approaches.
- Sensitivity analysis to interpret how each predictor affects the predicted probability of poor self-perceived health.
- Generation of country-specific plots and summary tables used in the thesis.

The models are based on a distance-based approach because the SHARE data include continuous, binary and categorical variables related to physical health, functional limitations, cognition, emotional well-being, social relationships and access to healthcare. The `dbglm` framework allows these heterogeneous predictors to be represented through distance matrices while modelling the binary outcome using a logistic structure.

## Data availability

The original SHARE microdata are not included in this repository.

According to SHARE Conditions of use (Section 7: Confidentiality of use), we are not allowed to make copies of the data available to others and/or enable any third party access to the database. Access to the SHARE data is only granted on an individual basis. This information is available at :https://share-eric.eu/data/data-access/conditions-of-use

Therefore, this repository only contains the R scripts used for the analysis. Users who wish to reproduce the results must request individual access to SHARE data through the official SHARE-ERIC data access procedure.

## Purpose
The purpose of this repository is to support transparency and reproducibility of the statistical workflow developed in the thesis, while respecting SHARE data access restrictions.

