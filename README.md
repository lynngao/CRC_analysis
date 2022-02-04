# Increasing prediction performance of colorectal cancer disease status using random forests classification based on metagenomic shotgun sequencing data

#### Yilin Gao, Zifan Zhu, Fengzhu Sun, Increasing prediction performance of colorectal cancer disease status using random forests classification based on metagenomic shotgun sequencing data, Synthetic and Systems Biotechnology, Volume 7, Issue 1, 2022, Pages 574-585, ISSN 2405-805X, https://doi.org/10.1016/j.synbio.2022.01.005.

This repository contains R scripts for running the analysis described in our manuscript. 
The R scripts can be divided into 4 basic components:

#### 1. Training Random forests, LASSO and SVM classifiers on within-dataset, cross-dataset and leave-one-dataset-out experimental settings
- [3classifiers.R](https://github.com/lynngao/CRC_analysis/blob/main/3classifiers.R): helper file contains functions for theses three classifiers.
- [3classifiers_example.R](https://github.com/lynngao/CRC_analysis/blob/main/3classifiers_example.R): one example about how to run the three classifiers using the helper file.

#### 2. Leave-One-Sample-Out(LODO) model stacking algorithm on within-dataset, cross-dataset and leave-one-dataset-out experimental settings
- [LOSO_model_stacking.R](https://github.com/lynngao/CRC_analysis/blob/main/LOSO_model_stacking.R): helper file contains functions for implementing LOSO algorithm.
- [LOSO_model_stacking_example.R](https://github.com/lynngao/CRC_analysis/blob/main/LOSO_model_stacking_example.R): one example about how to run LOSO model stacking method using the helper file.

#### 3. Calculating AUPRC
- [AUPRC.R](https://github.com/lynngao/CRC_analysis/blob/main/AUPRC.R): helper file contains functions for calculating AUPRC values in different settings.
- [AUPRC_example.R](https://github.com/lynngao/CRC_analysis/blob/main/AUPRC_example.R): one example about how to calculate AUPRC values using the helper file.

#### 4. Removing batch effect by ComBat function before training random forests on cross-dataset setting
- [ComBat.R](https://github.com/lynngao/CRC_analysis/blob/main/ComBat.R): helper file contains functions for running ComBat before random forests.
- [ComBat_example.R](https://github.com/lynngao/CRC_analysis/blob/main/ComBat_example.R): one example about how to run ComBat using the helper file.

For running tasks 1 and 4, microbial species abundance profiles are required. See abundance profiles in [abundance](https://github.com/lynngao/CRC_analysis/tree/main/abundance) as example of the format of required profiles.

For For running tasks 2 and 3, both microbial species abundance profiles and prediction probability files are required. See prediction probability profiles in [pred_prob](https://github.com/lynngao/CRC_analysis/tree/main/pred_prob) as example of the format of required files.
