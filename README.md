# privateEC

## Differential privacy-based Evaporative Cooling feature selection and classification with Relief-F and Random Forests

Written by Trang Le with Bill White
[University of Tulsa McKinney Lab](http://insilico.utulsa.edu)

Methods are described in the following publication.

Trang Le, W. K. Simmons, M. Misaki, B.C. White, J. Savitz, J. Bodurka, and B. A. McKinney. 
“Differential privacy-based Evaporative Cooling feature selection and classification with Relief-F and Random Forests,” 
Bioinformatics. Accepted. https://doi.org/10.1093/bioinformatics/btx298. 2017.

# Installation

```
library(devtools)
install_github("insilico/privateEC") 
library(privateEC)
vignette("Example1")
```

# Abstract

## Motivation

Classification of individuals into disease or clinical categories from high-dimensional biological data with low prediction error is an important challenge of statistical learning in bioinformatics. Feature selection can improve classification accuracy but must be incorporated carefully into cross-validation to avoid overfitting. Recently, feature selection methods based on differential privacy, such as differentially private random forests and reusable holdout sets, have been proposed. However, for domains such as bioinformatics, where the number of features is much larger than number of observations, these differential privacy methods are susceptible to overfitting.

## Methods

We introduce private Evaporative Cooling, a stochastic privacy-preserving machine learning algorithm that uses Relief-F for feature selection and random forest for privacy preserving classification that also prevents overfitting. We relate the privacy-preserving threshold mechanism to a thermodynamic Maxwell-Boltzmann distribution, where the temperature represents the privacy threshold. We use the thermal statistical physics concept of Evaporative Cooling of atomic gases to perform backward stepwise privacy-preserving feature selection.

## Results

On simulated data with main effects and statistical interactions, we compare accuracies on holdout and validation sets for three privacy-preserving methods: the reusable holdout, reusable holdout with random forest, and private Evaporative Cooling, which uses Relief-F feature selection and random forest classification. In simulations where interactions exist between attributes, private Evaporative Cooling provides higher classification accuracy without overfitting based on an independent validation set. In simulations without interactions, thresholdout with random forest and private Evaporative Cooling give comparable accuracies. We also apply these privacy methods to human brain resting-state fMRI data from a study of major depressive disorder.
