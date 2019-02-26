
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
========

If multiple tissue or cell-type specific biological pathways underlie a phenotype, it can be possible to identify the tissue-specific genetic subtype of the phenotype which is defined by heterogeneous tissue-specific genetic architecture of the trait. For example, an individual's body mass index (BMI) can be regulated more by the set of genes higher expressed in brain compared to the set of genes higher expressed in adipose and vice versa. We propose a method eGST which learns about such subtype structure of a complex trait based on individual-level data of marginal phenotype and genotypes of sets of expression quantitative trait loci (eQTLs), each corresponding to the set of genes higher expressed in a tissue. The method is based on a finite mixture model under Bayesian framework and implements the maximum a posteriori (MAP) expectation-maximization (EM) algorithm for estimation of subtype posterior probability across individuals and other parameters in the model.

The package consists of the following function:

1.  eGST: It estimates the posterior probability that an individual's phenotype is a tissue-specific genetic subtype. The phenotype across individuals can be classified as a tissue-specific genetic subtype based on the estimated subtype posterior probability across individuals.

Installation
============

You can install eGST from CRAN.

``` r
install.packages("eGST")
library("eGST")
```

How to run eGST for two tissues.
================================

Get the path to the data.

``` r
library("eGST")
# Load the phenotype data vector
phenofile <- system.file("data", "ExamplePhenoData.rda", package = "eGST")
head(ExamplePhenoData)
```

Here ExamplePhenoData is the phenotype data vector for 1000 individuals.

``` r
library("eGST")
# Load the list containing genotype matrices of tissue-specific eQTLs. 
genofile <- system.file("data", "ExampleEQTLgenoData.rda", package = "eGST")
ExampleEQTLgenoData[[1]][1:5, 1:5]
```

Here ExampleEQTLgenoData is a list containing two elements corresponding to two tissues each containing a 1000 by 100 ordered genotype matrix. Each matrix provides the genotype data of 1000 individuals at 100 tissue-specific eQTLs for each tissue. Here we have displayed genotypes for first 5 individuals at first 5 eQTLs in the set of first tissue-specific eQTLs. We normalize each SNP's genotype data across all individuals in the sample before running eGST.

Next we specify the name of the tissues.

``` r
# Specify the name of the tissues.
tissues <- paste("tissue", 1:2, sep = "")
tissues
```

In this simulated example dataset, we have considered two tissues and corresponding sets of 100 tissue-specific eQTLs each. First half of 1000 individuals' phenotypes were simulated to have genetic effect from the first tissue specific eQTLs, but no effect from the second tissue-specific eQTLs. Hence first 500 individuals were assigned as the first tissue-specific genetic subtype. Similarly, second half of the 1000 individuals were simulated to have genetic effct from the second-tissue specific eQTLs and hence they are second tissue-specific genetic subtype.

Next for this toy example dataset, we run eGST for 10 iterations. However, we recommend at least 50 iterations in your application. There are more options of arguments to pass into the function (see the Arguments section of eGST in the eGST manual).

``` r
#Run eGST to estimate the subtype posterior probability across individuals.
result <- eGST(ExamplePhenoData, ExampleEQTLgenoData, tissues, nIter = 10)
```

So at each iteration, eGST prints the average improvement in log-likelihood of the data. Next we display an overall summary of the results obtained by eGST.

``` r
# Overall summary of the results produced by eGST.
str(result)
```

The main output of interest is contained in result$gamma matrix which provides the estimate of tissue-specific subtype posterior probability across individuals. So the first (second) column of result$gamma provides the posterior probability that an individual's phenotype is the first (second) tissue-specific genetic subtype. Individuals can be classified as tissue-specific genetic subtype of the trait based on a posterior probability threshold, e.g. 65%, 70%, etc. For example, individuals for whom the first tissue-specific posterior probability is &gt; 65% can be assigned as first tissue-specific genetic subtype. The list 'result' contains other outputs from eGST. For more details, please see the 'Value' section in eGST manual.

Getting more details
====================

For any questions, please send an email to <statgen.arunabha@gmail.com> or <pasaniuc@ucla.edu>. See our manuscript for more details: Majumdar A, Giambartolomei C, Cai N, Haldar T, Freund MK, Pasaniuc B (2019) Leveraging eQTLs to identify tissue-specific genetic subtype of complex trait. bioRxiv.
