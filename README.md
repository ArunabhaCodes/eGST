
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
========

Genetic predisposition for complex traits is often manifested through multiple tissues of interest at different time points in the development. As an example, the genetic predisposition for obesity could be manifested through inherited variants that control metabolism through regulation of genes expressed in the brain and/or through the control of fat storage in the adipose tissue by dysregulation of genes expressed in adipose tissue. We present a method eGST that integrates tissue-specific eQTLs with GWAS data for a complex trait to probabilistically assign a tissue of interest to the phenotype of each individual in the study. eGST estimates the posterior probability that an individual's phenotype can be assigned to a tissue based on individual-level genotype data of tissue-specific eQTLs and marginal phenotype data in a GWAS cohort. Under a Bayesian framework of mixture model, eGST employs a maximum a posteriori (MAP) expectation-maximization (EM) algorithm to estimate the tissue-specific posterior probability across individuals.


The package consists of the following function:

1.  eGST: It estimates the posterior probability that the genetic susceptibility of the phenotype of an individual in the study is mediated through eQTLs specific to a tissue of interest. The phenotype across individuals can be classified into tissues under consideration based on the estimated tissue-specific posterior probability across individuals.

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

Here ExampleEQTLgenoData is a list containing two elements corresponding to two tissues each containing a 1000 by 100 ordered genotype matrix. Each matrix provides the genotype data of 1000 individuals at 100 tissue-specific eQTLs for each tissue. To create sets of tissue-specific eQTLs in your context, please see our manuscript: Majumdar A, Giambartolomei C, Cai N, Freund MK, Haldar T, J Flint, Pasaniuc B (2019) Leveraging eQTLs to identify tissue-specific genetic subtype of complex trait, bioRxiv. Here we have displayed genotypes for first 5 individuals at first 5 eQTLs in the set of first tissue-specific eQTLs. We normalize each SNP's genotype data across all individuals in the sample before running eGST. 
Next we specify the name of the tissues.

``` r
# Specify the name of the tissues.
tissues <- paste("tissue", 1:2, sep = "")
tissues
```

In this simulated example dataset, we have considered two tissues
and corresponding sets of 100 tissue-specific eQTLs each.
First half of 1000 individuals' phenotypes were simulated
to have genetic effect from the first tissue
specific eQTLs, but no effect from the second
tissue-specific eQTLs. Hence the phenotype of first 500 individuals
were assigned to the first tissue. Similarly, second half of the 1000 individuals
were simulated to have genetic effect from the second-tissue
specific eQTLs.

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

The main output of interest is contained in result\$gamma matrix which provides the
estimate of tissue-specific subtype posterior probability across individuals.
So the first (second) column of result\$gamma provides the posterior probability that an 
individual's phenotype is the first (second) tissue-specific genetic subtype.
Individuals can be classified as tissue-specific genetic subtype of the trait
based on a posterior probability threshold, e.g. 65\%, 70\%, etc. For example,
individuals for whom the first tissue-specific posterior
probability is > 65\% can be assigned as first tissue-specific
genetic subtype. The list 'result' contains other 
outputs from eGST. For more details, please see the 'Value' section in eGST manual.

#Getting more details

For any questions, please send an email to statgen.arunabha@gmail.com or pasaniuc@ucla.edu.
See our manuscript for more details: A Majumdar, C Giambartolomei, N Cai, MK Freund, T Haldar, T Schwarz, J Flint, B Pasaniuc (2019) Leveraging eQTLs to identify tissue-specific genetic subtype of complex trait. bioRxiv.
