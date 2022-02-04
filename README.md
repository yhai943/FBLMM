# FBLMM
Bayesian linear mixed model with multiple random effects for the analysis of genomic data obtained from family-based studies. 


A Bayesian linear mixed model with multiple random effects (denoted as FBLMM) designed for family-based studies. The proposed model can capture the predictive effects from both common and rare variants, take the information from family design into consideration, and is robust against disease models. It can efficiently select isolated predictors with large effects and a group of predictors with small-to-large effects.


## Installation
This package can be installed from Github by using `devtools`.
```{r, eval=F}
devtools::install_github("senjoro/BLMM")
```
## Usage
Normalize the data before BLMM analysis
```{r, eval=F}
read_file(gene_file,y,kin.mat=F)
```
Compute allele frequenciesGenotypic
```{r, eval=F}
getAlleleFrequencies(...)
```
Take in genotype, phenotype and analyses the target sequence by using BLMM.
```{r, eval=F}
vb_fit_rarecommon(...)
vb_fit_family(...) #Additional need kinship covariance matrix
```
Bayesian linear mixed model with multiple random effects
```{r, eval=F}
vb_predictive(...)
```
Bayesian linear mixed model with multiple random effects for family data
```{r, eval=F}
vb_predictive_family(...)
```
