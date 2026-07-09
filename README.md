
# exposome - version 0.90

## Overview

*exposome* package provides an interface to fit univarate regression
models for exposome-wide and metabolite-wide association studies.

## How to use exposome

``` r
library(exposome)
```

### Alternatively install from tar.gz:
```r
install.packages(
  "build/ExpoNet_0.0.0.9000.tar.gz",
  repos = NULL,
  type = "source"
)
```

As a simple example, we use linear regression with a binary outcome
(case versus control) on the left side of the equation and define 5
exposures to run multiple regression models on, on the right side of the
equation:

``` r
data(exposome_example)
  data(pheno_example)
  pheno_test <- pheno_example
  pheno_test$Disease_binary <- ifelse(pheno_test$Disease == "Yes", 1, 0)

  fit <- exwas(
    Disease_binary ~ omic_features + age + sex,
    pheno = pheno_test,
    omic_features = exposome_example,
    id_col = "exp_id",
    family = "binomial"
  )

  head(fit$results)
  
# Edit above or add other functions here.
```

## Output

## Plot

## Citing of exposome

## FAQ
