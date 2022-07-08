# Evaluation and comparison of katdetectr

![license](https://img.shields.io/badge/license-GPL--3-blue.svg) [![GitHub issues](https://img.shields.io/github/issues/ErasmusMC-CCBC/katdetectr.svg)]() ![rversion](https://img.shields.io/badge/R%20version-%3E4.1.0-lightgrey.svg)

# Introduction

This repository contains all the custom scripts used in the evaluation and comparison of [katdetectr]('https://github.com/ErasmusMC-CCBC/katdetectr/') as described in the corresponding [Application Note](https://www.google.com) (under submission).

# Usage

All required files were deposited on [Zenodo](link).
These can directly be downloaded using `zen4R` and be used as input.

```R
# Increase the timeout (due to some large files).
options(timeout=5000)

# Download the required files into the data/ folder (~1GB).
zen4R::download_zenodo(doi = "10.5281/zenodo.6810477", path = 'data/')
```