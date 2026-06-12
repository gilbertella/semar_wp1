## About the project
This R project supports the analysis of sewage-based antimicrobial resistance surveillance data generated using two complementary approaches: an isolate-based approach and a gene-based metagenomic approach. The overall aim of the project is to directly compare, combine, and evaluate these approaches to determine how well they reflect clinical antimicrobial resistance rates in Escherichia coli.

The study uses municipal sewage samples collected across ten European countries. From these samples, antimicrobial resistance was assessed through two strategies. 
First, an isolate-based approach was applied using susceptibility testing of collected E. coli isolates. 
Second, a gene-based approach was applied using metagenomic sequencing to quantify antimicrobial resistance genes in sewage.

## Data Sources

This project uses data from three main sources:

#### Isolate-based sewage surveillance data
Antimicrobial susceptibility testing results from E. coli isolates recovered from municipal sewage samples.
#### Gene-based sewage surveillance data
Metagenomic sequencing outputs describing the abundance and distribution of antimicrobial resistance genes in sewage samples.
#### Clinical antimicrobial resistance data
Country-level oclinical resistance prevalence estimates for E. coli, covering aminopenicillins, fluoroquinolones, third-generation cephalosporins, and aminoglycosides.

## Analysis Approach

The analysis is implemented in R and focuses on data cleaning, integration, statistical modelling, and visualization. The core modelling framework uses beta regression, which is appropriate for proportional outcomes such as resistance prevalence.

The workflow includes:  
- Importing and cleaning sewage isolate, metagenomic, and clinical datasets.  
- Aggregating resistance indicators by country, sample, antimicrobial class, or surveillance approach.  
- Matching sewage-derived indicators to corresponding clinical resistance outcomes.
- Fitting beta regression models to quantify associations between sewage-based indicators and clinical resistance prevalence.
- Comparing model performance across isolate-based, gene-based, and combined approaches.
- Generating figures and summary tables for interpretation and reporting.

## Software
R (version 4.3 and above)

## Usage
To rerun this analysis, follow the steps below.  
1. Clone the repo

``` 
git clone https://github.com/gilbertella/semar_wp1.git
```
2. Install the following R packages
``` 
library(ggplot2)
library(pals)
library(knitr)
library(cowplot)
library(betareg)
library(tidyverse)
```
3. Load the ```analysis_file.Rmd``` and run all chunks
4. View outputs in the ```tables``` and ```figures``` folders
