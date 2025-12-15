# Trans-Ancestry GWAS: Aortic Stenosis

Genome-wide association study analyzing aortic stenosis genetics in South Asian populations, identifying novel population-specific genetic effects.

## Getting Started

These instructions will help you replicate the statistical genetics analysis.

### Prerequisites

Before running this project, you need:

* R 4.3.0 or higher
* PLINK v1.9
* Required R packages: data.table, ggplot2, dplyr, tidyr
* Command line access

### Installation

Install required R packages:

```
$ R
> install.packages(c("data.table", "ggplot2", "dplyr", "tidyr"))
> quit()
```

Install PLINK:
```
$ wget https://www.cog-genomics.org/plink/1.9/plink_linux_x86_64.zip
$ unzip plink_linux_x86_64.zip
```

## Usage

###Step 1- Heritability Analysis

```
$ Rscript statistical_genetics_AS_pipeline.R
```

###Step 2- Quality Control

```
$ bash plink_qc_commands.sh
```

###Step 3- Association Testing

```
$ bash plink_association_commands.sh
```

###Step 4- Results Analysis

```
$ Rscript association_analysis_pipeline.R
```

##Key Findings

* Novel variant rs11166276 shows opposite effect in South Asians (OR=0.83) vs Europeans (OR=1.07)
* Linked to serum phosphate levels (p=5×10⁻³⁷)
* Demonstrates need for ancestry-specific genetic risk scores

##Technical Details

* Sample size: 7,000 individuals (1,100 cases, 5,900 controls)
* Population: South Asian (Indian/Pakistani)
* SNPs tested: 15 variants from European discovery cohort
* Statistical threshold: Bonferroni-corrected p<0.00333

## Results Structure

```
results/
├── heritability_estimates.txt
├── association_results.assoc.logistic
├── power_calculations.txt
└── manhattan_plot.png
```

## Additional Information

* Analysis framework: Trans-ancestry replication study
* Discovery cohort: European (13,765 cases)
* Method: Case-control logistic regression
