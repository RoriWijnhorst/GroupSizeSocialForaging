# Foraging group size affects repeatable individual variation and strategic social plasticity

This repository contains the data-processing, statistical modelling, and analysis code for the study:

> **Foraging group size affects repeatable individual variation and strategic social plasticity**

The study investigates how social group size affects individual variation in foraging behaviour and social responsiveness in a producer–scrounger game using wild house sparrows (*Passer domesticus*).

## Repository structure

```text
GroupSizeSocialForaging/
│
├── Data/
│   ├── Raw_data/
│   └── Wrangled_data/
│       └── STAN_input_data.rds
│
├── Figures/
│
├── Model_output/
│
├── Models/
│   ├── Mod1.stan
│   └── Mod2.stan
│
├── Scripts/
│   ├── Data_wrangling/
│   └── Analyses/
│   └── Figures/
│   └── Tables/
│
└── README.md
```

### Data

`Data/Raw_data/` contains the raw data used for the analyses, including RFID derived behavioural data and additional assay meta-data. 

`Data/Wrangled_data/` contains processed data prepared for model fitting. The file `STAN_input_data.rds` contains the data object required by the Stan models.

### Models

The Bayesian statistical models are implemented in Stan:

* `Mod1.stan`
* `Mod2.stan`

### Scripts

`Scripts/Data_wrangling/` contains scripts used to process the raw data and prepare the data for analysis.

`Scripts/Analyses/` contains the scripts used to fit the models and generate the analyses, figures, and tables presented in the manuscript.

### Figures

`Figures/` contains figures generated from the analysis scripts.

### Model output

Model output is not included in the repository because of file size. The models can be fitted locally using the provided Stan code and analysis scripts.

## Reproducibility

The analysis workflow is:

```text
Data wrangling
   ↓
Analyses (stan model fitting)
   ↓
Figures and tables
```

To reproduce the analyses, first fit the Stan models using the scripts in `Scripts/Analyses/`. Model results are saved in `Model_output/` and subsequently used to generate the figures and tables.
