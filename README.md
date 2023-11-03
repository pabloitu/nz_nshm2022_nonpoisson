# NZ - NSHM2022. Non-Poisson Forecasts and Hazard - Results
Non-Poisson forecasts and hazard input files for the New Zealand 2022 National Seismic Hazard Model.

The repository contains the forecast's files for the Distributed Seismicity Model of the NZNSHM2022,
as well as figures, and the Paraview files to explore them in the software interactively.
{https://www.paraview.org/}. It also contains the Openquake source files (https://github.com/gem/oq-engine) to run a simplified NZ-NSHM2022 model, for single non-Poisson forecasts branches.



## Installation instructions

For reproducibility, this package should install OpenQuake (https://github.com/gem/oq-engine) in its
version v3.16.4. However, Openquake should remain backward compatible for the Negative Binomial formulation.
To install the version 3.16.4, a virtual environment can be created used Anaconda/Miniconda/Micromamba (the latter is recommended, see installation instructions https://mamba.readthedocs.io/en/latest/installation.html) by using:

```shell
conda env create -f environment.yml
```
This environment should already contain the Openquake version.

If the Openquake software should be installed manually into an environment created by the user:

```shell
source activate {user_env}
git clone https://github.com/gem/oq-engine --depth=1 --branch=v3.16.4
cd oq-engine
pip install -e .
```

## Contents
The forecasts from the Hybrid, Poisson Floor and Negative Binomial Floor model types are found at `./forecasts/{model_type}/src` , along with their respective figures in `{model_type}/forecasts/{model_type}/figures`.

The hazard map results for the complete models (all NZ), hazard curves, openquake source files and calculation databases are found in `./hazard_results/maps/{model_type}`, `./hazard_results/curves/{cities}`, `./hazard_results/oq_source/{model_type}` and `./hazard_results/{model_type}/oq_calc`, respectively. For these, a value of N=5.1 and b=0.929 was used to construc the models, with a depth distribution of 15 and 30 km, with equal probability.

An example model construction, calculation and results visualization (for only one site) can be performed by running the file `make_examples.py`

```shell
conda activate nonpoisson_nznshm2022
python make_examples.py
```

The full hazard results can be run in openquake as:

```shell
conda activate nonpoisson_nznshm2022
cd hazard/oq_sources/{model_class}/{model_name}
oq engine --run job.ini
```

## Model types nomenclature

### Hybrid models (Rastin et al., 2022)

* Multiplicative 1346 Grünthal:  m
    * Raw forecast: 'hybrids/forecasts/src/m.csv'
* Additive Total 346 Grünthal:  at
    * Raw forecast: 'hybrids/forecasts/src/at.csv'
* Additive Optimized Grünthal:  ao
    * Raw forecast: 'hybrids/forecasts/src/ao.csv'

### Poisson Floor Ensemble (FE) models

* Floor Multiplicative:  fe_m
    * Raw forecast: 'poisson_floors/forecasts/src/fe_m.csv'
* Floor Additive Total:  fe_at
    * Raw forecast: 'poisson_floors/forecasts/src/fe_at.csv'
* Floor Additive Optimized :  fe_ao
    * Raw forecast: 'poisson_floors/forecasts/src/fe_ao.csv'

### Non-Poisson (Negative Binomial) Floor Ensemble (NPFE) models

* Non-Poisson Floor Multiplicative:  npfe_m
    * Raw forecast: 'nb_floors/forecasts/src/npfe_m.csv'
* Non-Poisson Floor Additive Total:  npfe_at
    * Raw forecast: 'nb_floors/forecasts/src/npfe_at.csv'
* Non-Poisson Floor Additive Optimized :  npfe_ao
    * Raw forecast: 'nb_floors/forecasts/src/npfe_ao.csv'


### EEPAS

The forecasts containing EEPAS, has the eepas suffix attached to them.

## Formats

### Non-Eepas Forecasts
Forecasts are given in CSEP format, where the column `rate` represents the spatial field for N = 1 training earthquake events. Does not add up to 1 in a negative binomial due to non-stationarity.

| lon_min | lon_max | lat_min | lat_max | depth_min | depth_max | m_min | m_max |   rate   | dispersion | mask |
| ------- | ------- | ------- | ------- | --------- | --------- | ----- | ----- | -------- | ---------- | ---- |
| 165.7 | 165.8 | -46.2 | -46.1 | 0.0 | 40.0 | 5.0 | 10.5 | 2.893e-04 | 0. | 1 |
| ------- | ------- | ------- | ------- | --------- | --------- | ----- | ----- | -------- | ---------- | ---- |

### EEPAS forecasts
Forecasts are given in a modified CSEP format, where the column `rate m_i` represents the spatial field for the magnitude m_i, for N = 1 training earthquake events of such magnitude. Does not add up to 1 in a negative binomial due to non-stationarity.

| lon_min | lon_max | lat_min | lat_max | depth_min | depth_max |  rate M=5.0  |   rate M=5.1  | ... |   rate M=8.0  | dispersion |
| ------- | ------- | ------- | ------- | --------- | --------- | ----- | ----- | -------- | ---------- | ---- |
| 165.7 | 165.8 | -46.2 | -46.1 | 0.0 | 40.0 | 3.152e-04 | 3.191e-04 | ... | 3.502e-04 | 0.0 |
| ------- | ------- | ------- | ------- | --------- | --------- | ----- | ----- | -------- | ---------- | ---- |
