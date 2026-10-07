---

## 1. `TEMPO_ACX`

This directory contains processed CMAQ output, satellite/model comparison products, TEMPO-related variables, and scripts used for classification and bias correction.

### `CMAQ_columns_second_run`

Contains **hourly CMAQ column products** generated from the second CMAQ simulation.

Major variables include:

- Tropospheric NOx columns
- Total HCHO columns
- Total ozone columns

These files are primarily used for satellite comparison, TEMPO-related analyses, and generation of gridded column products.

### `CMAQ_PA_columns_second_run`

Contains hourly vertically integrated or columnar **CMAQ process-analysis / physicochemical tendency diagnostics**.

The files include:

- Tropospheric physicochemical tendencies
- Total-column physicochemical tendencies

These outputs are useful for investigating the chemical and physical processes controlling trace-gas abundance and variability.

### `CMAQ_validation_TROPOMI`

Contains scripts used to validate CMAQ simulations against **TROPOMI satellite observations**.

### `CMAQ_CONC`

Contains the **hourly three-dimensional CMAQ concentration fields**.

These files include modeled trace gases and other CMAQ species on the native CMAQ vertical grid.

They represent the primary 3-D chemical concentration output from CMAQ and are used to generate column products and other derived datasets.

### `scripts_classification`

Contains scripts used to extract variables and features needed for classification.


### `TEMPO_TROPOMI_bias_correction`

Contains scripts used to pair TEMPO and TROPOMI observations and estimate systematic differences between the two satellite datasets.

The scripts calculate both:

- **Multiplicative bias**
- **Additive bias**

These corrections can then be applied when harmonizing TEMPO and TROPOMI observations.

### `TEMPO_V*`

Directories beginning with `TEMPO_V` contain processed **TEMPO Level-3 and model-related variables**.

Depending on the processing experiment/version, these directories contain monthly or approximately biweekly products.

---

## 2. `WRFV4.5`

Contains the **Weather Research and Forecasting (WRF) Model version 4.5** and associated simulation files.

WRF provides the meteorological fields required to drive CMAQ.

There are many subdirectories associated with the WRF source code, compilation, utilities, and configuration.

The most important directory for this project is:

### `WRFV4.5/run/`

Contains the WRF runtime environment and the **hourly WRF output files** used for subsequent CMAQ simulations.

These outputs are subsequently processed through MCIP before being used by CMAQ.

---

## 3. `WPS`

Contains the **WRF Preprocessing System (WPS)**.

WPS prepares the meteorological inputs needed to initialize WRF.

The resulting meteorological preprocessing files are then used by WRF to create the simulation initial and boundary conditions.

---

## 4. `CMAQv5.5`

Contains **CMAQ version 5.5** together with the utilities and workflows required to prepare inputs and run the chemical transport model.

This directory includes tools and scripts for several major components of the CMAQ workflow.

### Initial and Boundary Conditions

Utilities are included for generating CMAQ chemical initial and boundary conditions using **GEOS-CF** atmospheric composition fields.

These products provide the chemical state surrounding and initializing the regional CMAQ domain.

### FINN Emissions

Contains tools/workflows used to process **FINN biomass-burning emissions** for CMAQ.

These emissions account for fires and other biomass-burning sources within or affecting the model domain.

### CMAQ Simulation

Contains the configuration files, executables, and scripts used to run CMAQ.

### MCIP

The **Meteorology-Chemistry Interface Processor (MCIP)** converts WRF meteorological output into the meteorological input format required by CMAQ.

MCIP output is also used in the post-processing workflow to convert the three-dimensional `CMAQ_CONC` fields into **columnar quantities**.

The current MCIP outputs are a subset of the whole thing to reduce the space.

These meteorological fields are needed for calculations involving vertical coordinates, atmospheric layer thicknesses, air density, and conversion from model mixing ratios/concentrations into vertically integrated columns.

---

## 5. `OI-SAT-GMI`

`OI-SAT-GMI` is a toolbox developed to match **satellite observations and atmospheric models** in a consistent manner.

Its main purpose is to enable a fair **apples-to-apples comparison** between satellite retrievals and model simulations.

This toolbox was used as part of the **TEMPO Level-3 generation workflow**.

### TEMPO Run Directory

The TEMPO-specific configuration and processing workflow can be found under:

`OI-SAT-GMI/run_CMAQ_TEMPO/`

This directory contains the scripts and configuration required to process CMAQ simulations together with TEMPO observations and generate corresponding Level-3 products.

---

## 6. `CMAQ_scripts`

This directory contains a large collection of miscellaneous scripts accumulated during the WRF-CMAQ and TEMPO projects.

It is effectively a general-purpose working directory containing scripts for many stages of the modeling and analysis workflow.

Examples include scripts for:

- Creating and modifying emissions
- Processing CMAQ inputs
- Validating CMAQ output
- Comparing CMAQ with satellite observations
- Modifying or troubleshooting model inputs
- Post-processing CMAQ output
- Extracting model variables
- Preparing satellite/model comparisons
- Running diagnostic analyses

The directory also contains configuration files and example setups for:

- WRF
- WPS
- CMAQ
- Associated preprocessing tools

Because this directory contains scripts developed for many different experiments, some scripts may be project-specific, experimental, or superseded by newer versions.

