# eDNAssay: a learned model of qPCR cross-amplification
We used supervised machine learning to enhance the prediction of qPCR assay specificity. Our training data were produced via two reaction chemistries: 
SYBR Green intercalating dye and TaqMan MGB probes. Separate models were trained for each and the full-assay model (TaqMan probe-based results) is also available online as [eDNAssay](https://nationalgenomicscenter.shinyapps.io/eDNAssay/). Degenerate bases with IUPAC ambiguity codes are accepted. Indels are treated as N (i.e., any base) as a conservative estimate of assay specificity.

The impetus for building these models was to streamline development of environmental DNA (eDNA) assays. Environmental DNA assays need to discriminate among 
suites of sequences that may be very similar. To ensure assay specificity, eDNA practitioners typically evaluate sequences from all closely related taxa 
(e.g., confamilials) within a pre-defined geographic area. Any taxa that are not deemed "different enough" in computer-based *in silico* testing must be 
put through time- and resource-intensive, laboratory-based *in vitro* testing. However, the determination that an assay is "different enough" *in silico* 
is often dubious. Instead of relying on thermodynamic models and simple mismatch heuristics (as do the vast majority of existing *in silico* tools), our 
models have been trained on empirical data and are therefore highly accurate. Results from model training scripts in this repository will vary slightly 
from those in Kronenberger et al. (2022) depending on the seed selected.

IMPORTANT: eDNAssay was trained using qPCR results generated at the National Genomics Center for Wildlife and Fish Conservation and using it to predict specificity under different reaction conditions may be less accurate. For optimal performance, we recommend users either follow the reaction conditions outlined in Kronenberger et al. (2022) or reassess/retrain the model using their conditions of choice.

## File guide
- **SYBR_training_data.csv** - Empirical dataset containing information on base-pair mismatches, oligonucleotide characteristics, and the results of SYBR Green-based qPCR tests. These data were used to test the primer-only model. See Training_variable_definitions.xlsx for more information.
- **SYBR_testing_data.csv** - Empirical dataset containing SYBR model predictions and qPCR results of assay-template combinations used for testing. See Testing_variable_definitions.xlsx for more information.
- **SYBR_model_training.R** - Script used to train a random forest model to predict cross-amplification of SYBR Green-based qPCR assays.
- **SYBR_trained_model.RData** - The learned primer-only model (SYBR Green results; produced using the SYBR_model_training.R script).
- **SYBR_optimal_thresholds.R** - Script used to calculate optimal class assignment probability thresholds for the primer-only model (SYBR Green results) and a range of false negative (FN) to false positive (FP) cost ratios. For a given FN:FP cost ratio, the threshold that results in the lowest total error cost is optimal.
- **SYBR_specificity_prediction** - Script used to calculate base-pair mismatches between assay oligonucleotides and templates, and then assign templates
probabilities of belonging to the "amplify" class via the learned primer-only model (SYBR Green results).
- **TaqMan_training_data.csv** - Empirical dataset containing information on base-pair mismatches, oligonucleotide characteristics, and the results of TaqMan probe-based qPCR tests. These data were used to train the full-assay model. See Training_variable_definitions.xlsx for more information.
- **TaqMan_testing_data.csv** - Empirical dataset containing SYBR model predictions and qPCR results of assay-template combinations used for testing. See Testing_variable_definitions.xlsx for more information.
- **TaqMan_model_training.R** - Script used to train a random forest model to predict cross-amplification of TaqMan probe-based qPCR assays.
- **TaqMan_trained_model.RData** - The learned full-assay model (TaqMan probe-based results; produced using the TaqMan_model_training.R script), referred to as eDNAssay.
- **TaqMan_optimal_thresholds.R** - Script used to calculate optimal class assignment probability thresholds for the full-assay model (TaqMan probe-based results) and a range of false negative (FN) to false positive (FP) cost ratios. For a given FN:FP cost ratio, the threshold that results in the lowest total error cost is optimal.
- **eDNAssay_offline_version.R** - Script used to calculate base-pair mismatches between assay oligonucleotides and templates, and then assign templates
probabilities of belonging to the "amplify" class via the learned full-assay model (TaqMan probe-based results). This may be used as an alternative to the eDNAssay Shiny app.
- **Testing_variable_definitions.xlsx** - Spreadsheet defining the column names in SYBR_testing_data.csv and TaqMan_testing_data.csv.
- **Training_variable_definitions.xlsx** - Spreadsheet defining the column names in SYBR_training_data.csv and TaqMan_training_data.csv.
- **app.R** - Script behind the eDNAssay Shiny app.
- **eDNAssay_alignment_example.fas** - An example sequence alignment file for use with eDNAssay.
- **eDNAssay_AP_stats.R** - Script used to calculate summary statistics (minimum, maximum, mean, and standard deviation of the mean) 
of assignment probabilties when multiple sequences are included per taxon.

## Contact information
Please reach out to us at the [National Genomics Center for Wildlife and Fish Conservation](https://www.fs.usda.gov/rmrs/ngc) with any questions or comments. 
Scripts and models were created by John Kronenberger at john.kronenberger@usda.gov and Taylor Wilcox at taylor.wilcox@usda.gov.
