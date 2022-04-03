# eDNAssay: a learned model of qPCR cross-amplification
We used supervised machine learning to enhance the prediction of qPCR assay specificity. Our training data were produced via two reaction chemistries: 
SYBR Green (primers only) and TaqMan (primers and probe). Separate models were learned for each and the TaqMan model is also available online as [eDNAssay](https://nationalgenomicscenter.shinyapps.io/eDNAssay/). Degenerate bases with IUPAC ambiguity codes are accepted. Indels are treated as N (i.e., any
base) as a conservative estimate of assay specificity.

The impetus for building these models was to streamline development of environmental DNA (eDNA) assays. Environmental DNA assays need to discriminate among 
suites of sequences that may very similar. To ensure assay specificity, eDNA practitioners typically evaluate sequences from all closely related taxa 
(e.g., confamilials) within a pre-defined geographic area. Any taxa that are not deemed "different enough" in computer-based *in silico* testing must be 
put through time- and resource-intensive, laboratory-based *in vitro* testing. However, the determination that an assay is "different enough" *in silico* 
is often dubious. Instead of relying on thermodynamic models and simple mismatch heuristics (as do the vast majority of existing *in silico* tools) our 
models have been trained on empirical data, and are exceptionally accurate as a result. Results from the scripts in this repository will vary slightly from
those in Kronenberger et al. (in review) depending on the seed selected.

For optimal performance, we recommend users either 1) develop assays under the reaction conditions used to train these models or 2) test the accuracy 
of models under other reaction conditions before relying on them to declare specificity. See Kronenberger et al. (in review) for details.

## File guide
- **SYBR_training_data.csv** - Empirical dataset containing information on base-pair mismatches, oligonucleotide characteristics, and the results of SYBR Green-based qPCR tests. These data were used to train the SYBR Green (primer-only) model.
- **SYBR_model_training.R** - Script used to train a random forest model to predict cross-amplification of SYBR Green-based qPCR assays.
- **SYBR_trained_model.RData** - A learned model produced through SYBR Green model training (SYBR_model_training.R script).
- **SYBR_specificity_prediction** - Script used to calculate base-pair mismatches between assay oligonucleotides and templates, and then assign templates
probabilities of belonging to either the "amplify" or "non-amplify" class via the learned SYBR Green model (SYBR_trained_model.RData). A metadata input file is not necessary for this script to run.
- **TaqMan_training_data.csv** - Empirical dataset containing information on base-pair mismatches, oligonucleotide characteristics, and the results of TaqMan-based qPCR tests. These data were used to train the TaqMan (full-assay) model.
- **TaqMan_model_training.R** - Script used to train a random forest model to predict cross-amplification of TaqMan-based qPCR assays.
- **TaqMan_trained_model.RData** - A learned model produced through TaqMan model training (TaqMan_model_training.R script), referred to as eDNAssay.
- **TaqMan_optimal_thresholds.R** - Script used to calculate optimal class assignment probability thresholds for a range of false negative (FN)
to false positive (FP) cost ratios. For a given FN:FP cost ratio, the threshold that results in the lowest total error cost is optimal.
- **TaqMan_optimal_thresholds.RData** - A vector of optimal thresholds for a range of FN:FP cost ratios, produced by the TaqMan_optimal_thresholds.R script.
- **eDNAssay_offline_version.R** - Script used to calculate base-pair mismatches between assay oligonucleotides and templates, and then assign templates
probabilities of belonging to either the "amplify" or "non-amplify" class via the learned TaqMan model (TaqMan_trained_model.RData). This may be used as an 
alternative to the eDNAssay Shiny app. A metadata input file is not necessary for this script to run.
- **app.R** - Script behind the eDNAssay Shiny app.
- **eDNAssay_alignment_example.fas** - An example sequence alignment file for use with the eDNAssay script and app.
- **eDNAssay_metadata_example.csv** - An example metadata file for use with the eDNAssay script and app.
- **eDNAssay_metadata_parse.R** - Script used to parse a .fas file into a .csv file formatted for use with eDNAssay. Some post-processing may be desired 
to ensure species names are consistant.
- **eDNAssay_assignment_probability_stats.R** - Script used to calculate summary statistics (minimum, maximum, mean, and standard deviation of the mean) 
of assignment probailties when multiple sequences are included per taxon.

## Contact information
Please reach out to us at the [National Genomics Center for Wildlife and Fish Conservation](https://www.fs.usda.gov/rmrs/ngc) with any questions or comments. 
Scripts and models were created by Taylor Wilcox at taylor.wilcox@usda.gov and John Kronenberger at john.kronenberger@usda.gov.

