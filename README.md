## eDNAssay: a learned model of qPCR cross-amplification

Our machine learning classifiers predict qPCR cross-amplification via to two reaction chemistries: SYBR Green (primers only) and TaqMan (primers and probe). 
Separate models were trained for each and the TaqMan model is also available online as eDNAssay (https://nationalgenomicscenter.shinyapps.io/eDNAssay/). 

The impetus for building these models was to streamline development of environmental DNA (eDNA) assays. Environmental DNA assays need to discriminate among 
suites of sequences that may very similar. To ensure assay specificity, eDNA practitioners typically evaluate sequences from all closely related taxa 
(e.g., confamilials) within a pre-defined geographic area. Any taxa that are not deemed "different enough" in computer-based in silico testing must be 
put through time- and resource-intensive, laboratory-based in vitro testing. However, the determination that an assay is "different enough" in silico 
is often dubious. Instead of relying on thermodynamic models and simple mismatch heuristics (as do the vast majority of existing in silico tools) our 
models have been trained on empirical data, and are exceptionally accurate as a result.

For optimal performance, we recommend users either 1) develop assays under the reaction conditions used to train these models or 2) test the accuracy 
of models under other reaction conditions before relying on them to declare specificity. See Kronenberger et al. (2022) for details.
