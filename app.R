##################################################################################################
### eDNAssay: A learned model for evaluating qPCR assay specificity
##################################################################################################
library(shiny)
#setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.8/bioc"))
options(repos = BiocManager::repositories())
library(Biostrings)
library(dplyr)
library(caret)
library(randomForest)
library(DT)
library(shinycssloaders)

### Load trained model
load("eDNAssay_trained_model.RData")

### Used below to format slider
sliderInput2 <-
    function(inputId,
             id,
             label,
             min,
             max,
             value,
             step = NULL,
             from_min,
             to_max) {
        x <- sliderInput(inputId, label, min, max, value, step)
        x$children[[2]]$attribs <- c(x$children[[2]]$attribs,
                                     "data-from-min" = from_min,
                                     "data-to-max" = to_max)
        x
    }

##################################################################################################
### Define user interface
ui <- tagList(
    tags$head(tags$style(
        HTML(
            "body {padding-top: 70px; padding-bottom: 20px; padding-left: 20px; padding-right: 20px;}",
            "#big-heading {display: flex; flex-direction: row; justify-content: start; align: left; align-items: center; padding-left: 15px;
            padding-right: 15px; padding-top: 0px; padding-bottom: 0px; color: #024f94; font-size: 34pt; font-weight: bold; width: 100%;}",
            "h1 {font-size: 18pt; font-weight: bold;}",
            "h2 {font-size: 22pt; font-weight: bold; margin-top: 10px; margin-bottom: 0px;}",
            "p {font-size: 11.5pt;}",
            "a {font-size: 11.5pt;}",
            "code {color: black; background-color: #F5F5F5;}",
            ".progress-bar{background-color: #024f94; border-color: #024f94;}",
            ".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar, js-irs-0 .irs-line
            {background: #024f94; border-top: #024f94; border-bottom: #024f94;}",
            ".js-irs-0 .irs-single {font-size: 9.5pt; font-weight: bold; top: 0px; bottom: -20px;}",
            ".irs-grid-text {font-size: 11.5pt; font-weight: bold;}",
            ".irs-from, .irs-to, .irs-min, .irs-max {visibility: hidden !important;}",
            "label {font-size: 11.5pt; margin-bottom: 6px;}",
            ".form-group, .selectize-control {margin-bottom: -10px;}",
            ".shiny-split-layout>div {overflow: hidden;}"
        )
    )),
    
    navbarPage(
        "eDNAssay",
        position = "fixed-top",
        header = div(
            id = "big-heading",
            "eDNAssay: A learned model of qPCR cross-amplification",
            windowTitle = "eDNAssay"
        ),
        
        ### Predict amplification page
        tabPanel(
            "Predict Amplification",
            hr(),
            p(
                paste(
                    'Welcome to eDNAssay! This tool uses supervised machine learning to predict qPCR assay specificity, particularly
        as applied to environmental samples. Simply input aligned sequences and brief metadata. For each template, a
        binary classification model \U2012 trained on TaqMan qPCR results \U2012 outputs the probability of being assigned
        to the "amplify" class. In addition, the Optimize Threshold page can be used to lower the assignment probability
        threshold from the default of 0.5 to better hedge against false negative errors. See the Learn More page, Kronenberger
        et al. (2022), and our'
                ),
                a(href = "https://github.com/NationalGenomicsCenter/eDNAssay", "GitHub repository", target = "_blank"),
                "for more information."
            ),
            p(
                "For optimal performance: 1) reaction conditions from Kronenberger et al. (2022) should be used or the model
            retrained using the conditions of choice, 2) melting temperatures should be ~58-60 C
            for primers and ~68-70 C for probes, and 3) forward primers and probes should anneal to the antisense strand
              and reverse primers to the sense strand. Accuracy has not been tested outside these parameters."
            ),
            br(),
            
            ### Sidebar for inputting and outputting data
            sidebarLayout(
                sidebarPanel(
                    h1("Input aligned sequences"),
                    p(
                        ("\U2022 Input a FASTA file containing aligned sequences."),
                        p(
                            "\U2022 Primer and probe sequences must appear first,
                        ordered as forward primer, reverse primer, then probe."
                        ),
                        p(
                            '\U2022 Name using identical four-letter codes followed
                        by a space and a single-digit oligonucleotide signifier ("XXXX F", "XXXX R", and "XXXX P").'
                        ),
                        p(
                            "\U2022 Only IUPAC-approved characters are allowed. Dashes are treated as Ns (any base) for a
                          conservative estimate of assay specificity."
                        ),
                        p(
                            "\U2022 See"
                            ,
                            a(
                                href = "FVIR_alignment.fas",
                                "example FASTA",
                                download = NA,
                                target = "_blank"
                            )
                        )
                    ),
                    
                    fileInput("alignment", "Aligned sequence file (.fas)"),
                    hr(style = "border-top: 1px solid #A9A9A9;"),
                    h1("Input metadata"),
                    p(
                        (
                            "\U2022 Input a CSV file containing metadata, with rows ordered as in the FASTA file."
                        ),
                        p(
                            '\U2022 There must be three columns: Taxon = taxon name, Name = sequence name (for oligos, name
                            as in the FASTA file) and Type = sequence type ("Oligo" or "Template").'
                        ),
                        p(
                            "\U2022 See "
                            ,
                            a(
                                href = "FVIR_metadata.csv",
                                "example CSV",
                                download = NA,
                                target = "_blank"
                            )
                        )
                    ),
                    fileInput("metadata", "Metadata file (.csv)"),
                    hr(style = "border-top: 1px solid #A9A9A9;"),
                    h1("Input melting temperatures"),
                    numericInput(
                        "F_Tm",
                        HTML(paste0("Forward primer T", tags$sub("m"))),
                        60,
                        min = 45,
                        max = 75,
                        step = 0.1
                    ),
                    br(),
                    numericInput(
                        "R_Tm",
                        HTML(paste0("Reverse primer T", tags$sub("m"))),
                        60,
                        min = 45,
                        max = 75,
                        step = 0.1
                    ),
                    br(),
                    numericInput(
                        "P_Tm",
                        HTML(paste0("Probe T", tags$sub("m"))),
                        70,
                        min = 55,
                        max = 85,
                        step = 0.1
                    ),
                    hr(style = "border-top: 1px solid #A9A9A9;"),
                    h1("Predict amplification"),
                    actionButton("goButton", "PREDICT",
                                 class =
                                     "btn-primary",
                                 style = "margin-top: 20px; background: #024f94"),
                    br(),
                    br(),
                    hr(style = "border-top: 1px solid #A9A9A9;"),
                    
                    h1("Download results"),
                    p(strong("Prediction file (.csv)")),
                    downloadButton("downloadData", "Download"),
                    br(),
                    br(),
                    p(strong("Formatted input file (.csv)")),
                    p(em("For mismatch information, if of interest.")),
                    downloadButton("downloadData2", "Download"),
                    br(),
                ),
                
                # Main panel for displaying outputs ----
                mainPanel(tableOutput("table") %>% withSpinner(color =
                                                                   "#024f94"))
            ),
            br()
        ),
        
        ### Optimize Threshold page
        tabPanel(
            "Optimize Threshold",
            hr(),
            h1(
                "Adjust false negative error cost to determine testing requirements"
            ),
            p(
                paste(
                    'Binary classification models like eDNAssay utilize a threshold value to assign class labels to
            probabilities. The default threshold is typically 0.5. Here, this would mean that templates with
            assignment probabilities < 0.5 are predicted to not amplify and those > 0.5 are predicted to amplify.
            However, this threshold may be moved (or "optimized") to make certain types of errors less likely. This
            is known as a type of cost-sensitive learning, which appreciates that the consequences of false negative
            and false positive errors are often unequal. For an eDNA assay, false negatives (predicting a
            template will not amplify when it does) are typically more "costly" than false positives (predicting a
            template will amplify when it does not). The strength of false negatives relative to false
            positives may change depending on the assay and application. For example, false negative costs may be
            low for a target species that is unprotected and easy to detect using traditional
            methods (e.g., 1X the false positive cost), but high for a target that is protected and highly cryptic
            (e.g., 100X the false positive cost).'
                ),
                em(
                    "Move the slider to see how false negative error tolerance influences
            the optimal assignment threshold (left panel) and the number of taxa requiring in vitro specificity testing (right
            panel; visible following analysis)."
                )
            ),
            br(),
            
            sidebarPanel(
                width = 12,
                height = 12,
                sliderInput2(
                    inputId = "costratio",
                    label = "How many times more impactful is a false negative than a false positive for your assay?",
                    br(),
                    min = 0,
                    max = 100,
                    value = 1,
                    step = 1,
                    from_min = 1,
                    to_max = 100
                )
            ),
            
            mainPanel(splitLayout(
                cellWidths = c("75%", "75%"),
                plotOutput(outputId = "opt"),
                plotOutput(outputId = "hist")
            ))
        ),
        
        ### Information page
        tabPanel(
            "Learn More",
            hr(),
            h1("What this tool is for"),
            p(
                'Quantitative PCR (qPCR) assays applied to environmental samples need to discriminate among suites of
            sequences that may very similar. To ensure assay specificity, environmental DNA (eDNA) practitioners
            typically evaluate sequences from all closely related taxa (e.g., confamilials) within a pre-defined
            geographic area. Any taxa that are not deemed "different enough" in computer-based in silico testing
            must be put through time- and resource-intensive, laboratory-based in vitro testing. However, the
            determination that an assay is "different enough" in silico is often dubious. Practitioners would
            benefit greatly from more accurate and reliable in silico testing methods.'
            ),
            p(
                "To address this need, we developed eDNAssay \U2012 an online tool that uses supervised machine learning to
            predict the probability of qPCR cross-amplification. Instead of relying on thermodynamic models and
            simple mismatch heuristics (as do the vast majority of existing in silico tools) our model uses a random forest
            algorithm to directly incorporate empirical data and achieves exceptional performance as a result. Assignment
            probabilities can be evaluated at different false negative error tolerances to determine which nontarget taxa
            require in vitro testing."
            ),
            hr(),
            
            h1("How to use this tool"),
            p(
                "Users provide three inputs on the Predict Amplification page: 1) a sequence alignment file, 2) a sequence
            metadata file, and 3) estimated oligonucleotide melting temperatures. Input files have specific formatting
            requirements."
            ),
            p(
                "Users receive three outputs: 1) a prediction file with assignment probabilities for assay
            cross-amplification, 2) a data file displaying mismatch count information for off-app analysis \U2012 both
            on the Predict Amplification page, and 3) a suggested optimal threshold based on a user-defined cost
            ratio \U2012 on the Optimize Threshold page. The figure below indicates the relationship between user
            inputs and tool outputs."
            ),
            img(
                src = "eDNAssay_workflow.jpeg",
                height = 450,
                width = 720,
            ),
            p(
                paste(
                    "eDNAssay was trained using qPCR results as detailed in Kronenberger et al. (2022). Using this tool to
                    predict specificity under different reaction conditions will likely produce less accurate results."
                ),
                hr(),
                
                h1("How to cite this tool"),
                p(
                    paste(
                        "Kronenberger JA, Wilcox TM, Mason DH, Franklin TW, McKelvey KS, Young MK, and Schwartz MK (in press).
                  eDNAssay: a machine learning tool that accurately predicts qPCR cross-amplification."
                    ),
                    em("Molecular Ecology Resources.")
                ),
                hr(),
                
                h1("Where to access scripts and previous model versions"),
                p(
                    paste("See our"),
                    a(href = "https://github.com/NationalGenomicsCenter/eDNAssay", "GitHub repository", target = "_blank"),
                    "for the training data and code. This repository also contains a SYBR Green-based primer-only model."
                ),
                hr(),
                
                h1("Published model performance"),
                
                p(
                    "The performance of machine learning classifiers is commonly assessed either through cross-validation and/or using a
          seperate test dataset. Testing the model on new data helps ensure that high accuracy is not a result of overfitting. See these
            metrics below, including confusion matrices and frequency distributions of class assignment probabilities, with incorrect
            predictions in red."
                ),
                
                p(
                    'Cross-validation accuracy was estimated using the results of 268 specificity tests, produced by pairing 10 qPCR assays
            with synthetic gene fragments from 82 nontarget species. "Amplify"" and "nonamplify" classes were unbalanced, so new
            instances of the minority class ("nonamplify") were created for 386 effective tests.'
                ),
                br(),
                img(
                    src = "Training_CM_AP.jpeg",
                    height = 350,
                    width = 700,
                    style = "align: center; padding-left: 0px; padding-right: 0px"
                ),
                br(),
                br(),
                p(
                    "Model performance was also tested via an additional 144 specificity tests, produced by pairing 6 new qPCR assays
              (not used in model training) with genomic DNA from 25 nontarget species. These templates were not sequenced.
              Rather, for each species, we used the mean assigment probability for all sequences on GenBank. Because the actual
              sequences tested may differ from those on GenBank, accuracy estimates may be biased."
                ),
                br(),
                img(
                    src = "Testing_CM_AP.jpeg",
                    height = 350,
                    width = 700,
                    style = "align: center; padding-left: 0px; padding-right: 0px"
                ),
                hr(),
                h1("Unpublished model performance"),
                p("UPDATED JULY 2022"),
                p(
                    "Additional tests of model accuracy have been generated since the publication of Kronenberger et al. (2022) during
              routine specificity testing. These include 171 specificity tests, produced by pairing seven assays with 74 synthetic
              gene fragments. We are only including templates with known sequences as they will be incorporated into the larger
              training dataset for subsequent versions of eDNAssay. Some of these assay oligonucleotides include a degenerate base
              and/or have melting temperatures outside the range of the current training data, indicating that the model may be
              robust to novel assay parameters."
                ),
                br(),
                img(
                    src = "Training_2.0_CM_AP.jpeg",
                    height = 350,
                    width = 700,
                    style = "align: center; padding-left: 0px; padding-right: 0px"
                ),
                hr(),
                
                h1("Questions or comments?"),
                p(
                    paste("Feel free to contact us at the"),
                    a(
                        href = "https://www.fs.usda.gov/rmrs/ngc",
                        "National Genomics Center for Wildlife and Fish Conservation",
                        target = "_blank"
                    ),
                    "with any feedback."
                ),
                p(
                    "eDNAssay was created by John Kronenberger at ",
                    a(
                        href = "john.kronenberger@usda.gov",
                        "john.kronenberger@usda.gov",
                        target = "_blank"
                    ),
                    "and Taylor Wilcox at ",
                    a(href = "taylor.wilcox@usda.gov",
                      "taylor.wilcox@usda.gov",
                      target = "_blank"),
                    br()
                )
            )
        )
    )
)

##################################################################################################
### Define server logic
server <- function(input, output) {
    ### Reactive inputs
    F_Tm <- eventReactive(input$goButton, {
        input$F_Tm
    })
    R_Tm <- eventReactive(input$goButton, {
        input$R_Tm
    })
    P_Tm <- eventReactive(input$goButton, {
        input$P_Tm
    })
    
    threshold <- reactive(input$threshold)
    
    ### Combine matrices and save as a dataframe
    input_metadata <-
        eventReactive(input$goButton, {
            as.matrix(read.csv(input$metadata$datapath))
        })
    input_seqs <-
        eventReactive(input$goButton, {
            as.matrix(readDNAStringSet(input$alignment$datapath))
        })
    
    out_matrix <- eventReactive(input$goButton, {
        req(input$metadata)
        req(input$alignment)
        
        input_matrix <- cbind(input_metadata(), input_seqs())
        input_matrix <-
            as.data.frame(input_matrix, stringsAsFactors = FALSE)
        input_matrix <- na_if(input_matrix, "-")
        
        ### Convert "NA" in template sequences to "N" to accommodate indels
        input_matrix_oligos <- input_matrix[1:3,]
        input_matrix_templates <-
            input_matrix[4:nrow(input_matrix),]
        input_matrix_templates[is.na(input_matrix_templates)] <- "N"
        input_matrix <-
            rbind(input_matrix_oligos, input_matrix_templates)
        
        ### Define length variables and relabel rows
        length_type <- length(input_matrix$Type)
        length_oligo <- length(which(input_matrix$Type == "Oligo"))
        length_template <-
            length(which(input_matrix$Type == "Template"))
        row.names(input_matrix) <- c(1:length_type)
        
        ### Create separate dataframes for each oligo and remove nucleotides outside oligo binding sites
        oligo_matrices <-
            lapply(1:length_oligo, function(x)
                input_matrix[input_matrix$Type == "Template",-which(is.na(input_matrix[x, ]))])
        names(oligo_matrices) = input_matrix$Name[1:length_oligo]
        
        ### Match oligo and template sequences (accounting for IUPAC ambiguity codes) and count the number
        ### of mismatches (=TRUE) in different positions; total, the last five bp on 3' and 5' ends,
        ### and the last bp on the 3' end (terminus)
        ot_match <-
            list() # Returns list of dataframes with sequences of oligo and templates at binding sites
        mm_list <-
            list() # Returns a list of lists with mismatches (accounting for IUPAC ambiguities)
        mm_unlist <-
            list() # Converts mm_list to list of dataframes and finds mismatches in different positions
        
        for (i in 1:length_oligo) {
            ot_match[[i]] <-
                rbind(input_matrix[i, -which(is.na(input_matrix[i, ]))],
                      oligo_matrices[[which(names(oligo_matrices) ==
                                                input_matrix[i, 2])]])
            mm_list[[i]] <-
                lapply(2:(length_template + 1), function(x)
                    ifelse(
                        ot_match[[i]][1, -(1:3)] == ot_match[[i]][x, -(1:3)] |
                            ot_match[[i]][1, -(1:3)] ==
                            "R" & ot_match[[i]][x, -(1:3)] == "A" |
                            ot_match[[i]][1, -(1:3)] ==
                            "R" & ot_match[[i]][x, -(1:3)] == "G" |
                            ot_match[[i]][1, -(1:3)] ==
                            "Y" & ot_match[[i]][x, -(1:3)] == "C" |
                            ot_match[[i]][1, -(1:3)] ==
                            "Y" & ot_match[[i]][x, -(1:3)] == "T" |
                            ot_match[[i]][1, -(1:3)] ==
                            "M" & ot_match[[i]][x, -(1:3)] == "A" |
                            ot_match[[i]][1, -(1:3)] ==
                            "M" & ot_match[[i]][x, -(1:3)] == "C" |
                            ot_match[[i]][1, -(1:3)] ==
                            "K" & ot_match[[i]][x, -(1:3)] == "G" |
                            ot_match[[i]][1, -(1:3)] ==
                            "K" & ot_match[[i]][x, -(1:3)] == "T" |
                            ot_match[[i]][1, -(1:3)] ==
                            "S" & ot_match[[i]][x, -(1:3)] == "G" |
                            ot_match[[i]][1, -(1:3)] ==
                            "S" & ot_match[[i]][x, -(1:3)] == "C" |
                            ot_match[[i]][1, -(1:3)] ==
                            "W" & ot_match[[i]][x, -(1:3)] == "A" |
                            ot_match[[i]][1, -(1:3)] ==
                            "W" & ot_match[[i]][x, -(1:3)] == "T" |
                            ot_match[[i]][1, -(1:3)] ==
                            "B" & ot_match[[i]][x, -(1:3)] == "C" |
                            ot_match[[i]][1, -(1:3)] ==
                            "B" & ot_match[[i]][x, -(1:3)] == "G" |
                            ot_match[[i]][1, -(1:3)] ==
                            "B" & ot_match[[i]][x, -(1:3)] == "T" |
                            ot_match[[i]][1, -(1:3)] ==
                            "D" & ot_match[[i]][x, -(1:3)] == "A" |
                            ot_match[[i]][1, -(1:3)] ==
                            "D" & ot_match[[i]][x, -(1:3)] == "G" |
                            ot_match[[i]][1, -(1:3)] ==
                            "D" & ot_match[[i]][x, -(1:3)] == "T" |
                            ot_match[[i]][1, -(1:3)] ==
                            "H" & ot_match[[i]][x, -(1:3)] == "A" |
                            ot_match[[i]][1, -(1:3)] ==
                            "H" & ot_match[[i]][x, -(1:3)] == "C" |
                            ot_match[[i]][1, -(1:3)] ==
                            "H" & ot_match[[i]][x, -(1:3)] == "T" |
                            ot_match[[i]][1, -(1:3)] ==
                            "V" & ot_match[[i]][x, -(1:3)] == "A" |
                            ot_match[[i]][1, -(1:3)] ==
                            "V" & ot_match[[i]][x, -(1:3)] == "C" |
                            ot_match[[i]][1, -(1:3)] ==
                            "V" & ot_match[[i]][x, -(1:3)] == "G" |
                            ot_match[[i]][1, -(1:3)] ==
                            "N" & ot_match[[i]][x, -(1:3)] == "A" |
                            ot_match[[i]][1, -(1:3)] ==
                            "N" & ot_match[[i]][x, -(1:3)] == "C" |
                            ot_match[[i]][1, -(1:3)] ==
                            "N" & ot_match[[i]][x, -(1:3)] == "G" |
                            ot_match[[i]][1, -(1:3)] ==
                            "N" & ot_match[[i]][x, -(1:3)] == "T" |
                            ot_match[[i]][1, -(1:3)] ==
                            "A" & ot_match[[i]][x, -(1:3)] == "R" |
                            ot_match[[i]][1, -(1:3)] ==
                            "G" & ot_match[[i]][x, -(1:3)] == "R" |
                            ot_match[[i]][1, -(1:3)] ==
                            "C" & ot_match[[i]][x, -(1:3)] == "Y" |
                            ot_match[[i]][1, -(1:3)] ==
                            "T" & ot_match[[i]][x, -(1:3)] == "Y" |
                            ot_match[[i]][1, -(1:3)] ==
                            "A" & ot_match[[i]][x, -(1:3)] == "M" |
                            ot_match[[i]][1, -(1:3)] ==
                            "C" & ot_match[[i]][x, -(1:3)] == "M" |
                            ot_match[[i]][1, -(1:3)] ==
                            "G" & ot_match[[i]][x, -(1:3)] == "K" |
                            ot_match[[i]][1, -(1:3)] ==
                            "T" & ot_match[[i]][x, -(1:3)] == "K" |
                            ot_match[[i]][1, -(1:3)] ==
                            "G" & ot_match[[i]][x, -(1:3)] == "S" |
                            ot_match[[i]][1, -(1:3)] ==
                            "C" & ot_match[[i]][x, -(1:3)] == "S" |
                            ot_match[[i]][1, -(1:3)] ==
                            "A" & ot_match[[i]][x, -(1:3)] == "W" |
                            ot_match[[i]][1, -(1:3)] ==
                            "T" & ot_match[[i]][x, -(1:3)] == "W" |
                            ot_match[[i]][1, -(1:3)] ==
                            "C" & ot_match[[i]][x, -(1:3)] == "B" |
                            ot_match[[i]][1, -(1:3)] ==
                            "G" & ot_match[[i]][x, -(1:3)] == "B" |
                            ot_match[[i]][1, -(1:3)] ==
                            "T" & ot_match[[i]][x, -(1:3)] == "B" |
                            ot_match[[i]][1, -(1:3)] ==
                            "A" & ot_match[[i]][x, -(1:3)] == "D" |
                            ot_match[[i]][1, -(1:3)] ==
                            "G" & ot_match[[i]][x, -(1:3)] == "D" |
                            ot_match[[i]][1, -(1:3)] ==
                            "T" & ot_match[[i]][x, -(1:3)] == "D" |
                            ot_match[[i]][1, -(1:3)] ==
                            "A" & ot_match[[i]][x, -(1:3)] == "H" |
                            ot_match[[i]][1, -(1:3)] ==
                            "C" & ot_match[[i]][x, -(1:3)] == "H" |
                            ot_match[[i]][1, -(1:3)] ==
                            "T" & ot_match[[i]][x, -(1:3)] == "H" |
                            ot_match[[i]][1, -(1:3)] ==
                            "A" & ot_match[[i]][x, -(1:3)] == "V" |
                            ot_match[[i]][1, -(1:3)] ==
                            "C" & ot_match[[i]][x, -(1:3)] == "V" |
                            ot_match[[i]][1, -(1:3)] ==
                            "G" & ot_match[[i]][x, -(1:3)] == "V" |
                            ot_match[[i]][1, -(1:3)] ==
                            "A" & ot_match[[i]][x, -(1:3)] == "N" |
                            ot_match[[i]][1, -(1:3)] ==
                            "C" & ot_match[[i]][x, -(1:3)] == "N" |
                            ot_match[[i]][1, -(1:3)] ==
                            "G" & ot_match[[i]][x, -(1:3)] == "N" |
                            ot_match[[i]][1, -(1:3)] ==
                            "T" & ot_match[[i]][x, -(1:3)] == "N",
                        TRUE,
                        FALSE
                    ))
            mm_unlist[[i]] <-
                data.frame(matrix(
                    unlist(mm_list[[i]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            mm_unlist[c(FALSE, TRUE, FALSE)] = lapply(mm_unlist[c(FALSE, TRUE, FALSE)], function(x)
                rev(x))
            mm_unlist[[i]]$mm_term <-
                as.numeric(mm_unlist[[i]][, ncol(mm_unlist[[i]])] == "FALSE")
            mm_unlist[[i]]$mm_total <-
                rowSums(mm_unlist[[i]][, -ncol(mm_unlist[[i]])] == "FALSE", na.rm = TRUE)
            mm_unlist[[i]]$mm_end3p <-
                rowSums(mm_unlist[[i]][, ((ncol(mm_unlist[[i]]) - 2 - 4):(ncol(mm_unlist[[i]]) -
                                                                              2))] == "FALSE")
            mm_unlist[[i]]$mm_end5p <-
                rowSums(mm_unlist[[i]][, 1:5] == "FALSE")
        }
        
        names(mm_unlist) = input_matrix$Name[1:length_oligo]
        
        ### Create a dataframe for total mismatch number and position
        mm_final <-
            as.data.frame(
                cbind(
                    rep(sub(" .*", "", names(mm_unlist)), each = length_template),
                    rep(input_matrix[1:length_oligo, 2], each =
                            length_template),
                    rep(input_matrix[(length_oligo + 1):length_type, 1], length_oligo),
                    rep(input_matrix[(length_oligo + 1):length_type, 2], length_oligo),
                    rep(input_matrix[(length_oligo + 1):length_type, 3], length_oligo),
                    as.numeric(sapply(mm_unlist, `[[`, "mm_total")),
                    as.numeric(sapply(mm_unlist, `[[`, "mm_end5p")),
                    as.numeric(sapply(mm_unlist, `[[`, "mm_end3p")),
                    as.numeric(sapply(mm_unlist, `[[`, "mm_term"))
                )
            )
        
        colnames(mm_final) <-
            c(
                "Assay",
                "Oligo",
                "Taxon",
                "Name",
                "Type",
                "Total_mm",
                "End5p_mm",
                "End3p_mm",
                "Term_mm"
            )
        mm_final$Total_mm <-
            as.numeric(as.character(mm_final$Total_mm))
        mm_final$End5p_mm <-
            as.numeric(as.character(mm_final$End5p_mm))
        mm_final$End3p_mm <-
            as.numeric(as.character(mm_final$End3p_mm))
        mm_final$Term_mm <-
            as.numeric(as.character(mm_final$Term_mm))
        
        ### Repeat the above loop for all mismatch types (AA, AG/GA, AC/CA, TT, TG/GT, TC/CT, CC, and GG);
        ### loop i used for forward primers (seq(1,length_oligo,2)), j for reverse primers (i+1), and k
        ### for probes (j+1)
        aa_list <-
            list() # Returns a list of lists with AA mismatches
        aa_unlist <-
            list() # Converts aa_list to dataframes and counts AA mismatches in different positions
        
        for (i in seq(1, length_oligo, 3)) {
            for (j in (i + 1)) {
                for (k in (j + 1)) {
                    aa_list[[i]] <- lapply(2:(length_template + 1), function(x)
                        ifelse(
                            ot_match[[i]][1, -(1:3)] == "A" &
                                ot_match[[i]][x, -(1:3)] == "T",
                            TRUE,
                            FALSE
                        ))
                }
                aa_list[[j]] <-
                    lapply(2:(length_template + 1), function(x)
                        ifelse(ot_match[[j]][1, -(1:3)] == "T" &
                                   ot_match[[j]][x, -(1:3)] == "A", TRUE, FALSE))
            }
            aa_list[[k]] <-
                lapply(2:(length_template + 1), function(x)
                    ifelse(ot_match[[k]][1, -(1:3)] == "A" &
                               ot_match[[k]][x, -(1:3)] == "T", TRUE, FALSE))
            aa_unlist[[i]] <-
                data.frame(matrix(
                    unlist(aa_list[[i]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            aa_unlist[[j]] <-
                rev(data.frame(
                    matrix(
                        unlist(aa_list[[j]]),
                        nrow = length_template,
                        byrow = TRUE
                    )
                ))
            aa_unlist[[k]] <-
                data.frame(matrix(
                    unlist(aa_list[[k]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            aa_unlist[[i]]$aa_term <-
                as.numeric(aa_unlist[[i]][, ncol(aa_unlist[[i]])] == "TRUE")
            aa_unlist[[j]]$aa_term <-
                as.numeric(aa_unlist[[j]][, ncol(aa_unlist[[j]])] == "TRUE")
            aa_unlist[[k]]$aa_term <-
                as.numeric(aa_unlist[[k]][, ncol(aa_unlist[[k]])] == "TRUE")
            aa_unlist[[i]]$aa_total <-
                rowSums(aa_unlist[[i]][, -ncol(aa_unlist[[i]])] == "TRUE", na.rm = TRUE)
            aa_unlist[[j]]$aa_total <-
                rowSums(aa_unlist[[j]][, -ncol(aa_unlist[[j]])] == "TRUE", na.rm = TRUE)
            aa_unlist[[k]]$aa_total <-
                rowSums(aa_unlist[[k]][, -ncol(aa_unlist[[k]])] == "TRUE", na.rm = TRUE)
            aa_unlist[[i]]$aa_end3p <-
                rowSums(aa_unlist[[i]][, ((ncol(aa_unlist[[i]]) - 2 - 4):(ncol(aa_unlist[[i]]) -
                                                                              2))] == "TRUE")
            aa_unlist[[j]]$aa_end3p <-
                rowSums(aa_unlist[[j]][, ((ncol(aa_unlist[[j]]) - 2 - 4):(ncol(aa_unlist[[j]]) -
                                                                              2))] == "TRUE")
            aa_unlist[[k]]$aa_end3p <-
                rowSums(aa_unlist[[k]][, ((ncol(aa_unlist[[k]]) - 2 - 4):(ncol(aa_unlist[[k]]) -
                                                                              2))] == "TRUE")
            aa_unlist[[i]]$aa_end5p <-
                rowSums(aa_unlist[[i]][, 1:5] == "TRUE")
            aa_unlist[[j]]$aa_end5p <-
                rowSums(aa_unlist[[j]][, 1:5] == "TRUE")
            aa_unlist[[k]]$aa_end5p <-
                rowSums(aa_unlist[[k]][, 1:5] == "TRUE")
        }
        
        names(aa_unlist) = input_matrix$Name[1:length_oligo]
        
        #-------------------------------------
        ag_list <-
            list() # Returns a list of lists with AG and GA mismatches
        ag_unlist <-
            list() # Converts ag_list to dataframes and counts AG and GA mismatches in different positions
        
        for (i in seq(1, length_oligo, 3)) {
            for (j in (i + 1)) {
                for (k in (j + 1)) {
                    ag_list[[i]] <- lapply(2:(length_template + 1), function(x)
                        ifelse(
                            ot_match[[i]][1, -(1:3)] == "A" & ot_match[[i]][x, -(1:3)] == "C" |
                                ot_match[[i]][1, -(1:3)] == "G" &
                                ot_match[[i]][x, -(1:3)] == "T",
                            TRUE,
                            FALSE
                        ))
                }
                ag_list[[j]] <-
                    lapply(2:(length_template + 1), function(x)
                        ifelse(
                            ot_match[[j]][1, -(1:3)] == "T" & ot_match[[j]][x, -(1:3)] == "G" |
                                ot_match[[j]][1, -(1:3)] == "C" &
                                ot_match[[j]][x, -(1:3)] == "A",
                            TRUE,
                            FALSE
                        ))
            }
            ag_list[[k]] <-
                lapply(2:(length_template + 1), function(x)
                    ifelse(
                        ot_match[[k]][1, -(1:3)] == "A" & ot_match[[k]][x, -(1:3)] == "C" |
                            ot_match[[k]][1, -(1:3)] == "G" &
                            ot_match[[k]][x, -(1:3)] == "T",
                        TRUE,
                        FALSE
                    ))
            ag_unlist[[i]] <-
                data.frame(matrix(
                    unlist(ag_list[[i]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            ag_unlist[[j]] <-
                rev(data.frame(
                    matrix(
                        unlist(ag_list[[j]]),
                        nrow = length_template,
                        byrow = TRUE
                    )
                ))
            ag_unlist[[k]] <-
                data.frame(matrix(
                    unlist(ag_list[[k]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            ag_unlist[[i]]$ag_term <-
                as.numeric(ag_unlist[[i]][, ncol(ag_unlist[[i]])] == "TRUE")
            ag_unlist[[j]]$ag_term <-
                as.numeric(ag_unlist[[j]][, ncol(ag_unlist[[j]])] == "TRUE")
            ag_unlist[[k]]$ag_term <-
                as.numeric(ag_unlist[[k]][, ncol(ag_unlist[[k]])] == "TRUE")
            ag_unlist[[i]]$ag_total <-
                rowSums(ag_unlist[[i]][, -ncol(ag_unlist[[i]])] == "TRUE", na.rm = TRUE)
            ag_unlist[[j]]$ag_total <-
                rowSums(ag_unlist[[j]][, -ncol(ag_unlist[[j]])] == "TRUE", na.rm = TRUE)
            ag_unlist[[k]]$ag_total <-
                rowSums(ag_unlist[[k]][, -ncol(ag_unlist[[k]])] == "TRUE", na.rm = TRUE)
            ag_unlist[[i]]$ag_end3p <-
                rowSums(ag_unlist[[i]][, ((ncol(ag_unlist[[i]]) - 2 - 4):(ncol(ag_unlist[[i]]) -
                                                                              2))] == "TRUE")
            ag_unlist[[j]]$ag_end3p <-
                rowSums(ag_unlist[[j]][, ((ncol(ag_unlist[[j]]) - 2 - 4):(ncol(ag_unlist[[j]]) -
                                                                              2))] == "TRUE")
            ag_unlist[[k]]$ag_end3p <-
                rowSums(ag_unlist[[k]][, ((ncol(ag_unlist[[k]]) - 2 - 4):(ncol(ag_unlist[[k]]) -
                                                                              2))] == "TRUE")
            ag_unlist[[i]]$ag_end5p <-
                rowSums(ag_unlist[[i]][, 1:5] == "TRUE")
            ag_unlist[[j]]$ag_end5p <-
                rowSums(ag_unlist[[j]][, 1:5] == "TRUE")
            ag_unlist[[k]]$ag_end5p <-
                rowSums(ag_unlist[[k]][, 1:5] == "TRUE")
        }
        
        names(ag_unlist) = input_matrix$Name[1:length_oligo]
        
        #-------------------------------------
        ac_list <-
            list() # Returns a list of lists with AC and CA mismatches
        ac_unlist <-
            list() # Converts ac_list to dataframes and counts AC and CA mismatches in different positions
        
        for (i in seq(1, length_oligo, 3)) {
            for (j in (i + 1)) {
                for (k in (j + 1)) {
                    ac_list[[i]] <- lapply(2:(length_template + 1), function(x)
                        ifelse(
                            ot_match[[i]][1, -(1:3)] == "A" & ot_match[[i]][x, -(1:3)] == "G" |
                                ot_match[[i]][1, -(1:3)] == "C" &
                                ot_match[[i]][x, -(1:3)] == "T",
                            TRUE,
                            FALSE
                        ))
                }
                ac_list[[j]] <-
                    lapply(2:(length_template + 1), function(x)
                        ifelse(
                            ot_match[[j]][1, -(1:3)] == "T" & ot_match[[j]][x, -(1:3)] == "C" |
                                ot_match[[j]][1, -(1:3)] == "G" &
                                ot_match[[j]][x, -(1:3)] == "A",
                            TRUE,
                            FALSE
                        ))
            }
            ac_list[[k]] <-
                lapply(2:(length_template + 1), function(x)
                    ifelse(
                        ot_match[[k]][1, -(1:3)] == "A" & ot_match[[k]][x, -(1:3)] == "G" |
                            ot_match[[k]][1, -(1:3)] == "C" &
                            ot_match[[k]][x, -(1:3)] == "T",
                        TRUE,
                        FALSE
                    ))
            ac_unlist[[i]] <-
                data.frame(matrix(
                    unlist(ac_list[[i]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            ac_unlist[[j]] <-
                rev(data.frame(
                    matrix(
                        unlist(ac_list[[j]]),
                        nrow = length_template,
                        byrow = TRUE
                    )
                ))
            ac_unlist[[k]] <-
                data.frame(matrix(
                    unlist(ac_list[[k]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            ac_unlist[[i]]$ac_term <-
                as.numeric(ac_unlist[[i]][, ncol(ac_unlist[[i]])] == "TRUE")
            ac_unlist[[j]]$ac_term <-
                as.numeric(ac_unlist[[j]][, ncol(ac_unlist[[j]])] == "TRUE")
            ac_unlist[[k]]$ac_term <-
                as.numeric(ac_unlist[[k]][, ncol(ac_unlist[[k]])] == "TRUE")
            ac_unlist[[i]]$ac_total <-
                rowSums(ac_unlist[[i]][, -ncol(ac_unlist[[i]])] == "TRUE", na.rm = TRUE)
            ac_unlist[[j]]$ac_total <-
                rowSums(ac_unlist[[j]][, -ncol(ac_unlist[[j]])] == "TRUE", na.rm = TRUE)
            ac_unlist[[k]]$ac_total <-
                rowSums(ac_unlist[[k]][, -ncol(ac_unlist[[k]])] == "TRUE", na.rm = TRUE)
            ac_unlist[[i]]$ac_end3p <-
                rowSums(ac_unlist[[i]][, ((ncol(ac_unlist[[i]]) - 2 - 4):(ncol(ac_unlist[[i]]) -
                                                                              2))] == "TRUE")
            ac_unlist[[j]]$ac_end3p <-
                rowSums(ac_unlist[[j]][, ((ncol(ac_unlist[[j]]) - 2 - 4):(ncol(ac_unlist[[j]]) -
                                                                              2))] == "TRUE")
            ac_unlist[[k]]$ac_end3p <-
                rowSums(ac_unlist[[k]][, ((ncol(ac_unlist[[k]]) - 2 - 4):(ncol(ac_unlist[[k]]) -
                                                                              2))] == "TRUE")
            ac_unlist[[i]]$ac_end5p <-
                rowSums(ac_unlist[[i]][, 1:5] == "TRUE")
            ac_unlist[[j]]$ac_end5p <-
                rowSums(ac_unlist[[j]][, 1:5] == "TRUE")
            ac_unlist[[k]]$ac_end5p <-
                rowSums(ac_unlist[[k]][, 1:5] == "TRUE")
        }
        
        names(ac_unlist) = input_matrix$Name[1:length_oligo]
        
        #-------------------------------------
        tt_list <-
            list() # Returns a list of lists with TT mismatches
        tt_unlist <-
            list() # Converts tt_list to dataframes and counts TT mismatches in different positions
        
        for (i in seq(1, length_oligo, 3)) {
            for (j in (i + 1)) {
                for (k in (j + 1)) {
                    tt_list[[i]] <- lapply(2:(length_template + 1), function(x)
                        ifelse(
                            ot_match[[i]][1, -(1:3)] == "T" &
                                ot_match[[i]][x, -(1:3)] == "A",
                            TRUE,
                            FALSE
                        ))
                }
                tt_list[[j]] <-
                    lapply(2:(length_template + 1), function(x)
                        ifelse(ot_match[[j]][1, -(1:3)] == "A" &
                                   ot_match[[j]][x, -(1:3)] == "T", TRUE, FALSE))
            }
            tt_list[[k]] <-
                lapply(2:(length_template + 1), function(x)
                    ifelse(ot_match[[k]][1, -(1:3)] == "T" &
                               ot_match[[k]][x, -(1:3)] == "A", TRUE, FALSE))
            tt_unlist[[i]] <-
                data.frame(matrix(
                    unlist(tt_list[[i]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            tt_unlist[[j]] <-
                rev(data.frame(
                    matrix(
                        unlist(tt_list[[j]]),
                        nrow = length_template,
                        byrow = TRUE
                    )
                ))
            tt_unlist[[k]] <-
                data.frame(matrix(
                    unlist(tt_list[[k]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            tt_unlist[[i]]$tt_term <-
                as.numeric(tt_unlist[[i]][, ncol(tt_unlist[[i]])] == "TRUE")
            tt_unlist[[j]]$tt_term <-
                as.numeric(tt_unlist[[j]][, ncol(tt_unlist[[j]])] == "TRUE")
            tt_unlist[[k]]$tt_term <-
                as.numeric(tt_unlist[[k]][, ncol(tt_unlist[[k]])] == "TRUE")
            tt_unlist[[i]]$tt_total <-
                rowSums(tt_unlist[[i]][, -ncol(tt_unlist[[i]])] == "TRUE", na.rm = TRUE)
            tt_unlist[[j]]$tt_total <-
                rowSums(tt_unlist[[j]][, -ncol(tt_unlist[[j]])] == "TRUE", na.rm = TRUE)
            tt_unlist[[k]]$tt_total <-
                rowSums(tt_unlist[[k]][, -ncol(tt_unlist[[k]])] == "TRUE", na.rm = TRUE)
            tt_unlist[[i]]$tt_end3p <-
                rowSums(tt_unlist[[i]][, ((ncol(tt_unlist[[i]]) - 2 - 4):(ncol(tt_unlist[[i]]) -
                                                                              2))] == "TRUE")
            tt_unlist[[j]]$tt_end3p <-
                rowSums(tt_unlist[[j]][, ((ncol(tt_unlist[[j]]) - 2 - 4):(ncol(tt_unlist[[j]]) -
                                                                              2))] == "TRUE")
            tt_unlist[[k]]$tt_end3p <-
                rowSums(tt_unlist[[k]][, ((ncol(tt_unlist[[k]]) - 2 - 4):(ncol(tt_unlist[[k]]) -
                                                                              2))] == "TRUE")
            tt_unlist[[i]]$tt_end5p <-
                rowSums(tt_unlist[[i]][, 1:5] == "TRUE")
            tt_unlist[[j]]$tt_end5p <-
                rowSums(tt_unlist[[j]][, 1:5] == "TRUE")
            tt_unlist[[k]]$tt_end5p <-
                rowSums(tt_unlist[[k]][, 1:5] == "TRUE")
        }
        
        names(tt_unlist) = input_matrix$Name[1:length_oligo]
        
        #-------------------------------------
        tg_list <-
            list() # Returns a list of lists with TG and GT mismatches
        tg_unlist <-
            list() # Converts tg_list to dataframes and counts TG and GT mismatches in different positions
        
        for (i in seq(1, length_oligo, 3)) {
            for (j in (i + 1)) {
                for (k in (j + 1)) {
                    tg_list[[i]] <- lapply(2:(length_template + 1), function(x)
                        ifelse(
                            ot_match[[i]][1, -(1:3)] == "T" & ot_match[[i]][x, -(1:3)] == "C" |
                                ot_match[[i]][1, -(1:3)] == "G" &
                                ot_match[[i]][x, -(1:3)] == "A",
                            TRUE,
                            FALSE
                        ))
                }
                tg_list[[j]] <-
                    lapply(2:(length_template + 1), function(x)
                        ifelse(
                            ot_match[[j]][1, -(1:3)] == "A" & ot_match[[j]][x, -(1:3)] == "G" |
                                ot_match[[j]][1, -(1:3)] == "C" &
                                ot_match[[j]][x, -(1:3)] == "T",
                            TRUE,
                            FALSE
                        ))
            }
            tg_list[[k]] <-
                lapply(2:(length_template + 1), function(x)
                    ifelse(
                        ot_match[[k]][1, -(1:3)] == "T" & ot_match[[k]][x, -(1:3)] == "C" |
                            ot_match[[k]][1, -(1:3)] == "G" &
                            ot_match[[k]][x, -(1:3)] == "A",
                        TRUE,
                        FALSE
                    ))
            tg_unlist[[i]] <-
                data.frame(matrix(
                    unlist(tg_list[[i]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            tg_unlist[[j]] <-
                rev(data.frame(
                    matrix(
                        unlist(tg_list[[j]]),
                        nrow = length_template,
                        byrow = TRUE
                    )
                ))
            tg_unlist[[k]] <-
                data.frame(matrix(
                    unlist(tg_list[[k]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            tg_unlist[[i]]$tg_term <-
                as.numeric(tg_unlist[[i]][, ncol(tg_unlist[[i]])] == "TRUE")
            tg_unlist[[j]]$tg_term <-
                as.numeric(tg_unlist[[j]][, ncol(tg_unlist[[j]])] == "TRUE")
            tg_unlist[[k]]$tg_term <-
                as.numeric(tg_unlist[[k]][, ncol(tg_unlist[[k]])] == "TRUE")
            tg_unlist[[i]]$tg_total <-
                rowSums(tg_unlist[[i]][, -ncol(tg_unlist[[i]])] == "TRUE", na.rm = TRUE)
            tg_unlist[[j]]$tg_total <-
                rowSums(tg_unlist[[j]][, -ncol(tg_unlist[[j]])] == "TRUE", na.rm = TRUE)
            tg_unlist[[k]]$tg_total <-
                rowSums(tg_unlist[[k]][, -ncol(tg_unlist[[k]])] == "TRUE", na.rm = TRUE)
            tg_unlist[[i]]$tg_end3p <-
                rowSums(tg_unlist[[i]][, ((ncol(tg_unlist[[i]]) - 2 - 4):(ncol(tg_unlist[[i]]) -
                                                                              2))] == "TRUE")
            tg_unlist[[j]]$tg_end3p <-
                rowSums(tg_unlist[[j]][, ((ncol(tg_unlist[[j]]) - 2 - 4):(ncol(tg_unlist[[j]]) -
                                                                              2))] == "TRUE")
            tg_unlist[[k]]$tg_end3p <-
                rowSums(tg_unlist[[k]][, ((ncol(tg_unlist[[k]]) - 2 - 4):(ncol(tg_unlist[[k]]) -
                                                                              2))] == "TRUE")
            tg_unlist[[i]]$tg_end5p <-
                rowSums(tg_unlist[[i]][, 1:5] == "TRUE")
            tg_unlist[[j]]$tg_end5p <-
                rowSums(tg_unlist[[j]][, 1:5] == "TRUE")
            tg_unlist[[k]]$tg_end5p <-
                rowSums(tg_unlist[[k]][, 1:5] == "TRUE")
        }
        
        names(tg_unlist) = input_matrix$Name[1:length_oligo]
        
        #-------------------------------------
        tc_list <-
            list() # Returns a list of lists with TC and CT mismatches
        tc_unlist <-
            list() # Converts tc_list to dataframes and counts TC and CT mismatches in different positions
        
        for (i in seq(1, length_oligo, 3)) {
            for (j in (i + 1)) {
                for (k in (j + 1)) {
                    tc_list[[i]] <- lapply(2:(length_template + 1), function(x)
                        ifelse(
                            ot_match[[i]][1, -(1:3)] == "T" & ot_match[[i]][x, -(1:3)] == "G" |
                                ot_match[[i]][1, -(1:3)] == "C" &
                                ot_match[[i]][x, -(1:3)] == "A",
                            TRUE,
                            FALSE
                        ))
                }
                tc_list[[j]] <-
                    lapply(2:(length_template + 1), function(x)
                        ifelse(
                            ot_match[[j]][1, -(1:3)] == "A" & ot_match[[j]][x, -(1:3)] == "C" |
                                ot_match[[j]][1, -(1:3)] == "G" &
                                ot_match[[j]][x, -(1:3)] == "T",
                            TRUE,
                            FALSE
                        ))
            }
            tc_list[[k]] <-
                lapply(2:(length_template + 1), function(x)
                    ifelse(
                        ot_match[[k]][1, -(1:3)] == "T" & ot_match[[k]][x, -(1:3)] == "G" |
                            ot_match[[k]][1, -(1:3)] == "C" &
                            ot_match[[k]][x, -(1:3)] == "A",
                        TRUE,
                        FALSE
                    ))
            tc_unlist[[i]] <-
                data.frame(matrix(
                    unlist(tc_list[[i]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            tc_unlist[[j]] <-
                rev(data.frame(
                    matrix(
                        unlist(tc_list[[j]]),
                        nrow = length_template,
                        byrow = TRUE
                    )
                ))
            tc_unlist[[k]] <-
                data.frame(matrix(
                    unlist(tc_list[[k]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            tc_unlist[[i]]$tc_term <-
                as.numeric(tc_unlist[[i]][, ncol(tc_unlist[[i]])] == "TRUE")
            tc_unlist[[j]]$tc_term <-
                as.numeric(tc_unlist[[j]][, ncol(tc_unlist[[j]])] == "TRUE")
            tc_unlist[[k]]$tc_term <-
                as.numeric(tc_unlist[[k]][, ncol(tc_unlist[[k]])] == "TRUE")
            tc_unlist[[i]]$tc_total <-
                rowSums(tc_unlist[[i]][, -ncol(tc_unlist[[i]])] == "TRUE", na.rm = TRUE)
            tc_unlist[[j]]$tc_total <-
                rowSums(tc_unlist[[j]][, -ncol(tc_unlist[[j]])] == "TRUE", na.rm = TRUE)
            tc_unlist[[k]]$tc_total <-
                rowSums(tc_unlist[[k]][, -ncol(tc_unlist[[k]])] == "TRUE", na.rm = TRUE)
            tc_unlist[[i]]$tc_end3p <-
                rowSums(tc_unlist[[i]][, ((ncol(tc_unlist[[i]]) - 2 - 4):(ncol(tc_unlist[[i]]) -
                                                                              2))] == "TRUE")
            tc_unlist[[j]]$tc_end3p <-
                rowSums(tc_unlist[[j]][, ((ncol(tc_unlist[[j]]) - 2 - 4):(ncol(tc_unlist[[j]]) -
                                                                              2))] == "TRUE")
            tc_unlist[[k]]$tc_end3p <-
                rowSums(tc_unlist[[k]][, ((ncol(tc_unlist[[k]]) - 2 - 4):(ncol(tc_unlist[[k]]) -
                                                                              2))] == "TRUE")
            tc_unlist[[i]]$tc_end5p <-
                rowSums(tc_unlist[[i]][, 1:5] == "TRUE")
            tc_unlist[[j]]$tc_end5p <-
                rowSums(tc_unlist[[j]][, 1:5] == "TRUE")
            tc_unlist[[k]]$tc_end5p <-
                rowSums(tc_unlist[[k]][, 1:5] == "TRUE")
        }
        
        names(tc_unlist) = input_matrix$Name[1:length_oligo]
        
        #-------------------------------------
        gg_list <-
            list() # Returns a list of lists with GG mismatches
        gg_unlist <-
            list() # Converts gg_list to dataframes and counts GG mismatches in different positions
        
        for (i in seq(1, length_oligo, 3)) {
            for (j in (i + 1)) {
                for (k in (j + 1)) {
                    gg_list[[i]] <- lapply(2:(length_template + 1), function(x)
                        ifelse(
                            ot_match[[i]][1, -(1:3)] == "G" &
                                ot_match[[i]][x, -(1:3)] == "C",
                            TRUE,
                            FALSE
                        ))
                }
                gg_list[[j]] <-
                    lapply(2:(length_template + 1), function(x)
                        ifelse(ot_match[[j]][1, -(1:3)] == "C" &
                                   ot_match[[j]][x, -(1:3)] == "G", TRUE, FALSE))
            }
            gg_list[[k]] <-
                lapply(2:(length_template + 1), function(x)
                    ifelse(ot_match[[k]][1, -(1:3)] == "G" &
                               ot_match[[k]][x, -(1:3)] == "C", TRUE, FALSE))
            gg_unlist[[i]] <-
                data.frame(matrix(
                    unlist(gg_list[[i]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            gg_unlist[[j]] <-
                rev(data.frame(
                    matrix(
                        unlist(gg_list[[j]]),
                        nrow = length_template,
                        byrow = TRUE
                    )
                ))
            gg_unlist[[k]] <-
                data.frame(matrix(
                    unlist(gg_list[[k]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            gg_unlist[[i]]$gg_term <-
                as.numeric(gg_unlist[[i]][, ncol(gg_unlist[[i]])] == "TRUE")
            gg_unlist[[j]]$gg_term <-
                as.numeric(gg_unlist[[j]][, ncol(gg_unlist[[j]])] == "TRUE")
            gg_unlist[[k]]$gg_term <-
                as.numeric(gg_unlist[[k]][, ncol(gg_unlist[[k]])] == "TRUE")
            gg_unlist[[i]]$gg_total <-
                rowSums(gg_unlist[[i]][, -ncol(gg_unlist[[i]])] == "TRUE", na.rm = TRUE)
            gg_unlist[[j]]$gg_total <-
                rowSums(gg_unlist[[j]][, -ncol(gg_unlist[[j]])] == "TRUE", na.rm = TRUE)
            gg_unlist[[k]]$gg_total <-
                rowSums(gg_unlist[[k]][, -ncol(gg_unlist[[k]])] == "TRUE", na.rm = TRUE)
            gg_unlist[[i]]$gg_end3p <-
                rowSums(gg_unlist[[i]][, ((ncol(gg_unlist[[i]]) - 2 - 4):(ncol(gg_unlist[[i]]) -
                                                                              2))] == "TRUE")
            gg_unlist[[j]]$gg_end3p <-
                rowSums(gg_unlist[[j]][, ((ncol(gg_unlist[[j]]) - 2 - 4):(ncol(gg_unlist[[j]]) -
                                                                              2))] == "TRUE")
            gg_unlist[[k]]$gg_end3p <-
                rowSums(gg_unlist[[k]][, ((ncol(gg_unlist[[k]]) - 2 - 4):(ncol(gg_unlist[[k]]) -
                                                                              2))] == "TRUE")
            gg_unlist[[i]]$gg_end5p <-
                rowSums(gg_unlist[[i]][, 1:5] == "TRUE")
            gg_unlist[[j]]$gg_end5p <-
                rowSums(gg_unlist[[j]][, 1:5] == "TRUE")
            gg_unlist[[k]]$gg_end5p <-
                rowSums(gg_unlist[[k]][, 1:5] == "TRUE")
        }
        
        names(gg_unlist) = input_matrix$Name[1:length_oligo]
        
        #-------------------------------------
        cc_list <-
            list() # Returns a list of lists with CC mismatches
        cc_unlist <-
            list() # Converts cc_list to dataframes and counts CC mismatches in different positions
        
        for (i in seq(1, length_oligo, 3)) {
            for (j in (i + 1)) {
                for (k in (j + 1)) {
                    cc_list[[i]] <- lapply(2:(length_template + 1), function(x)
                        ifelse(
                            ot_match[[i]][1, -(1:3)] == "C" &
                                ot_match[[i]][x, -(1:3)] == "G",
                            TRUE,
                            FALSE
                        ))
                }
                cc_list[[j]] <-
                    lapply(2:(length_template + 1), function(x)
                        ifelse(ot_match[[j]][1, -(1:3)] == "G" &
                                   ot_match[[j]][x, -(1:3)] == "C", TRUE, FALSE))
            }
            cc_list[[k]] <-
                lapply(2:(length_template + 1), function(x)
                    ifelse(ot_match[[k]][1, -(1:3)] == "C" &
                               ot_match[[k]][x, -(1:3)] == "G", TRUE, FALSE))
            cc_unlist[[i]] <-
                data.frame(matrix(
                    unlist(cc_list[[i]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            cc_unlist[[j]] <-
                rev(data.frame(
                    matrix(
                        unlist(cc_list[[j]]),
                        nrow = length_template,
                        byrow = TRUE
                    )
                ))
            cc_unlist[[k]] <-
                data.frame(matrix(
                    unlist(cc_list[[k]]),
                    nrow = length_template,
                    byrow = TRUE
                ))
            cc_unlist[[i]]$cc_term <-
                as.numeric(cc_unlist[[i]][, ncol(cc_unlist[[i]])] == "TRUE")
            cc_unlist[[j]]$cc_term <-
                as.numeric(cc_unlist[[j]][, ncol(cc_unlist[[j]])] == "TRUE")
            cc_unlist[[k]]$cc_term <-
                as.numeric(cc_unlist[[k]][, ncol(cc_unlist[[k]])] == "TRUE")
            cc_unlist[[i]]$cc_total <-
                rowSums(cc_unlist[[i]][, -ncol(cc_unlist[[i]])] == "TRUE", na.rm = TRUE)
            cc_unlist[[j]]$cc_total <-
                rowSums(cc_unlist[[j]][, -ncol(cc_unlist[[j]])] == "TRUE", na.rm = TRUE)
            cc_unlist[[k]]$cc_total <-
                rowSums(cc_unlist[[k]][, -ncol(cc_unlist[[k]])] == "TRUE", na.rm = TRUE)
            cc_unlist[[i]]$cc_end3p <-
                rowSums(cc_unlist[[i]][, ((ncol(cc_unlist[[i]]) - 2 - 4):(ncol(cc_unlist[[i]]) -
                                                                              2))] == "TRUE")
            cc_unlist[[j]]$cc_end3p <-
                rowSums(cc_unlist[[j]][, ((ncol(cc_unlist[[j]]) - 2 - 4):(ncol(cc_unlist[[j]]) -
                                                                              2))] == "TRUE")
            cc_unlist[[k]]$cc_end3p <-
                rowSums(cc_unlist[[k]][, ((ncol(cc_unlist[[k]]) - 2 - 4):(ncol(cc_unlist[[k]]) -
                                                                              2))] == "TRUE")
            cc_unlist[[i]]$cc_end5p <-
                rowSums(cc_unlist[[i]][, 1:5] == "TRUE")
            cc_unlist[[j]]$cc_end5p <-
                rowSums(cc_unlist[[j]][, 1:5] == "TRUE")
            cc_unlist[[k]]$cc_end5p <-
                rowSums(cc_unlist[[k]][, 1:5] == "TRUE")
        }
        
        names(cc_unlist) = input_matrix$Name[1:length_oligo]
        
        ### Create a dataframe for mismatch type data
        mm_type_final <-
            as.data.frame(
                cbind(
                    rep(sub(" .*", "", names(mm_unlist)), each = length_template),
                    rep(input_matrix[1:length_oligo, 2], each =
                            length_template),
                    rep(input_matrix[(length_oligo +
                                          1):length_type, 1], length_oligo),
                    rep(input_matrix[(length_oligo +
                                          1):length_type, 2], length_oligo),
                    rep(input_matrix[(length_oligo +
                                          1):length_type, 3], length_oligo),
                    as.numeric(sapply(aa_unlist, `[[`, "aa_total")),
                    as.numeric(sapply(ag_unlist, `[[`, "ag_total")),
                    as.numeric(sapply(ac_unlist, `[[`, "ac_total")),
                    as.numeric(sapply(tt_unlist, `[[`, "tt_total")),
                    as.numeric(sapply(tg_unlist, `[[`, "tg_total")),
                    as.numeric(sapply(tc_unlist, `[[`, "tc_total")),
                    as.numeric(sapply(gg_unlist, `[[`, "gg_total")),
                    as.numeric(sapply(cc_unlist, `[[`, "cc_total"))
                )
            )
        
        colnames(mm_type_final) <-
            c(
                "Assay",
                "Oligo",
                "Taxon",
                "Name",
                "Type",
                "AA",
                "AG",
                "AC",
                "TT",
                "TG",
                "TC",
                "GG",
                "CC"
            )
        mm_type_final$AA <-
            as.numeric(as.character(mm_type_final$AA))
        mm_type_final$AG <-
            as.numeric(as.character(mm_type_final$AG))
        mm_type_final$AC <-
            as.numeric(as.character(mm_type_final$AC))
        mm_type_final$TT <-
            as.numeric(as.character(mm_type_final$TT))
        mm_type_final$TG <-
            as.numeric(as.character(mm_type_final$TG))
        mm_type_final$TC <-
            as.numeric(as.character(mm_type_final$TC))
        mm_type_final$GG <-
            as.numeric(as.character(mm_type_final$GG))
        mm_type_final$CC <-
            as.numeric(as.character(mm_type_final$CC))
        
        ### Combine dataframes, append new variables, and reshape
        testdata <- cbind(mm_final, mm_type_final)
        testdata$Oligo <- as.character(substring(testdata$Oligo, 6))
        testdata$Center_mm <-
            as.integer(paste(
                testdata$Total_mm - testdata$End5p_mm - testdata$End3p_mm
            ))
        
        testdata <- testdata[,-c(7, 10:14)]
        
        ### Reshape dataframe to wide format
        testdata <- reshape(
            data = testdata,
            idvar = c("Assay", "Name"),
            timevar = "Oligo",
            v.names = c(
                "Total_mm",
                "End3p_mm",
                "Term_mm",
                "Center_mm",
                "AA",
                "AG",
                "AC",
                "TT",
                "TG",
                "TC",
                "GG",
                "CC"
            ),
            direction = "wide"
        )
        
        testdata$FRmm_total <-
            as.integer(paste(testdata$Total_mm.F + testdata$Total_mm.R))
        testdata$FRmm_diff <-
            as.numeric(paste(
                abs(testdata$Total_mm.F - testdata$Total_mm.R) / (testdata$Total_mm.F + testdata$Total_mm.R)
            ))
        testdata$FRmm_3p <-
            as.numeric(paste((testdata$End3p_mm.F + testdata$End3p_mm.R) / (testdata$Total_mm.F + testdata$Total_mm.R)
            )) # Proportion
        testdata$FRmm_term <-
            as.numeric(paste((testdata$Term_mm.F + testdata$Term_mm.R) / (testdata$Total_mm.F + testdata$Total_mm.R)
            )) # Proportion
        testdata$FRmm_AA <-
            as.numeric(paste((testdata$AA.F + testdata$AA.R) / (testdata$Total_mm.F + testdata$Total_mm.R)
            )) # Proportion
        testdata$FRmm_AG <-
            as.numeric(paste((testdata$AG.F + testdata$AG.R) / (testdata$Total_mm.F + testdata$Total_mm.R)
            )) # Proportion
        testdata$FRmm_AC <-
            as.numeric(paste((testdata$AC.F + testdata$AC.R) / (testdata$Total_mm.F + testdata$Total_mm.R)
            )) # Proportion
        testdata$FRmm_TT <-
            as.numeric(paste((testdata$TT.F + testdata$TT.R) / (testdata$Total_mm.F + testdata$Total_mm.R)
            )) # Proportion
        testdata$FRmm_TG <-
            as.numeric(paste((testdata$TG.F + testdata$TG.R) / (testdata$Total_mm.F + testdata$Total_mm.R)
            )) # Proportion
        testdata$FRmm_TC <-
            as.numeric(paste((testdata$TC.F + testdata$TC.R) / (testdata$Total_mm.F + testdata$Total_mm.R)
            )) # Proportion
        testdata$FRmm_GG <-
            as.numeric(paste((testdata$GG.F + testdata$GG.R) / (testdata$Total_mm.F + testdata$Total_mm.R)
            )) # Proportion
        testdata$FRmm_CC <-
            as.numeric(paste((testdata$CC.F + testdata$CC.R) / (testdata$Total_mm.F + testdata$Total_mm.R)
            )) # Proportion
        names(testdata)[names(testdata) == "Total_mm.P"] <-
            "Pmm_total"
        testdata$Pmm_center <-
            as.numeric(paste(testdata$Center_mm.P / testdata$Pmm_total)) # Proportion
        testdata$Pmm_AA <-
            as.numeric(paste(testdata$AA.P / testdata$Pmm_total)) # Proportion
        testdata$Pmm_AG <-
            as.numeric(paste(testdata$AG.P / testdata$Pmm_total)) # Proportion
        testdata$Pmm_AC <-
            as.numeric(paste(testdata$AC.P / testdata$Pmm_total)) # Proportion
        testdata$Pmm_TT <-
            as.numeric(paste(testdata$TT.P / testdata$Pmm_total)) # Proportion
        testdata$Pmm_TG <-
            as.numeric(paste(testdata$TG.P / testdata$Pmm_total)) # Proportion
        testdata$Pmm_TC <-
            as.numeric(paste(testdata$TC.P / testdata$Pmm_total)) # Proportion
        testdata$Pmm_GG <-
            as.numeric(paste(testdata$GG.P / testdata$Pmm_total)) # Proportion
        testdata$Pmm_CC <-
            as.numeric(paste(testdata$CC.P / testdata$Pmm_total)) # Proportion
        
        testdata$F_length <-
            length(input_matrix[1, 4:ncol(input_matrix)]) - sum(is.na(input_matrix[1, ]))
        testdata$R_length <-
            length(input_matrix[2, 4:ncol(input_matrix)]) - sum(is.na(input_matrix[2, ]))
        testdata$P_length <-
            length(input_matrix[3, 4:ncol(input_matrix)]) - sum(is.na(input_matrix[3, ]))
        testdata$FR_length <-
            as.numeric(paste((
                testdata$F_length + testdata$R_length
            ) / 2))
        testdata$FR_Tm <- (input$F_Tm + input$R_Tm) / 2
        testdata$FR_Tmdiff <- abs(input$F_Tm - input$R_Tm)
        testdata$P_Tm <- input$P_Tm
        
        testdata <- subset(
            testdata,
            select = c(
                "Assay",
                "Taxon",
                "Name",
                "FRmm_total",
                "FRmm_diff",
                "FRmm_3p",
                "FRmm_term",
                "FRmm_AA",
                "FRmm_AG",
                "FRmm_AC",
                "FRmm_TT",
                "FRmm_TG",
                "FRmm_TC",
                "FRmm_GG",
                "FRmm_CC",
                "FR_length",
                "FR_Tm",
                "FR_Tmdiff",
                "Pmm_total",
                "Pmm_center",
                "Pmm_AA",
                "Pmm_AG",
                "Pmm_AC",
                "Pmm_TT",
                "Pmm_TG",
                "Pmm_TC",
                "Pmm_GG",
                "Pmm_CC",
                "P_length",
                "P_Tm"
            )
        )
        
        testdata <-
            testdata[order(testdata$Assay, testdata$Taxon), ]
        
        is.nan.data.frame <- function(x)
            do.call(cbind, lapply(x, is.nan))
        testdata[is.nan(testdata)] <- 0
        
        ### Predict amplification
        prediction <-
            predict(trainmodel_taqman,
                    newdata = testdata,
                    type = "prob") # Predict results of test data
        prediction <- cbind(testdata[, 1:3], prediction[, 1])
        names(prediction)[4] <- "Amp"
        prediction <- prediction[order(-prediction$Amp), ]
        
        return(list(data = testdata, pred = prediction))
        
    })
    
    ### Table of amplification predictions
    output$table <- renderTable({
        req(input$metadata)
        req(input$alignment)
        out <- out_matrix()
        df <- out$pred
    }, digits = 3)
    
    ### Table of amplification predictions for download
    output$downloadData <- downloadHandler(
        filename = "eDNAssay_output.csv",
        content = function(file) {
            out <- out_matrix()
            write.csv(out$pred, file, row.names = FALSE)
        }
    )
    
    ### Table of model input data for download
    output$downloadData2 <- downloadHandler(
        filename = "eDNAssay_input.csv",
        content = function(file) {
            out <- out_matrix()
            write.csv(out$data, file, row.names = FALSE)
        }
    )
    
    ### Used below
    load("eDNAssay_optimal_thresholds.RData")
    
    ### Optimal threshold plot
    output$opt <- renderPlot({
        par(mar = c(6, 5, 0.5, 0.5))
        plot(
            opt_threshold,
            pch = "",
            xlab = "Cost ratio (FN:FP)",
            ylab = "Optimal threshold",
            cex.lab = 1.5,
            cex.axis = 1.25
        )
        
        lines(opt_threshold, lwd = 2)
        points(
            input$costratio,
            opt_threshold[input$costratio],
            pch = 19,
            col = "#024f94",
            cex = 2
        )
        box(col = "black", lwd = 2)
        text(50,
             0.45,
             paste("Optimal threshold =", opt_threshold[input$costratio]),
             cex = 1.75)
    })
    
    ### In vitro testing histogram
    output$hist <- renderPlot({
        req(input$metadata)
        req(input$alignment)
        
        out <- out_matrix()
        h <- hist(out$pred[, 4], breaks = 100, plot = F)
        clr <-
            ifelse(h$breaks >= opt_threshold[input$costratio], "red", "grey")
        
        par(mar = c(6, 5, 0.5, 0.5))
        hist(
            out$pred[, 4],
            breaks = 100,
            col = clr,
            main = NULL,
            xlab = "Assignment probability",
            cex.lab = 1.5,
            cex.axis = 1.25
        )
        
        lines(rep(opt_threshold[input$costratio], 2), c(0, 100))
        text(0.7, 3, paste(sum(out$pred[, 4] > opt_threshold[input$costratio]), "above\nthreshold"), cex = 1.75)
        box(col = "black", lwd = 2)
    })
    
}

##################################################################################################
### Run the application
shinyApp(ui = ui, server = server)
