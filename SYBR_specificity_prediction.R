##################################################################################################
### Script estimates the specificity of SYBR Green-based qPCR assays. Mismatch portion of script
### adapted from So et al. (2020; Pitfalls during in silico prediction of primer specificity for
### eDNA surveillance; Ecosphere, e03193); outputs the probability of being assigned to the
### "amplify" class for each template provided; accepts IUPAC ambiguities
##################################################################################################
library(Biostrings)
library(dplyr)
library(stringr)

### Input a FASTA file containing aligned sequences. Primer sequences must appear first, ordered
### as forward primer and then reverse primer; only IUPAC-approved characters are allowed 
### (A, C, G, T, M, R, W, S, Y, K, V, H, D, B, N, -, +, and .); dashes (from indels or sequences 
### not fully overlapping with the assay) are treated as Ns (any base) for a conservative estimate 
### of assay specificity
input_seqs <- readDNAStringSet(file.choose())

### Specify assay name
Assay <- "Assay_name"

### Specify melting temperatures most closely matching your reaction conditions
F_Tm <- 60
R_Tm <- 60

### Specify output CSV file
output_mismatches <- paste(Assay, "MMs.csv") 
output_probabilities <- paste(Assay, "APs.csv") 

##################################################################################################
### Combine matrices and save as a dataframe
### Compile metadata matrix
taxon <- input_seqs@ranges@NAMES
taxon <- gsub("\\.", "", taxon) 
taxon <- gsub("\\+", " ", taxon)
taxon <- gsub("NC ", "NC", taxon)
taxon <- gsub("UNVERIFIED: ", "", taxon)
taxon <- trimws(gsub("\\w*[0-9]+\\w*\\s*", "", taxon)) 
taxon <- word(taxon, start = 1, end=2, sep = fixed(" ")) 
taxon[1:2] <- rep("Target", 2)

name <- input_seqs@ranges@NAMES

type <- c(rep("Oligo", 2), rep("Template", length(taxon) - 2))

input_metadata <- data.frame(Taxon = taxon, Name = name, Type = type)

input_metadata[1, 2] = "F" 
input_metadata[2, 2] = "R" 
input_metadata <- as.matrix(input_metadata)

input_seqs <- as.matrix(input_seqs)

### Combine matrices and save as a dataframe
input_matrix <- cbind(input_metadata, input_seqs)
input_matrix <-
  as.data.frame(input_matrix, stringsAsFactors = FALSE)
input_matrix <- replace(input_matrix, input_matrix == "-", NA)

### Convert "NA" in template sequences to "N" to accommodate indels
input_matrix_names <- input_matrix[, 1:3] 
input_matrix_names[is.na(input_matrix_names)] <- "Unspecified" 

input_matrix_seqs <- input_matrix[, 4:ncol(input_matrix)] 

input_matrix <- cbind(input_matrix_names, input_matrix_seqs) 

### Define length variables and relabel rows
length_type <- length(input_matrix$Type)
length_oligo <- length(which(input_matrix$Type == "Oligo"))
length_template <- length(which(input_matrix$Type == "Template"))
row.names(input_matrix) <- c(1:length_type)

### Create separate dataframes for each oligo and remove nucleotides outside oligo binding sites
oligo_matrices <-
  lapply(1:length_oligo, function(x)
    input_matrix[input_matrix$Type == "Template", -which(is.na(input_matrix[x,]))])
names(oligo_matrices) = input_matrix$Name[1:length_oligo]

##################################################################################################
### Match oligo and template sequences (accounting for IUPAC ambiguity codes) and count the number
### of mismatches (=TRUE) in different positions; total, the last five basepairs on 3' end,
### and the last basepair on the 3' end (terminus)
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
          ot_match[[i]][1, -(1:3)] == "R" &
          ot_match[[i]][x, -(1:3)] == "A" |
          ot_match[[i]][1, -(1:3)] == "R" &
          ot_match[[i]][x, -(1:3)] == "G" |
          ot_match[[i]][1, -(1:3)] == "Y" &
          ot_match[[i]][x, -(1:3)] == "C" |
          ot_match[[i]][1, -(1:3)] == "Y" &
          ot_match[[i]][x, -(1:3)] == "T" |
          ot_match[[i]][1, -(1:3)] == "M" &
          ot_match[[i]][x, -(1:3)] == "A" |
          ot_match[[i]][1, -(1:3)] == "M" &
          ot_match[[i]][x, -(1:3)] == "C" |
          ot_match[[i]][1, -(1:3)] == "K" &
          ot_match[[i]][x, -(1:3)] == "G" |
          ot_match[[i]][1, -(1:3)] == "K" &
          ot_match[[i]][x, -(1:3)] == "T" |
          ot_match[[i]][1, -(1:3)] == "S" &
          ot_match[[i]][x, -(1:3)] == "G" |
          ot_match[[i]][1, -(1:3)] == "S" &
          ot_match[[i]][x, -(1:3)] == "C" |
          ot_match[[i]][1, -(1:3)] == "W" &
          ot_match[[i]][x, -(1:3)] == "A" |
          ot_match[[i]][1, -(1:3)] == "W" &
          ot_match[[i]][x, -(1:3)] == "T" |
          ot_match[[i]][1, -(1:3)] == "B" &
          ot_match[[i]][x, -(1:3)] == "C" |
          ot_match[[i]][1, -(1:3)] == "B" &
          ot_match[[i]][x, -(1:3)] == "G" |
          ot_match[[i]][1, -(1:3)] == "B" &
          ot_match[[i]][x, -(1:3)] == "T" |
          ot_match[[i]][1, -(1:3)] == "D" &
          ot_match[[i]][x, -(1:3)] == "A" |
          ot_match[[i]][1, -(1:3)] == "D" &
          ot_match[[i]][x, -(1:3)] == "G" |
          ot_match[[i]][1, -(1:3)] == "D" &
          ot_match[[i]][x, -(1:3)] == "T" |
          ot_match[[i]][1, -(1:3)] == "H" &
          ot_match[[i]][x, -(1:3)] == "A" |
          ot_match[[i]][1, -(1:3)] == "H" &
          ot_match[[i]][x, -(1:3)] == "C" |
          ot_match[[i]][1, -(1:3)] == "H" &
          ot_match[[i]][x, -(1:3)] == "T" |
          ot_match[[i]][1, -(1:3)] == "V" &
          ot_match[[i]][x, -(1:3)] == "A" |
          ot_match[[i]][1, -(1:3)] == "V" &
          ot_match[[i]][x, -(1:3)] == "C" |
          ot_match[[i]][1, -(1:3)] == "V" &
          ot_match[[i]][x, -(1:3)] == "G" |
          ot_match[[i]][1, -(1:3)] == "N" &
          ot_match[[i]][x, -(1:3)] == "A" |
          ot_match[[i]][1, -(1:3)] == "N" &
          ot_match[[i]][x, -(1:3)] == "C" |
          ot_match[[i]][1, -(1:3)] == "N" &
          ot_match[[i]][x, -(1:3)] == "G" |
          ot_match[[i]][1, -(1:3)] == "N" &
          ot_match[[i]][x, -(1:3)] == "T" |
          ot_match[[i]][1, -(1:3)] == "A" &
          ot_match[[i]][x, -(1:3)] == "R" |
          ot_match[[i]][1, -(1:3)] == "G" &
          ot_match[[i]][x, -(1:3)] == "R" |
          ot_match[[i]][1, -(1:3)] == "C" &
          ot_match[[i]][x, -(1:3)] == "Y" |
          ot_match[[i]][1, -(1:3)] == "T" &
          ot_match[[i]][x, -(1:3)] == "Y" |
          ot_match[[i]][1, -(1:3)] == "A" &
          ot_match[[i]][x, -(1:3)] == "M" |
          ot_match[[i]][1, -(1:3)] == "C" &
          ot_match[[i]][x, -(1:3)] == "M" |
          ot_match[[i]][1, -(1:3)] == "G" &
          ot_match[[i]][x, -(1:3)] == "K" |
          ot_match[[i]][1, -(1:3)] == "T" &
          ot_match[[i]][x, -(1:3)] == "K" |
          ot_match[[i]][1, -(1:3)] == "G" &
          ot_match[[i]][x, -(1:3)] == "S" |
          ot_match[[i]][1, -(1:3)] == "C" &
          ot_match[[i]][x, -(1:3)] == "S" |
          ot_match[[i]][1, -(1:3)] == "A" &
          ot_match[[i]][x, -(1:3)] == "W" |
          ot_match[[i]][1, -(1:3)] == "T" &
          ot_match[[i]][x, -(1:3)] == "W" |
          ot_match[[i]][1, -(1:3)] == "C" &
          ot_match[[i]][x, -(1:3)] == "B" |
          ot_match[[i]][1, -(1:3)] == "G" &
          ot_match[[i]][x, -(1:3)] == "B" |
          ot_match[[i]][1, -(1:3)] == "T" &
          ot_match[[i]][x, -(1:3)] == "B" |
          ot_match[[i]][1, -(1:3)] == "A" &
          ot_match[[i]][x, -(1:3)] == "D" |
          ot_match[[i]][1, -(1:3)] == "G" &
          ot_match[[i]][x, -(1:3)] == "D" |
          ot_match[[i]][1, -(1:3)] == "T" &
          ot_match[[i]][x, -(1:3)] == "D" |
          ot_match[[i]][1, -(1:3)] == "A" &
          ot_match[[i]][x, -(1:3)] == "H" |
          ot_match[[i]][1, -(1:3)] == "C" &
          ot_match[[i]][x, -(1:3)] == "H" |
          ot_match[[i]][1, -(1:3)] == "T" &
          ot_match[[i]][x, -(1:3)] == "H" |
          ot_match[[i]][1, -(1:3)] == "A" &
          ot_match[[i]][x, -(1:3)] == "V" |
          ot_match[[i]][1, -(1:3)] == "C" &
          ot_match[[i]][x, -(1:3)] == "V" |
          ot_match[[i]][1, -(1:3)] == "G" &
          ot_match[[i]][x, -(1:3)] == "V" |
          ot_match[[i]][1, -(1:3)] == "A" &
          ot_match[[i]][x, -(1:3)] == "N" |
          ot_match[[i]][1, -(1:3)] == "C" &
          ot_match[[i]][x, -(1:3)] == "N" |
          ot_match[[i]][1, -(1:3)] == "G" &
          ot_match[[i]][x, -(1:3)] == "N" |
          ot_match[[i]][1, -(1:3)] == "T" &
          ot_match[[i]][x, -(1:3)] == "N",
        TRUE,
        FALSE
      ))
  mm_unlist[[i]] <-
    data.frame(matrix(unlist(mm_list[[i]]), nrow = length_template, byrow =
                        TRUE))
  mm_unlist[c(FALSE, TRUE)] = lapply(mm_unlist[c(FALSE, TRUE)], function(x)
    rev(x))
  mm_unlist[[i]]$mm_term <-
    as.numeric(mm_unlist[[i]][, ncol(mm_unlist[[i]])] == "FALSE")
  mm_unlist[[i]]$mm_total <-
    rowSums(mm_unlist[[i]][, -ncol(mm_unlist[[i]])] == "FALSE", na.rm = TRUE)
  mm_unlist[[i]]$mm_end3p <-
    rowSums(mm_unlist[[i]][, ((ncol(mm_unlist[[i]]) - 2 - 4):(ncol(mm_unlist[[i]]) -
                                                                2))] == "FALSE")
}

names(mm_unlist) = input_matrix$Name[1:length_oligo]

##################################################################################################
### Create a dataframe for mismatch number and position data
mm_final <-
  as.data.frame(cbind(
    rep(sub(" .*", "", names(mm_unlist)), each = length_template),
    rep(input_matrix[1:length_oligo, 2], each =
          length_template),
    rep(input_matrix[(length_oligo + 1):length_type, 1], length_oligo),
    rep(input_matrix[(length_oligo + 1):length_type, 2], length_oligo),
    rep(input_matrix[(length_oligo + 1):length_type, 3], length_oligo),
    as.numeric(sapply(mm_unlist, `[[`, "mm_total")),
    as.numeric(sapply(mm_unlist, `[[`, "mm_end3p")),
    as.numeric(sapply(mm_unlist, `[[`, "mm_term"))
  ))

colnames(mm_final) <-
  c("Assay",
    "Oligo",
    "Taxon",
    "Name",
    "Type",
    "Total_mm",
    "End3p_mm",
    "Term_mm")
mm_final$Total_mm <- as.numeric(as.character(mm_final$Total_mm))
mm_final$End3p_mm <- as.numeric(as.character(mm_final$End3p_mm))
mm_final$Term_mm <- as.numeric(as.character(mm_final$Term_mm))

##################################################################################################
### Repeat the above loop for all mismatch types (AA, AG-GA, AC-CA, TT, TG-GT, TC-CT, CC, and GG);
### loop i used for forward primers (seq(1,length_oligo,2)) and j for reverse primers (i+1)
aa_list <- list() # Returns a list of lists with AA mismatches
aa_unlist <-
  list() # Converts aa_list to dataframes and counts AA mismatches in different positions

for (i in seq(1, length_oligo, 3)) {
  for (j in (i + 1)) {
    aa_list[[i]] <- lapply(2:(length_template + 1), function(x)
      ifelse(ot_match[[i]][1, -(1:3)] == "A" &
               ot_match[[i]][x, -(1:3)] == "T", TRUE, FALSE))
  }
  aa_list[[j]] <- lapply(2:(length_template + 1), function(x)
    ifelse(ot_match[[j]][1, -(1:3)] == "T" &
             ot_match[[j]][x, -(1:3)] == "A", TRUE, FALSE))
  aa_unlist[[i]] <-
    data.frame(matrix(unlist(aa_list[[i]]), nrow = length_template, byrow =
                        TRUE))
  aa_unlist[[j]] <-
    rev(data.frame(matrix(
      unlist(aa_list[[j]]), nrow = length_template, byrow = TRUE
    )))
  aa_unlist[[i]]$aa_term <-
    as.numeric(aa_unlist[[i]][, ncol(aa_unlist[[i]])] == "TRUE")
  aa_unlist[[j]]$aa_term <-
    as.numeric(aa_unlist[[j]][, ncol(aa_unlist[[j]])] == "TRUE")
  aa_unlist[[i]]$aa_total <-
    rowSums(aa_unlist[[i]][, -ncol(aa_unlist[[i]])] == "TRUE", na.rm = TRUE)
  aa_unlist[[j]]$aa_total <-
    rowSums(aa_unlist[[j]][, -ncol(aa_unlist[[j]])] == "TRUE", na.rm = TRUE)
  aa_unlist[[i]]$aa_end3p <-
    rowSums(aa_unlist[[i]][, ((ncol(aa_unlist[[i]]) - 2 - 4):(ncol(aa_unlist[[i]]) -
                                                                2))] == "TRUE")
  aa_unlist[[j]]$aa_end3p <-
    rowSums(aa_unlist[[j]][, ((ncol(aa_unlist[[j]]) - 2 - 4):(ncol(aa_unlist[[j]]) -
                                                                2))] == "TRUE")
}

names(aa_unlist) = input_matrix$Name[1:length_oligo]

#-------------------------------------
ag_list <-
  list() # Returns a list of lists with AG and GA mismatches
ag_unlist <-
  list() # Converts ag_list to dataframes and counts AG and GA mismatches in different positions

for (i in seq(1, length_oligo, 3)) {
  for (j in (i + 1)) {
    ag_list[[i]] <- lapply(2:(length_template + 1), function(x)
      ifelse(
        ot_match[[i]][1, -(1:3)] == "A" & ot_match[[i]][x, -(1:3)] == "C" |
          ot_match[[i]][1, -(1:3)] == "G" &
          ot_match[[i]][x, -(1:3)] == "T",
        TRUE,
        FALSE
      ))
  }
  ag_list[[j]] <- lapply(2:(length_template + 1), function(x)
    ifelse(
      ot_match[[j]][1, -(1:3)] == "T" & ot_match[[j]][x, -(1:3)] == "G" |
        ot_match[[j]][1, -(1:3)] == "C" &
        ot_match[[j]][x, -(1:3)] == "A",
      TRUE,
      FALSE
    ))
  ag_unlist[[i]] <-
    data.frame(matrix(unlist(ag_list[[i]]), nrow = length_template, byrow =
                        TRUE))
  ag_unlist[[j]] <-
    rev(data.frame(matrix(
      unlist(ag_list[[j]]), nrow = length_template, byrow = TRUE
    )))
  ag_unlist[[i]]$ag_term <-
    as.numeric(ag_unlist[[i]][, ncol(ag_unlist[[i]])] == "TRUE")
  ag_unlist[[j]]$ag_term <-
    as.numeric(ag_unlist[[j]][, ncol(ag_unlist[[j]])] == "TRUE")
  ag_unlist[[i]]$ag_total <-
    rowSums(ag_unlist[[i]][, -ncol(ag_unlist[[i]])] == "TRUE", na.rm = TRUE)
  ag_unlist[[j]]$ag_total <-
    rowSums(ag_unlist[[j]][, -ncol(ag_unlist[[j]])] == "TRUE", na.rm = TRUE)
  ag_unlist[[i]]$ag_end3p <-
    rowSums(ag_unlist[[i]][, ((ncol(ag_unlist[[i]]) - 2 - 4):(ncol(ag_unlist[[i]]) -
                                                                2))] == "TRUE")
  ag_unlist[[j]]$ag_end3p <-
    rowSums(ag_unlist[[j]][, ((ncol(ag_unlist[[j]]) - 2 - 4):(ncol(ag_unlist[[j]]) -
                                                                2))] == "TRUE")
}

names(ag_unlist) = input_matrix$Name[1:length_oligo]

#-------------------------------------
ac_list <-
  list() # Returns a list of lists with AC and CA mismatches
ac_unlist <-
  list() # Converts ac_list to dataframes and counts AC and CA mismatches in different positions

for (i in seq(1, length_oligo, 3)) {
  for (j in (i + 1)) {
    ac_list[[i]] <- lapply(2:(length_template + 1), function(x)
      ifelse(
        ot_match[[i]][1, -(1:3)] == "A" & ot_match[[i]][x, -(1:3)] == "G" |
          ot_match[[i]][1, -(1:3)] == "C" &
          ot_match[[i]][x, -(1:3)] == "T",
        TRUE,
        FALSE
      ))
  }
  ac_list[[j]] <- lapply(2:(length_template + 1), function(x)
    ifelse(
      ot_match[[j]][1, -(1:3)] == "T" & ot_match[[j]][x, -(1:3)] == "C" |
        ot_match[[j]][1, -(1:3)] == "G" &
        ot_match[[j]][x, -(1:3)] == "A",
      TRUE,
      FALSE
    ))
  ac_unlist[[i]] <-
    data.frame(matrix(unlist(ac_list[[i]]), nrow = length_template, byrow =
                        TRUE))
  ac_unlist[[j]] <-
    rev(data.frame(matrix(
      unlist(ac_list[[j]]), nrow = length_template, byrow = TRUE
    )))
  ac_unlist[[i]]$ac_term <-
    as.numeric(ac_unlist[[i]][, ncol(ac_unlist[[i]])] == "TRUE")
  ac_unlist[[j]]$ac_term <-
    as.numeric(ac_unlist[[j]][, ncol(ac_unlist[[j]])] == "TRUE")
  ac_unlist[[i]]$ac_total <-
    rowSums(ac_unlist[[i]][, -ncol(ac_unlist[[i]])] == "TRUE", na.rm = TRUE)
  ac_unlist[[j]]$ac_total <-
    rowSums(ac_unlist[[j]][, -ncol(ac_unlist[[j]])] == "TRUE", na.rm = TRUE)
  ac_unlist[[i]]$ac_end3p <-
    rowSums(ac_unlist[[i]][, ((ncol(ac_unlist[[i]]) - 2 - 4):(ncol(ac_unlist[[i]]) -
                                                                2))] == "TRUE")
  ac_unlist[[j]]$ac_end3p <-
    rowSums(ac_unlist[[j]][, ((ncol(ac_unlist[[j]]) - 2 - 4):(ncol(ac_unlist[[j]]) -
                                                                2))] == "TRUE")
}

names(ac_unlist) = input_matrix$Name[1:length_oligo]

#-------------------------------------
tt_list <- list() # Returns a list of lists with TT mismatches
tt_unlist <-
  list() # Converts tt_list to dataframes and counts TT mismatches in different positions

for (i in seq(1, length_oligo, 3)) {
  for (j in (i + 1)) {
    tt_list[[i]] <- lapply(2:(length_template + 1), function(x)
      ifelse(ot_match[[i]][1, -(1:3)] == "T" &
               ot_match[[i]][x, -(1:3)] == "A", TRUE, FALSE))
  }
  tt_list[[j]] <- lapply(2:(length_template + 1), function(x)
    ifelse(ot_match[[j]][1, -(1:3)] == "A" &
             ot_match[[j]][x, -(1:3)] == "T", TRUE, FALSE))
  tt_unlist[[i]] <-
    data.frame(matrix(unlist(tt_list[[i]]), nrow = length_template, byrow =
                        TRUE))
  tt_unlist[[j]] <-
    rev(data.frame(matrix(
      unlist(tt_list[[j]]), nrow = length_template, byrow = TRUE
    )))
  tt_unlist[[i]]$tt_term <-
    as.numeric(tt_unlist[[i]][, ncol(tt_unlist[[i]])] == "TRUE")
  tt_unlist[[j]]$tt_term <-
    as.numeric(tt_unlist[[j]][, ncol(tt_unlist[[j]])] == "TRUE")
  tt_unlist[[i]]$tt_total <-
    rowSums(tt_unlist[[i]][, -ncol(tt_unlist[[i]])] == "TRUE", na.rm = TRUE)
  tt_unlist[[j]]$tt_total <-
    rowSums(tt_unlist[[j]][, -ncol(tt_unlist[[j]])] == "TRUE", na.rm = TRUE)
  tt_unlist[[i]]$tt_end3p <-
    rowSums(tt_unlist[[i]][, ((ncol(tt_unlist[[i]]) - 2 - 4):(ncol(tt_unlist[[i]]) -
                                                                2))] == "TRUE")
  tt_unlist[[j]]$tt_end3p <-
    rowSums(tt_unlist[[j]][, ((ncol(tt_unlist[[j]]) - 2 - 4):(ncol(tt_unlist[[j]]) -
                                                                2))] == "TRUE")
}

names(tt_unlist) = input_matrix$Name[1:length_oligo]

#-------------------------------------
tg_list <-
  list() # Returns a list of lists with TG and GT mismatches
tg_unlist <-
  list() # Converts tg_list to dataframes and counts TG and GT mismatches in different positions

for (i in seq(1, length_oligo, 3)) {
  for (j in (i + 1)) {
    tg_list[[i]] <- lapply(2:(length_template + 1), function(x)
      ifelse(
        ot_match[[i]][1, -(1:3)] == "T" & ot_match[[i]][x, -(1:3)] == "C" |
          ot_match[[i]][1, -(1:3)] == "G" &
          ot_match[[i]][x, -(1:3)] == "A",
        TRUE,
        FALSE
      ))
  }
  tg_list[[j]] <- lapply(2:(length_template + 1), function(x)
    ifelse(
      ot_match[[j]][1, -(1:3)] == "A" & ot_match[[j]][x, -(1:3)] == "G" |
        ot_match[[j]][1, -(1:3)] == "C" &
        ot_match[[j]][x, -(1:3)] == "T",
      TRUE,
      FALSE
    ))
  tg_unlist[[i]] <-
    data.frame(matrix(unlist(tg_list[[i]]), nrow = length_template, byrow =
                        TRUE))
  tg_unlist[[j]] <-
    rev(data.frame(matrix(
      unlist(tg_list[[j]]), nrow = length_template, byrow = TRUE
    )))
  tg_unlist[[i]]$tg_term <-
    as.numeric(tg_unlist[[i]][, ncol(tg_unlist[[i]])] == "TRUE")
  tg_unlist[[j]]$tg_term <-
    as.numeric(tg_unlist[[j]][, ncol(tg_unlist[[j]])] == "TRUE")
  tg_unlist[[i]]$tg_total <-
    rowSums(tg_unlist[[i]][, -ncol(tg_unlist[[i]])] == "TRUE", na.rm = TRUE)
  tg_unlist[[j]]$tg_total <-
    rowSums(tg_unlist[[j]][, -ncol(tg_unlist[[j]])] == "TRUE", na.rm = TRUE)
  tg_unlist[[i]]$tg_end3p <-
    rowSums(tg_unlist[[i]][, ((ncol(tg_unlist[[i]]) - 2 - 4):(ncol(tg_unlist[[i]]) -
                                                                2))] == "TRUE")
  tg_unlist[[j]]$tg_end3p <-
    rowSums(tg_unlist[[j]][, ((ncol(tg_unlist[[j]]) - 2 - 4):(ncol(tg_unlist[[j]]) -
                                                                2))] == "TRUE")
}

names(tg_unlist) = input_matrix$Name[1:length_oligo]

#-------------------------------------
tc_list <-
  list() # Returns a list of lists with TC and CT mismatches
tc_unlist <-
  list() # Converts tc_list to dataframes and counts TC and CT mismatches in different positions

for (i in seq(1, length_oligo, 3)) {
  for (j in (i + 1)) {
    tc_list[[i]] <- lapply(2:(length_template + 1), function(x)
      ifelse(
        ot_match[[i]][1, -(1:3)] == "T" & ot_match[[i]][x, -(1:3)] == "G" |
          ot_match[[i]][1, -(1:3)] == "C" &
          ot_match[[i]][x, -(1:3)] == "A",
        TRUE,
        FALSE
      ))
  }
  tc_list[[j]] <- lapply(2:(length_template + 1), function(x)
    ifelse(
      ot_match[[j]][1, -(1:3)] == "A" & ot_match[[j]][x, -(1:3)] == "C" |
        ot_match[[j]][1, -(1:3)] == "G" &
        ot_match[[j]][x, -(1:3)] == "T",
      TRUE,
      FALSE
    ))
  tc_unlist[[i]] <-
    data.frame(matrix(unlist(tc_list[[i]]), nrow = length_template, byrow =
                        TRUE))
  tc_unlist[[j]] <-
    rev(data.frame(matrix(
      unlist(tc_list[[j]]), nrow = length_template, byrow = TRUE
    )))
  tc_unlist[[i]]$tc_term <-
    as.numeric(tc_unlist[[i]][, ncol(tc_unlist[[i]])] == "TRUE")
  tc_unlist[[j]]$tc_term <-
    as.numeric(tc_unlist[[j]][, ncol(tc_unlist[[j]])] == "TRUE")
  tc_unlist[[i]]$tc_total <-
    rowSums(tc_unlist[[i]][, -ncol(tc_unlist[[i]])] == "TRUE", na.rm = TRUE)
  tc_unlist[[j]]$tc_total <-
    rowSums(tc_unlist[[j]][, -ncol(tc_unlist[[j]])] == "TRUE", na.rm = TRUE)
  tc_unlist[[i]]$tc_end3p <-
    rowSums(tc_unlist[[i]][, ((ncol(tc_unlist[[i]]) - 2 - 4):(ncol(tc_unlist[[i]]) -
                                                                2))] == "TRUE")
  tc_unlist[[j]]$tc_end3p <-
    rowSums(tc_unlist[[j]][, ((ncol(tc_unlist[[j]]) - 2 - 4):(ncol(tc_unlist[[j]]) -
                                                                2))] == "TRUE")
}

names(tc_unlist) = input_matrix$Name[1:length_oligo]

#-------------------------------------
gg_list <- list() # Returns a list of lists with GG mismatches
gg_unlist <-
  list() # Converts gg_list to dataframes and counts GG mismatches in different positions

for (i in seq(1, length_oligo, 3)) {
  for (j in (i + 1)) {
    gg_list[[i]] <- lapply(2:(length_template + 1), function(x)
      ifelse(ot_match[[i]][1, -(1:3)] == "G" &
               ot_match[[i]][x, -(1:3)] == "C", TRUE, FALSE))
  }
  gg_list[[j]] <- lapply(2:(length_template + 1), function(x)
    ifelse(ot_match[[j]][1, -(1:3)] == "C" &
             ot_match[[j]][x, -(1:3)] == "G", TRUE, FALSE))
  gg_unlist[[i]] <-
    data.frame(matrix(unlist(gg_list[[i]]), nrow = length_template, byrow =
                        TRUE))
  gg_unlist[[j]] <-
    rev(data.frame(matrix(
      unlist(gg_list[[j]]), nrow = length_template, byrow = TRUE
    )))
  gg_unlist[[i]]$gg_term <-
    as.numeric(gg_unlist[[i]][, ncol(gg_unlist[[i]])] == "TRUE")
  gg_unlist[[j]]$gg_term <-
    as.numeric(gg_unlist[[j]][, ncol(gg_unlist[[j]])] == "TRUE")
  gg_unlist[[i]]$gg_total <-
    rowSums(gg_unlist[[i]][, -ncol(gg_unlist[[i]])] == "TRUE", na.rm = TRUE)
  gg_unlist[[j]]$gg_total <-
    rowSums(gg_unlist[[j]][, -ncol(gg_unlist[[j]])] == "TRUE", na.rm = TRUE)
  gg_unlist[[i]]$gg_end3p <-
    rowSums(gg_unlist[[i]][, ((ncol(gg_unlist[[i]]) - 2 - 4):(ncol(gg_unlist[[i]]) -
                                                                2))] == "TRUE")
  gg_unlist[[j]]$gg_end3p <-
    rowSums(gg_unlist[[j]][, ((ncol(gg_unlist[[j]]) - 2 - 4):(ncol(gg_unlist[[j]]) -
                                                                2))] == "TRUE")
}

names(gg_unlist) = input_matrix$Name[1:length_oligo]

#-------------------------------------
cc_list <- list() # Returns a list of lists with CC mismatches
cc_unlist <-
  list() # Converts cc_list to dataframes and counts CC mismatches in different positions

for (i in seq(1, length_oligo, 3)) {
  for (j in (i + 1)) {
    cc_list[[i]] <- lapply(2:(length_template + 1), function(x)
      ifelse(ot_match[[i]][1, -(1:3)] == "C" &
               ot_match[[i]][x, -(1:3)] == "G", TRUE, FALSE))
  }
  cc_list[[j]] <- lapply(2:(length_template + 1), function(x)
    ifelse(ot_match[[j]][1, -(1:3)] == "G" &
             ot_match[[j]][x, -(1:3)] == "C", TRUE, FALSE))
  cc_unlist[[i]] <-
    data.frame(matrix(unlist(cc_list[[i]]), nrow = length_template, byrow =
                        TRUE))
  cc_unlist[[j]] <-
    rev(data.frame(matrix(
      unlist(cc_list[[j]]), nrow = length_template, byrow = TRUE
    )))
  cc_unlist[[i]]$cc_term <-
    as.numeric(cc_unlist[[i]][, ncol(cc_unlist[[i]])] == "TRUE")
  cc_unlist[[j]]$cc_term <-
    as.numeric(cc_unlist[[j]][, ncol(cc_unlist[[j]])] == "TRUE")
  cc_unlist[[i]]$cc_total <-
    rowSums(cc_unlist[[i]][, -ncol(cc_unlist[[i]])] == "TRUE", na.rm = TRUE)
  cc_unlist[[j]]$cc_total <-
    rowSums(cc_unlist[[j]][, -ncol(cc_unlist[[j]])] == "TRUE", na.rm = TRUE)
  cc_unlist[[i]]$cc_end3p <-
    rowSums(cc_unlist[[i]][, ((ncol(cc_unlist[[i]]) - 2 - 4):(ncol(cc_unlist[[i]]) -
                                                                2))] == "TRUE")
  cc_unlist[[j]]$cc_end3p <-
    rowSums(cc_unlist[[j]][, ((ncol(cc_unlist[[j]]) - 2 - 4):(ncol(cc_unlist[[j]]) -
                                                                2))] == "TRUE")
}

names(cc_unlist) = input_matrix$Name[1:length_oligo]

##################################################################################################
### Create a dataframe for mismatch type data
mm_type_final <-
  as.data.frame(
    cbind(
      rep(sub(" .*", "", names(mm_unlist)), each = length_template),
      rep(input_matrix[1:length_oligo, 2], each =
            length_template),
      rep(input_matrix[(length_oligo + 1):length_type, 1], length_oligo),
      rep(input_matrix[(length_oligo + 1):length_type, 2], length_oligo),
      rep(input_matrix[(length_oligo + 1):length_type, 3], length_oligo),
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
  c("Assay",
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
    "CC")
mm_type_final$AA <- as.numeric(as.character(mm_type_final$AA))
mm_type_final$AG <- as.numeric(as.character(mm_type_final$AG))
mm_type_final$AC <- as.numeric(as.character(mm_type_final$AC))
mm_type_final$TT <- as.numeric(as.character(mm_type_final$TT))
mm_type_final$TG <- as.numeric(as.character(mm_type_final$TG))
mm_type_final$TC <- as.numeric(as.character(mm_type_final$TC))
mm_type_final$GG <- as.numeric(as.character(mm_type_final$GG))
mm_type_final$CC <- as.numeric(as.character(mm_type_final$CC))

##################################################################################################
### Combine dataframes, append new variables, and reshape
testdata <- cbind(mm_final, mm_type_final)
testdata <- testdata[,-c(9:13)]

testdata_names <- testdata[, 1:5] 
testdata_counts <- testdata[, 6:16] 
testdata_counts[is.na(testdata_counts)] <- 0 
testdata <- cbind(testdata_names, testdata_counts) 

### Reshape dataframe to wide format
testdata <- reshape(
  data = testdata,
  idvar = c("Name"),
  timevar = "Oligo",
  v.names = c(
    "Total_mm",
    "End3p_mm",
    "Term_mm",
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

testdata$Assay <- Assay 
testdata$FRmm_total <-
  as.integer(paste(testdata$Total_mm.F + testdata$Total_mm.R))
testdata$FRmm_diff <-
  as.numeric(paste(abs(
    testdata$Total_mm.F - testdata$Total_mm.R
  )/(testdata$Total_mm.F + testdata$Total_mm.R)))
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

F_length <-
  length(input_matrix[1, 5:ncol(input_matrix)]) - sum(is.na(input_matrix[1, ]))
R_length <-
  length(input_matrix[2, 5:ncol(input_matrix)]) - sum(is.na(input_matrix[2, ]))
FR_length <- as.numeric(paste((F_length + R_length) / 2))
FR_Tm <- as.numeric(paste((F_Tm + R_Tm) / 2))
FR_Tmdiff <- as.numeric(paste(abs(F_Tm - R_Tm)))
testdata <- cbind(testdata, FR_length, FR_Tm, FR_Tmdiff)

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
    "FR_Tmdiff"
  )
)

testdata <- testdata[order(testdata$Assay, testdata$Taxon), ]

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
testdata[is.nan(testdata)] <- 0

write.csv(testdata, output_mismatches, row.names = FALSE)

##################################################################################################
### Load training model and predict amplification
load("SYBR_trained_model.RData")

prediction <-
  predict(trainmodel_sybr, newdata = testdata, type = "prob") # Predict results of test data
prediction <- cbind(testdata[, 1:3], prediction[, 1])
names(prediction)[4] <- "Amp"
prediction <- prediction[order(-prediction$Amp), ]

write.csv(prediction, output_probabilities, row.names = FALSE)
print("Finished!")
