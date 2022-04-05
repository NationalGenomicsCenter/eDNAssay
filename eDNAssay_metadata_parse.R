##################################################################################################
### Script parses a metadata file for eDNAssay from an alignment that uses GenBank names
##################################################################################################
library(Biostrings)
library(stringr)

alignment <- readDNAStringSet(file.choose(), format = "fasta") # FASTA file
name <- alignment@ranges@NAMES
name <- gsub("NC ", "NC", name)

taxon <- word(name, start = 2, end=3, sep = fixed(" "))
taxon[1:3] <- rep("Target", 3)

type <- c(rep("Oligo", 3), rep("Template", length(name) - 3))

metadata <- data.frame(Taxon = taxon, Name = name, Type = type)

write.csv(metadata, "eDNAssay_metadata.csv", row.names = FALSE)
