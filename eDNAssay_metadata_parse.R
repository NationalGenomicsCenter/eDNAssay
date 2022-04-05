##################################################################################################
### Script parses a metadata file for eDNAssay from an alignment that uses GenBank names
##################################################################################################
library(Biostrings)
library(stringr)

alignment <- readDNAStringSet(file.choose(), format = "fasta") # FASTA file

taxon <- alignment@ranges@NAMES
taxon <- gsub("NC ", "NC", taxon)
taxon <- gsub("UNVERIFIED: ", "", taxon)
taxon <- word(taxon, start = 2, end=3, sep = fixed(" "))
taxon[1:3] <- rep("Target", 3)

name <- alignment@ranges@NAMES

type <- c(rep("Oligo", 3), rep("Template", length(name) - 3))

metadata <- data.frame(Taxon = taxon, Name = name, Type = type)

write.csv(metadata, "eDNAssay_metadata.csv", row.names = FALSE)
