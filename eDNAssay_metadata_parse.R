##################################################################################################
### Script parses a metadata file from an alignment, formatted for eDNAssay
##################################################################################################
library(Biostrings)

target_taxon <- "Target"

fasta <-
  readDNAStringSet(file.choose(), format = "fasta") # FASTA file

names <- fasta@ranges@NAMES

extract_name <- function(x) {
  print(gsub(
    "^[A-Z]+[^\\s]+\\s+([A-Z]+\\w+\\s+\\w+)\\s+.+",
    replacement = "\\1",
    x
  ))
}

Taxon <- extract_name(names)
Taxon[1:3] <- rep(target_taxon, 3)

Name <- names

Type <- c(rep("Oligo", 3), rep("Template", length(names) - 3))

metadata <- data.frame(Taxon = Taxon, Name = Name, Type = Type)
write.csv(metadata, "Metadata.csv", row.names = FALSE)
