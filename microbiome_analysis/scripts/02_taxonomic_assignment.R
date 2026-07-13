###SECONDARY SEQUENCE ANALYSIS###

# Taxonomic assignment
# Make sure the SILVA reference files are downloaded
# and located in the working directory
tax_bact <- assignTaxonomy(seqtab.nochim_bact, 
                           "silva_nr99_v138.1_train_set.fa.gz", 
                           outputBootstraps = TRUE, minBoot = 80)

tax_bact[["tax"]] <- addSpecies(tax_bact[["tax"]], 
                                "silva_species_assignment_v138.1.fa.gz", n = 100)

# Save output files to the working directory
write.table(tax_bact[["tax"]], 
            "taxa_bact.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

write.table(seqtab.nochim_bact, 
            "ASV_bact.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

# Export representative ASV sequences in FASTA format
uniquesToFasta(seqtab.nochim_bact, 
               fout = "rep.seqs_bact.fna", ids = colnames(seqtab.nochim_bact))

# Save the taxonomic assignment object as an RDS file
saveRDS(tax_bact, "tax_bact.rds")