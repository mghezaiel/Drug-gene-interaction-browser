library(optparse)
library("biomaRt")
library("UniProt.ws")
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input file path", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file path", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input <- read.csv(opt$input,sep = ",", header = TRUE)
results <- data.frame()

columns <- c("gene_names","annotation_score","organism_id","organism_name","cc_function","cc_cofactor","cc_tissue_specificity","cc_interaction","protein_existence","cc_allergen","protein_families","cc_biotechnology","redox_potential")
kt <- c("UniProtKB")

for (protein_id in 1:length(input$uniprot_accession))
    {
    # Protein name
    uniprot_accession <- input$uniprot_accession[protein_id]
    taxid <- UniProt.ws(as.integer(input$taxid[protein_id]))
    
    # Try to fetch uniprot infos
    res <- select(taxid, uniprot_accession, columns, kt)
    # Add the uniprotkb identifier
    res$uniprotkb <- uniprot_accession
    
    #print(dim(results))
    # Concatenate the results to the output dataframe
    results <- rbind(results,res)
}

print(results)
write.csv(results,opt$output)
