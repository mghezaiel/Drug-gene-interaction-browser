library(optparse)
library("biomaRt")

# ARGV
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input file path", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file path", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# INPUT GENE NAMES
input <- read.csv(opt$input,sep = ",", header = TRUE)

# BIOMART QUERY
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
attributes <- c('hgnc_symbol',"uniprotswissprot")

fetched_attributes <- getBM(attributes = attributes,
      filters = c('hgnc_symbol'),
      values = input$hgnc_symbol, 
      mart = ensembl)

# CLEAN DF
fetched_attributes$uniprotswissprot <- sapply(fetched_attributes$uniprotswissprot,
                                              function(x){ifelse(nchar(x)==0,NA,x)})
fetched_attributes <- na.omit(fetched_attributes)

# SAVE
write.csv(fetched_attributes,opt$output,row.names=FALSE)
