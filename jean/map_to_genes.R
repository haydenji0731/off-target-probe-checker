run_blat <- function(dna_sequence) {
  library(httr)
  library(rvest)
  blat_url <- "https://genome.ucsc.edu/cgi-bin/hgBlat"
  response <- POST(blat_url, 
                   body = list(
                     userSeq = dna_sequence,
                     organism = "Human",
                     db = "hg38",
                     type = "BLAT's guess",
                     sort = "query, score",
                     output = "psl"
                   ),
                   encode = "form")
  blat_results <- content(response, as = "text")
  
  lines <- strsplit(blat_results, "\n")[[1]]
  psl_data <- lines[grep("^\\d+", lines)]  # Find result lines
  
  # Parse relevant fields (assuming PSL format)
  out <- do.call(rbind, lapply(1:length(psl_data), function(i) {
    fields <- strsplit(psl_data[i], "\\s+")[[1]]
    fields
  }))
  
  return(out)
}

dna_sequence <- 'CCAGTACTCCAATCATGATGCTGACAGTGGCTCTAGCTGA'
out <- run_blat(dna_sequence) ## consistent with online blat

out

####### map to gene 
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

i = 5
out[i,]
seqname = gsub('chr', '', out[i,14])
start = as.numeric(out[i,16])+1
end = as.numeric(out[i,17])
strand = out[i,9]
c(seqname, strand, start, end) ## to compare with online

result <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "ensembl_gene_id"),
  filters = c("chromosome_name", "start", "end"),
  values = list(seqname, start, end),
  mart = ensembl
)
print(result)


