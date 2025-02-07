library(Biostrings)
fasta_file <- "~/Downloads/xenium_human_breast_gene_expression_panel_probe_sequences.fasta"
dna_sequences_all <- readDNAStringSet(fasta_file)

## subset to a gene
g <- 'CD8A'
dna_sequences <- dna_sequences_all[grepl(g, names(dna_sequences_all)),]
dna_sequences

results_all <- lapply(1:length(dna_sequences), function(k) {
  print(k)
  target <- names(dna_sequences[k,])
  print(target)
  
  seq <- as.character(dna_sequences[k,])
  
  tryCatch(
  {
    out <- run_blat(seq) 
  },
  error = function(cond) {
    message(conditionMessage(cond))
    NA
  })
  
  ## by probe
  results <- do.call(rbind, lapply(1:nrow(out), function(i) {
    #print(i)
    seqname = gsub('chr', '', out[i,14])
    start = as.numeric(out[i,16])+1
    end = as.numeric(out[i,17])
    strand = out[i,9]
    result <- getBM(
      attributes = c("chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "ensembl_gene_id"),
      filters = c("chromosome_name", "start", "end"),
      values = list(seqname, start, end),
      mart = ensembl
    )
    
    return(result)
  }))
  tryCatch(
    {
    results$target = target
    },
    error = function(cond) {
      message(conditionMessage(cond))
  })
  
  return(results)
  
})
#save(results_all, file="~/Desktop/xenium_offtarget/results_all.RData")

## summarize
final <- do.call(rbind, results_all)
final$target_gene <- do.call(rbind, strsplit(final$target, '[|]'))[,1]
final

## matches
sum(final$ensembl_gene_id == unique(final$target_gene))

## not matches
final[final$ensembl_gene_id != unique(final$target_gene),]

## how many not matches are expressed
badhits <- na.omit(final[final$ensembl_gene_id != unique(final$target_gene),]$external_gene_name)
badhits.have <- unique(badhits[badhits %in% rownames(visium$counts)])
badhits.have

## visualize composite
visium_df <- data.frame(visium$pos, gexp=colSums(visium$counts[c(g, badhits.have),]))
ggplot(visium_df, aes(x=x, y=y, col=gexp)) + geom_point(size=1, alpha=1) + coord_fixed() + 
  scale_color_gradient(low='lightgrey', high='red') +
  theme_void() + ggtitle('composite')

