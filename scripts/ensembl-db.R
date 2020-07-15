# Build sqllite db from ensembl gtf
in.gtf <- as.character(snakemake@input[["gtf_file"]])
out.db <- as.character(snakemake@output[["gene_db"]])

library("ensembldb")
ensembl.db <- ensDbFromGtf(gtf = in.gtf, outfile=out.db)
