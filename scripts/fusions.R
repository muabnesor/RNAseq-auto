#Get input and output from snakemake

fusion.files <- unlist(snakemake@input[["fusions"]])
sample.data <- as.character(snakemake@input[["sample_data"]])

output.file <- as.character(snakemake@output[["genes"]])



# Load the count data and the meta data for the samples
fusions_dir <- "/home/adaros/Documents/fusions//"
col_data_file <- "/home/adaros/Documents/samples.tsv"
colData <- read.table(sample.data, header = TRUE)
colData$group <- sub("-","_",colData$group)

sampleFiles <- fusion.files
sampleNames <- unlist(lapply(strsplit(sampleFiles, "/"), function(x){tail(x, n = 1)}))
sampleNames <- sub(".tsv", "", sampleNames)

sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = colData$condition,
                          group = colData$group)

#row.names(sampleTable) <- sampleTable$sampleName

sample_fusions <- list()
arriba.colnames = c("gene1", "gene2", "strand1(gene/fusion)", "strand2(gene/fusion)",
                    "breakpoint1", "breakpoint2", "site1", "site2", "type", "direction1", "direction2",
                    "split_reads1", "split_reads2", "discordant_mates", "coverage1", "coverage2", "confidence",
                    "closest_genomic_breakpoint1", "closest_genomic_breakpoint2", "filters", "fusion_transcript",
                    "reading_frame", "peptide_sequence", "read_identifiers")
for (row in 1:nrow(sampleTable)){
    #sample_fusions[[sample]] <- read.table(sampleTable[sample,"fileName"])
    file.name <- as.character(sampleTable[row,2])
    sample.name <- as.character(sampleTable[row, "sampleName"])
    sample_fusions[[sample.name]] <- read.table(file.name, col.names=arriba.colnames)
    }

# Gather fusions for all samples into one table
library("plyr")
fusions <- data.frame()

for (sample in names(sample_fusions)){
    fusionsPerSample <- data.frame(matrix(NA, nrow = 1, ncol = 0))

    for (row in sample_fusions[sample]){

        if ("high"%in%row$confidence){
            gene1 <- gsub("\\s*\\([^\\)]+\\)","",row$gene1)
            gene2 <- gsub("\\s*\\([^\\)]+\\)","",row$gene2)
            fusion.id <- paste0(gene1,"_", gene2)
            fusionsPerSample[fusion.id] <- TRUE
        }


    }

    fusions <- rbind.fill(fusions, fusionsPerSample)

}

row.names(fusions) <- names(sample_fusions)

fusions <- replace(fusions, is.na(fusions), FALSE)

#Count fusions for each group
perCondition <- list()
for (condition in unique(as.character(sampleTable$condition))) {
    samples <- which(sampleTable$condition == condition)
    fusionsPerCondition <- fusions[samples,]
    perCondition[[condition]] <- fusionsPerCondition
    }

fraction <- data.frame(sapply(perCondition, function(x){apply(x,2,sum)}))

# Order fusions by difference between groups
fraction$diff <- abs(fraction[,2]-fraction[,1])
fraction <- fraction[order(fraction$diff, decreasing=TRUE),]

#Get gene names of top 20 fusions
genes <- c()
for (names in row.names(fraction[1:20,])){
    gene.names <- unlist(strsplit(names,"_"))
    for (gene.name in gene.names){
        genes <- c(genes,unlist(strsplit(gene.name,",")))
    }
}

#Write gene names to file
write.table(as.data.frame(unique(genes)), file=output.file, quote=FALSE, sep=",", row.names=F)
