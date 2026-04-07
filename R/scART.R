#code for functions for scReadthrough R package
#File defines three functions: AddReadthroughAssay, GeneReadthrough, and ReadthroughGTF

#' Add readthrough assay
#'
#' @param seurat_object A seurat object
#' @param data.dir A directory
#' @param gene.column A number
#' @param cell.column A number
#' @param unique.features A boolean
#' @param strip.suffix A boolean
#' @return A seurat object with a readthrough assay and overall cell readthrough level
#' @export
AddReadthroughAssay <- function(seurat_object, data.dir, gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE) {
    seurat_object[["percent.Readthrough"]] <-
        PercentageFeatureSet(CreateSeuratObject(counts =
            Read10X(data.dir = data.dir, gene.column = gene.column, cell.column = cell.column,
            unique.features = unique.features, strip.suffix = strip.suffix)), pattern = "^Readthrough-")
    seurat_object[["percent.GeneCount"]] <-
        PercentageFeatureSet(CreateSeuratObject(counts =
            Read10X(data.dir = data.dir, gene.column = gene.column, cell.column = cell.column,
            unique.features = unique.features, strip.suffix = strip.suffix)), pattern = "^GeneCount")
    seurat_object$logRatio.Readthrough <-
        log((seurat_object$percent.Readthrough + 0.1) / (seurat_object$percent.GeneCount + 0.1))
    seurat_object <-
        seurat_object[ , colnames(Read10X(data.dir = data.dir, gene.column = gene.column,
        cell.column = cell.column, unique.features = unique.features, strip.suffix = strip.suffix))]
    seurat_object[['Readthrough']] <-
        CreateAssayObject(Read10X(data.dir = data.dir, gene.column = gene.column, cell.column = cell.column,
        unique.features = unique.features, strip.suffix = strip.suffix)[ , colnames(seurat_object@assays$RNA)])
    if (min(seurat_object[["percent.Readthrough"]]) == 0) {
        print('Some cells have no reads aligning to readthrough regions')
        print('Consider checking cells for low read counts')
    }
    return(seurat_object)
}

#' Calculate gene readthrough levels
#'
#' @param seurat_object A seurat object with a readthrough assay
#' @param group.by A grouping criterion
#' @param sort.by A sorting criterion
#' @return sorted data frame
#' @export
GeneReadthrough <- function(seurat_object, group.by, sort.by = mean) {
    print(paste('Number of cells in each group when grouping by', group.by))
    print(table(seurat_object[[group.by]]))
    if (min(table(seurat_object[[group.by]])) < 100) {
    print('It is not recommended to calculate gene readthrough levels for groups without many cells')
    }
    averages <- AggregateExpression(seurat_object, assays = 'Readthrough', features = NULL, return.seurat = FALSE,
        group.by = group.by, add.ident = NULL, normalization.method = "LogNormalize", scale.factor = 10000,
        margin = 1, verbose = TRUE)
    mean.data <- data.frame(averages)
    row.names(mean.data) <- gsub("-", "_", row.names(mean.data))
    mean.data$sort.by <- apply(mean.data, 1, sort.by, na.rm=TRUE)
    rc = length(grep("Readthrough_", row.names(mean.data)))
    cc = ncol(mean.data)
    m <- matrix(ncol = cc, nrow = rc)
    df <- data.frame(m, row.names = row.names(mean.data)[c(grep("Readthrough_", row.names(mean.data)))])
    colnames(df) <- colnames(mean.data)
    corr <- row.names(table(seurat_object[[group.by]]))
    corr2 <- corr
    for (i in 1:length(corr)) {
        corr[i] <- mean(seurat_object$nCount_RNA[which(seurat_object[[group.by]] == corr2[i])]) +
        mean(seurat_object$nCount_Readthrough[which(seurat_object[[group.by]] == corr2[i])])
    }
    corr <- as.numeric(corr)

    for (j in 1:length(colnames(mean.data))) {
      if (j == length(colnames(mean.data))) {
        for (i in 1:length(row.names(mean.data))) {
            r = length(grep("GeneCounts_", row.names(mean.data)[i]))
            if (r) {
                Mean = mean.data[i, j]
                rowname = gsub("GeneCounts", "Readthrough", row.names(mean.data)[i])
                df[rowname, j] = Mean
            }
        }
        } else {
        for (i in 1:length(row.names(mean.data))) {
          g = length(grep("Readthrough_", row.names(mean.data)[i]))
          if (g) {
            Readthrough = mean.data[i, j]
            rowname = gsub("Readthrough", "GeneCounts", row.names(mean.data)[i])
            GeneCount = mean.data[rowname, j]
            if (is.na(Readthrough)) {
                Readthrough = (table(seurat_object[[group.by]])[j]/100) * corr[j]/1000
            } else {
                Readthrough = Readthrough + (table(seurat_object[[group.by]])[j]/100) * corr[j]/1000
            }
            if (is.na(GeneCount)) {
                GeneCount = (table(seurat_object[[group.by]])[j]/100) * corr[j]/1000
            } else {
                GeneCount = GeneCount + (table(seurat_object[[group.by]])[j]/100) * corr[j]/1000
            }
            rowname = row.names(mean.data)[i]
            df[rowname, j] = log( Readthrough / GeneCount)
          }
        }
      }
    }
    commonrows = c()
     for (rowname in rownames(mean.data)) {
        if (!grepl("Readthrough", rowname)) next;
        otherrow = gsub("Readthrough", "GeneCounts", rowname)
        if (otherrow %in% rownames(mean.data)) {
           commonrows = append(commonrows, rowname)
        }
     }
    df_filtered = df[commonrows, ]
    df_sorted <- df_filtered[order(df_filtered$sort.by ,decreasing = TRUE) , ]
    return(df_sorted)
}

#' Create GTF file annotating readthrough regions
#'
#' @param input A path to a gtf file
#' @param output A path to an output gtf file name
#' @param chr.sizes A path to a chromosome sizes file
#' @param readthrough.gap A number
#' @param readthrough.window A number
#' @param min.window A number
#' @return void
#' @export
ReadthroughGTF <- function(input, output, chr.sizes = "", readthrough.gap = 2000, readthrough.window = 150000, min.window = 15000) {
    #data <- gtf_file
    data <- read.delim(input, stringsAsFactors = FALSE, header = FALSE, comment.char = "#")
    if (chr.sizes != "") {
        chr <- read.table(chr.sizes, stringsAsFactors = FALSE, header = FALSE)
    }
    data <- data[which(data$V3 == 'gene') ,]
    ngene <- nrow(data)
    print(paste('GTF file lists', ngene, 'genes'))
    data$V9 <- gsub("\"", "", data$V9)
    maxgap <- readthrough.gap + readthrough.window + 1000
    readdata <- data
    editdata <- data

    for (i in 1:nrow(data)) {
       if (data[i , 7] == "+") {
           suppressWarnings(after <- min(data[which(data$V4 > data[i, 5] & data$V1 == data[i , 1] & data$V7 == data[i , 7]) , 4]))
           suppressWarnings(afterbeg <- min(data[which(data$V5 > data[i, 5] & data$V1 == data[i , 1] & data$V7 == data[i , 7]) , 4]))
           if (afterbeg < data[i , 5]) {
               readdata[i , 4] <- 0
               readdata[i , 5] <- 0
           } else if (after - data[i, 5] > maxgap) {
               readdata[i , 4] <- data[i , 5] + readthrough.gap
               readdata[i , 5] <- data[i , 5] + readthrough.gap + readthrough.window
           } else {
               readdata[i , 4] <- data[i , 5] + readthrough.gap
               readdata[i , 5] <- data[i , 5] + (after - data[i , 5] - 1000)
           }
       } else {
           suppressWarnings(after <- max(data[which(data$V5 < data[i, 5] & data$V1 == data[i , 1] & data$V7 == data[i , 7]) , 5]))
           suppressWarnings(afterbeg <- max(data[which(data$V4 < data[i, 4] & data$V1 == data[i , 1] & data$V7 == data[i , 7]) , 5]))
           if (afterbeg > data[i , 4]) {
               readdata[i , 5] <- 0
               readdata[i , 4] <- 0
           } else if (data[i, 4] - after > maxgap) {
               readdata[i , 5] <- data[i , 4] - readthrough.gap
               if (readdata[i , 5] < 0) readdata[i , 5] = 0
               readdata[i , 4] <- data[i , 4] - readthrough.gap - readthrough.window
               if (readdata[i , 4] < 0) readdata[i , 4] = 0
           } else {
               readdata[i , 5] <- data[i , 4] - readthrough.gap
               if (readdata[i , 5] < 0) readdata[i , 5] = 0
               readdata[i , 4] <- data[i , 4] - (data[i, 4] - after - 1000)
               if (readdata[i , 4] < 0) readdata[i , 4] = 0
           }
       }
    }
    if (chr.sizes != "") {
        for (i in 1:nrow(readdata)) {
            if (readdata[i , 5] > chr[which(chr$V1 == readdata[i , 1]) , 2]) {
                readdata[i , 5] <- chr[which(chr$V1 == readdata[i , 1]) , 2]
            }
        }
    }
    editdata <- editdata[which(readdata$V5 - readdata$V4 > min.window), ]
    readdata <- readdata[which(readdata$V5 - readdata$V4 > min.window), ]
    nreadthrough <- nrow(readdata)
    print(paste(nreadthrough, 'genes will be used for readthrough calculations'))
    editdata$V9 <- gsub("gene_id ", "exon_number 1; gene_id GeneCounts-", editdata$V9)
    editdata$V9 <- gsub("gene_name ", "gene_name GeneCounts-", editdata$V9)
    readdata$V9 <- gsub("gene_id ", "exon_number 1; gene_id Readthrough-", readdata$V9)
    readdata$V9 <- gsub("gene_name ", "gene_name Readthrough-", readdata$V9)
    editdataexon <- editdata
    editdataexon$V3 <- gsub("gene", "exon", editdataexon$V3)
    editdatatranscript <- editdata
    editdatatranscript$V3 <- gsub("gene", "transcript", editdatatranscript$V3)
    readdataexon <- readdata
    readdataexon$V3 <- gsub("gene", "exon", readdataexon$V3)
    readdatatranscript <- readdata
    readdatatranscript$V3 <- gsub("gene", "transcript", readdatatranscript$V3)
    final <- rbind(editdata, editdataexon, editdatatranscript, readdata, readdataexon, readdatatranscript)
    #return(final)
    write.table(final, file = output, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}
