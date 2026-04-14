# scART
scART is a tool for analyzing readthrough transcription levels in single cell RNA sequencing (scRNA-seq) data.
Included are instructions for System Setup/Installation and an example Workflow.

## Step 0: System Setup and Installation
Make sure to have star, samtools, sra-tools, r-base, and r-seurat installed. You can do this using conda with:

```
conda install star samtools sra-tools
```
```
conda install r-base r-patchwork r-tidyverse r-seurat
```

You can also install using mamba with:

```
conda install mamba
```
```
mamba install star samtools sra-tools
```
```
mamba install r-base r-patchwork r-tidyverse r-seurat
```

Next, start R and install the 'remotes' package.

```
install.packages("remotes")
```

Finally, download the scART package:

```
remotes::install_github("nvshende/scART")
```

## Step 1: Generating a STAR index for your organism of interest
In order to generate a STAR index, you will need a genome file and gene anotation (gtf) file for the organism of interest.
You can download a human genome and gtf file using:

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
gunzip hg38.fa.gz
gunzip gencode.v39.annotation.gtf.gz
```

Next, create a directory for your STAR index:

```
mkdir Human-STAR-Index
```

Create a STAR index with:

```
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir Human-STAR-index/ --genomeFastaFiles hg38.fa --sjdbGTFfile gencode.v39.annotation.gtf --genomeSAindexNbases 10
```

You may increase or decrease --runThreadN depending on availability of computational resources.

## Step 2: Generating a readthrough gtf file
In order to generate readthrough gtf file, you will need a gtf file for your organism and genome.
This gtf file must have the chromosome name in the first column, the feature type in the third column,
the start coordinate in the fourth column, the end coordinate in the fifth column,
the strand (+ or -) in the seventh column, and additional information, such as gene_name and gene_id in the ninth column.
Readthrough regions will be defined based on features labelled genes in the third column.
Optionally, you can also have a chromosome sizes file in order to ensure that readthrough detection regions
do not extend beyond chromosome boundaries. If you do not have a chromosome sizes file, you can generate one using
samtools and the following code:

```
samtools faidx hg38.fa
cut -f1,2 hg38.fa.fai > chr.sizes
```

You can create a readthrough gtf file using the ReadthroughGTF function.
For example, you can use the following code in R:

```
library(scART)
ReadthroughGTF(input = "gencode.v39.annotation.gtf", output = "readthrough.gtf", chr.sizes = "chr.sizes")
```

The input and output parameters are required. The chr.sizes parameter is optional.
Other optional parameters include:
readthrough.gap the gap between the end of a gene and the start of the readthrough region (default = 2000)
readthrough.window the maximum size of the readthrough region (default = 150000)
min.window the minimum size of the readthrough region (default = 15000)
Exit R and use the new gtf file to generate a new alignment index:

```
mkdir Human-Readthrough-STAR-Index
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir Human-Readthrough-STAR-index/ --genomeFastaFiles hg38.fa --sjdbGTFfile readthrough.gtf --genomeSAindexNbases 10
```


## Step 3: Performing an Alignment
Download scRNA-seq data. Influenza infection scRNA-seq data will be used for this example:

```
fastq-dump --gzip --split-files SRR10832416
```

Downloading data with fastq-dump may take a long time.
Optionally, you can install parallel-fatq-dump to do this faster:

```
conda install parallel-fastq-dump
parallel-fastq-dump --sra-id SRR10832416 --threads 4 --outdir . --split-files --gzip
```

You may increase or decrease --threads depending on availability of computational resources.
Next, align your scRNA seq data to your genome using CellRanger or STARsolo.
In order to perform an scRNA-seq alignment, you will need to download a whitelist of cell barcodes.
This whitelist will be different for different scRNA-seq protocols.
The data in this example was generated using 10X Chromium Single Cell 3' V3 chemistry.
Download the corresponding whitelist using:

```
wget https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz
gunzip 3M-february-2018.txt.gz
```

Do the alignmnet using STAR.
Perform one alignment using the regular and one using the readthrough alignment index:

```
STAR --runThreadN 4 --genomeDir Human-STAR-index/ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --soloType Droplet --soloCBwhitelist 3M-february-2018.txt --soloBarcodeReadLength 28 --readFilesIn SRR10832416_2.fastq.gz SRR10832416_1.fastq.gz --outFileNamePrefix human_alignment.
STAR --runThreadN 4 --genomeDir Human-Readthrough-STAR-index/ --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --soloType Droplet --soloCBwhitelist 3M-february-2018.txt --soloBarcodeReadLength 28 --readFilesIn SRR10832416_2.fastq.gz SRR10832416_1.fastq.gz --outFileNamePrefix readthrough_alignment.
```
You may increase or decrease --runThreadN depending on availability of computational resources.


## Step 4: Loading your data into a Seurat Object
This requires the dplyr, Seurat, and patchwork R packages to be installed.
Each alignment should give a matrix.mtx file, a features.tsv file, and a barcodes.tsv file in a separate directory.
Ensure that each of these files is gzipped:

```
gzip human_alignment.Solo.out/Gene/filtered/*
gzip readthrough_alignment.Solo.out/Gene/filtered/*
```

The final file names should be matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz
You should have one directory with all three files from your regular alignment and
another directory with all three files from your reathrough alignment.
Start R, and load the dplyr, Seurat, and patchwork R packages
Load the data from your regular alignment into Seurat using the Read10X or ReadMtx function.
Create a Seurat object using the CreateSeuratObject command.
Add an assay with readthrough data to your Seurat object using the AddReadthroughAssay function.
For example, you can use the following code in R:

```
library(dplyr)
library(Seurat)
library(patchwork)
library(scART)
data <- Read10X(data.dir = "human_alignment.Solo.out/Gene/filtered/")
old_seurat_object <- CreateSeuratObject(counts = data, project = "your_project")
new_seurat_object <- AddReadthroughAssay(old_seurat_object, data.dir = "readthrough_alignment.Solo.out/Gene/filtered/")
```

This should create a Seurat object with regular scRNA seq data and readthrough data in a separate assay.
There will also be a metadata column called logRatio.Readthrough which has the overall readthrough levels for each cell.


## Step 5: Visualizing Readthrough Data
You can see a violin plot of the overall readthrough levels using the following code:

```
VlnPlot(new_seurat_object, features = "logRatio.Readthrough")
```

You can run UMAP on your data using:

```
new_seurat_object <- NormalizeData(new_seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
new_seurat_object <- FindVariableFeatures(new_seurat_object, selection.method = "vst", nfeatures = 2000)
new_seurat_object <- ScaleData(new_seurat_object)
new_seurat_object <- RunPCA(new_seurat_object, features = VariableFeatures(object = new_seurat_object))
new_seurat_object <- FindNeighbors(new_seurat_object, dims = 1:30)
new_seurat_object <- FindClusters(new_seurat_object, resolution = 0.3)
new_seurat_object <- RunUMAP(new_seurat_object, dims = 1:30, n.neighbors=7, n.components=2, min.dist=0.5)
```

After you run UMAP on your data, you can visualize the overall readthrough levels using:

```
FeaturePlot(new_seurat_object, features = "logRatio.Readthrough")
```

You can also get a data frame containing individual readthrough levels for each gene using the GeneReadthrough function.
This function takes a seurat object with a 'Readthrough' assay and a category by which to group the cells.
In this example, the cells will be grouped by seurat cluster.
The optional sort.by parameter lets you sort the genes by a parameter. The default value is mean.
Alternatives include median, min, and max.

```
df_sorted <- GeneReadthrough(new_seurat_object, group.by = 'seurat_clusters', sort.by = mean)
```

You can plot curves showing the distribution of readthrough levels for the top 1000 most highly expressed genes using the following code:

```
df_final <- df_sorted[1:1000, -ncol(df_sorted)]
plot(density(df_final$Readthrough.0), xlab = "Log Ratio Downstream Region vs. Gene", main = "Distribution of Readthrough Levels for Top 1000 Most Highly Expressed Genes", col = "#FF6666", ylim = c(0,0.8))
lines(density(df_final$Readthrough.1), col = "#CC9933")
lines(density(df_final$Readthrough.2), col = "#00CC33")
lines(density(df_final$Readthrough.3), col = "#0099FF")
lines(density(df_final$Readthrough.4), col = "#FF66FF")
legend("topleft", c("Seurat Cluster 0","Seurat Cluster 1", "Seurat Cluster 2", "Seurat Cluster 3", "Seurat Cluster 4"), col =c("#FF6666","#CC9933", "#00CC33", "#0099FF", "#FF66FF"), lty=1)
```
