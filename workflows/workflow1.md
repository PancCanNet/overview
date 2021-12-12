# Workflow 1: Exploration PancCanNet pathway collection

* Download .Rmd file for the workflow [here](PancCanNet-workflow1.Rmd). 
* Open the .Rmd file in RStudio (see installation guide here) and walk through one section after the other. 

## Introduction

In this workflow, you will explore the PancCanNet pathway resource in
detail and study the pathways related to pancreatic cancer. You will
integrate knowledge from the Human Protein Atlas and TCGA about
pancreatic cancer prognostic markers and pancreatic cancer specific
genes. Lastly, you can visualize your gene or protein expression data in
the context of pancreatic cancer related processes.

The PancCanNet project is a collaboration between the Erasmus MC,
Maastricht University and Omnigen. The pathway portal is integrated in
the collaborative pathway database WikiPathways
(<http://panccannet.wikipathways.org>) and all analysis workflows are
available on Github (<https://panccannet.github.io>).

Pathway curation is an ongoing effort and we are continuously improving
current pathway models and adding new pathway models. The workflow will
always use the most-up-to-date information and results can slightly
differ to the information on the project website.

## R environment setup

First, we need to make sure all required R-packages are installed. Run
the following code snippet to install and load the libraries.

``` r
# check if libraries are installed > otherwise install
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")
if(!"rWikiPathways" %in% installed.packages()) BiocManager::install("rWikiPathways")
if(!"RCy3" %in% installed.packages()) BiocManager::install("RCy3")
if(!"data.table" %in% installed.packages()) BiocManager::install("data.table")
if(!"heatmaply" %in% installed.packages()) BiocManager::install("heatmaply")
if(!"igraph" %in% installed.packages()) BiocManager::install("igraph")
if(!"RColorBrewer" %in% installed.packages()) BiocManager::install("RColorBrewer")

# load required libraries
library(rstudioapi)
library(rWikiPathways)
library(RCy3)
library(data.table)
library(heatmaply)
library(igraph)
library(RColorBrewer)

# set working environment to source file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
```

## Explore PancCanNet pathway resources

The following section will download the current PancCanNet pathway
collection and provide basic statistics and information about the
pathways.

``` r
# retrieve pathway list
pwys <- rWikiPathways::getPathwayIdsByCurationTag("Curation:PancCanNet")
print(paste0("The PancCanNet pathway portal currently contains ",length(pwys)," pathways."))

# retrieve gene list per pathway
mylist <- list()
x <- 1
for(p in pwys) {
  genes <- rWikiPathways::getXrefList(p,"En")
  for(g in genes) {
    mylist[[x]] <- c(p,g)
    x <- x+1
  }
}
pwy.genes <- as.data.frame(do.call("rbind",mylist))
rm(mylist,x,genes,p,g)

# report pathway size distribution
hist(table(pwy.genes$V1), breaks=15, xlab="Pathway size (num of genes)", ylab = "Num of pathways", main="Pathway size distribution", col="#FFAAAA",border="#550000")

# report unique number of genes in PancCanNet pathways
print(paste0("The pathways in the PancCanNet collection currently contain ",length(unique(pwy.genes$V2))," unique genes."))
```

---

Let’s explore the overlap and crosstalk between the pathways in the
collection. We will first calculate the Jaccard similarity score (based
on gene content) and additionally, we will visualize the pathway-gene
network in Cytoscape which allows us to see which genes are present in
many pathways.

``` r
# calculate Jaccard index
jaccard <- matrix(ncol=length(pwys), nrow=length(pwys))
for(i in 1:length(pwys)) {
  for(j in i:length(pwys)) {
    genes1 <- pwy.genes[pwy.genes$V1 == pwys[i],2]
    genes2 <- pwy.genes[pwy.genes$V1 == pwys[j],2]
    intersection = length(intersect(genes1, genes2))
    union = length(genes1) + length(genes2) - intersection
    score <- (intersection/union)
    jaccard[i,j] <- score
    jaccard[j,i] <- score
  }
}

colnames(jaccard) <- pwys
rownames(jaccard) <- pwys

heatmaply(jaccard,main="Similarity plot PancCanNet pathways",plot_method = "plotly", colors = viridis(n = 256,  option = "magma"))

rm(genes1,genes2,i,j,intersection,union, score)
```

**Jaccard similarity heatmap:** 

Lighter areas in the heatmap show pathways with similar gene content (high similarity score). Overall, the
PancCanNet collections contains pathways with very different gene
content (low similarity score).

For the next step, please make sure that you have Cytoscape v3.8.0+
running. We want to explore the cross-talk between the pathway and see
which central genes are present in many different pathways.

``` r
# create pathway-gene network
graph <- igraph::graph_from_edgelist(as.matrix(pwy.genes), directed = TRUE)
suid <- RCy3::createNetworkFromIgraph(graph, title = "PancCanNet pathway gene network", collection = "PancCanNet")

table <- RCy3::getTableColumns(columns = c("name"))
table$type <- ifelse(startsWith(table$name, "WP"), 'Pathway', 'Gene')
RCy3::loadTableData(table, data.key.column = "name", table.key.column = "name")

RCy3::mapTableColumn("name",species="Human",map.from="Ensembl", map.to = "HGNC")

# analyse network properties
RCy3::analyzeNetwork(directed=TRUE)

# adapt visual style
toggleGraphicsDetails()
RCy3::copyVisualStyle("default","crosstalk")
RCy3::setVisualStyle("crosstalk")

RCy3::deleteStyleMapping("crosstalk","NODE_LABEL")
RCy3::lockNodeDimensions(TRUE, style.name="crosstalk")
RCy3::setNodeShapeMapping('type', c('Gene','Pathway'), c("ellipse","hexagon"), style.name="crosstalk")

RCy3::setNodeSizeMapping('type', c('Gene','Pathway'), c(25,40), mapping.type = "d", style.name = "crosstalk")

setNodeColorMapping("type", c('Gene','Pathway'), c("#EBCD9E","#70A689"), mapping.type = 'd', default.color = "#99FF99", style.name = "crosstalk")

# visualize hub nodes in pathway gene network
suid.hubs <- RCy3::cloneNetwork(suid)
RCy3::setCurrentNetwork(suid.hubs)
RCy3::renameNetwork("PancCanNet pathway-gene network: Hubs")

# adapt visual style
toggleGraphicsDetails()
RCy3::copyVisualStyle("crosstalk","hubs")
RCy3::setVisualStyle("hubs")

# identify most connected genes
degree <- RCy3::getTableColumns(columns = c("name","HGNC","Indegree"))
degree <- degree[degree$Indegree > 0,]
colnames(degree) <- c("Gene","GeneName","NumPathways")
degree <- degree[order(-degree$NumPathways),]

hist(degree$NumPathways, breaks=14, xlab="Node degree (genes)", ylab = "Num of genes", main="Distribution number of pathways", col="#FFAAAA",border="#550000",xlim=c(1,15))

hub.node.cutoff <- 10

hubs <- degree[degree$NumPathways > hub.node.cutoff,]
hubs

RCy3::setNodeBorderColorBypass(hubs$Gene, new.colors = "#052549")
RCy3::setNodeBorderWidthBypass(hubs$Gene, new.sizes = 10)

rm(runningRemote,graph,table,notebookIsRunning)
```

Explore the network in Cytoscape. You can see a lot of crosstalk and
overlap between the pathways in the portal. The hub nodes are
highlighted with a thick border. You can change the hub degree cutoff if
you want to be more strict or less strict. Even though the Jaccard index
indicated that there are very few similar pathways, there is still a lot
of crosstalk between the pathways in the collection.

## Integrate basic pancreatic cancer knowledge

In this section, you will integrate current knowledge about pancreatic
cancer from the Human Protein Atlas and TCGA. Let’s first investigate
known prognostic markers for pancreatic cancer and assess their presence
in the PancCanNet pathways.

``` r
# favorable prognostic markers
link <- "https://www.proteinatlas.org/api/search_download.php?search=prognostic%3Apancreatic%20cancer%3BFavorable%20AND%20sort_by%3Aprognostic%20pancreatic%20cancer&columns=g,eg,up,rnacas,rnacad,rnacass,rnacasm,prognostic_pancreatic_cancer,t_RNA_pancreas,sc_RNA_Pancreatic_endocrine_cells&compress=no&format=tsv"
fav.markers <- read.table(file = link, sep = '\t', header = TRUE)

cat(paste0("There are currently ", length(fav.markers$Ensembl), " known favorable prognostic markers for pancreatic cancer."))

# unfavorable prognostic markers
link <- "https://www.proteinatlas.org/api/search_download.php?search=prognostic%3Apancreatic%20cancer%3BUnfavorable%20AND%20sort_by%3Aprognostic%20pancreatic%20cancer&columns=g,eg,up,rnacas,rnacad,rnacass,rnacasm,prognostic_pancreatic_cancer,t_RNA_pancreas,sc_RNA_Pancreatic_endocrine_cells&compress=no&format=tsv"
unfav.markers <- read.table(file = link, sep = '\t', header = TRUE)

cat(paste0("There are currently ", length(unfav.markers$Ensembl), " known unfavorable prognostic markers for pancreatic cancer."))

# how many pathways contain prognostic markers
fav.markers.pwy <- fav.markers[fav.markers$Ensembl %in% pwy.genes$V2,]
unfav.markers.pwy <- unfav.markers[unfav.markers$Ensembl %in% pwy.genes$V2,]

cat(paste0("Number of favorable markers in PancCanNet pathways: ",length(fav.markers.pwy$Ensembl),"\n","Number of unfavorable markers in PancCanNet pathways: ",length(unfav.markers.pwy$Ensembl)))

# visualize prognostic markers in pathway gene network
suid.markers <- RCy3::cloneNetwork(suid)
RCy3::setCurrentNetwork(suid.markers)
RCy3::renameNetwork("PancCanNet pathway-gene network: Markers")

# adapt visual style
toggleGraphicsDetails()
RCy3::copyVisualStyle("crosstalk","markers")
RCy3::setVisualStyle("markers")

RCy3::setNodeBorderColorBypass(fav.markers.pwy$Ensembl, new.colors = "#902D10")
RCy3::setNodeBorderWidthBypass(fav.markers.pwy$Ensembl, new.sizes = 10)
RCy3::setNodeBorderColorBypass(unfav.markers.pwy$Ensembl, new.colors = "#052549")
RCy3::setNodeBorderWidthBypass(unfav.markers.pwy$Ensembl, new.sizes = 10)
```

Explore the network in Cytoscape. You will see that most pathways
contain at least a couple of marker genes, both favorable (red border
color) and unfavorable (blue border color).

### 

Next, we will use the knowledge from the pathology information from the
Human Protein Atlas and TCGA to identify and visualize pancreatic
cancer-enriched genes and cancer-enriched genes.

``` r
# genes enriched only in pancreatic cancer
link <- "https://www.proteinatlas.org/api/search_download.php?search=cancer_category_rna%3Apancreatic%20cancer%3BDetected%20in%20single&columns=g,eg&compress=no&format=tsv"
pancreatic.cancer.only <- read.table(file = link, sep = '\t', header = TRUE)

# genes enriched in only a few cancers incl. pancreatic cancer
link <- "https://www.proteinatlas.org/api/search_download.php?search=cancer_category_rna%3Apancreatic%20cancer%3BDetected%20in%20some&columns=g,eg&compress=no&format=tsv"
pancreatic.cancer.some <- read.table(file = link, sep = '\t', header = TRUE)

# genes enriched in many cancers
link <- "https://www.proteinatlas.org/api/search_download.php?search=cancer_category_rna%3Apancreatic%20cancer%3BDetected%20in%20many&columns=g,eg&compress=no&format=tsv"
pancreatic.cancer.many <- read.table(file = link, sep = '\t', header = TRUE)

# genes enriched in all cancers
link <- "https://www.proteinatlas.org/api/search_download.php?search=cancer_category_rna%3Apancreatic%20cancer%3BDetected%20in%20all&columns=g,eg&compress=no&format=tsv"
pancreatic.cancer.all <- read.table(file = link, sep = '\t', header = TRUE)

pancreatic.cancer.only.pwy <- pancreatic.cancer.only[pancreatic.cancer.only$Ensembl %in% pwy.genes$V2,]
pancreatic.cancer.some.pwy <- pancreatic.cancer.some[pancreatic.cancer.some$Ensembl %in% pwy.genes$V2,]
pancreatic.cancer.many.pwy <- pancreatic.cancer.many[pancreatic.cancer.many$Ensembl %in% pwy.genes$V2,]
pancreatic.cancer.all.pwy <- pancreatic.cancer.all[pancreatic.cancer.all$Ensembl %in% pwy.genes$V2,]

cat(paste0("There are ",length(pancreatic.cancer.only$Ensembl), " genes only detected in pancreatic cancer of which ",length(pancreatic.cancer.only.pwy$Ensembl), " are in PancCanNet pathways.\n\n","There are ",length(pancreatic.cancer.some$Ensembl), " genes detected in some cancers incl. pancreatic cancer of which ",length(pancreatic.cancer.some.pwy$Ensembl), " are in PancCanNet pathways.\n\n","There are ",length(pancreatic.cancer.many$Ensembl), " genes detected in many cancers incl. pancreatic cancer of which ",length(pancreatic.cancer.many.pwy$Ensembl), " are in PancCanNet pathways.\n\n","There are ",length(pancreatic.cancer.all$Ensembl), " detected in all cancers of which ",length(pancreatic.cancer.all.pwy$Ensembl), " are in PancCanNet pathways."))

# visualize specificity in pathway gene network
suid.spec <- RCy3::cloneNetwork(suid)
RCy3::setCurrentNetwork(suid.spec)
RCy3::renameNetwork("PancCanNet pathway-gene network: Cancer specificity")

genes <- RCy3::getTableColumns(columns = c("name","type"))

genes <- genes[genes$type=="Gene",]
genes[,"Specificity"] <- "none"
genes[genes$name %in% pancreatic.cancer.only$Ensembl,"Specificity"] <- "only"
genes[genes$name %in% pancreatic.cancer.some$Ensembl,"Specificity"] <- "some"
genes[genes$name %in% pancreatic.cancer.many$Ensembl,"Specificity"] <- "many"
genes[genes$name %in% pancreatic.cancer.all$Ensembl,"Specificity"] <- "all"
RCy3::loadTableData(genes,data.key.column = "name", table.key.column = "name")

# adapt visual style
RCy3::copyVisualStyle("crosstalk","specificity")
RCy3::setVisualStyle("specificity")

RCy3::toggleGraphicsDetails()
RCy3::setNodeColorMapping("Specificity", c('only','some','many','all','none'), c("#0594F4","#17285E","#4D5B88","#737C9D","#000000"), mapping.type = 'd', default.color = "#FFF7BC", style.name = "specificity")
```

Explore the network in Cytoscape. The darker the color the more specific
the expression is pancreatic cancer.

## Visualize your own data

Now let’s load and explore your differential gene or protein expression
dataset and visualize it on the pathway-gene network in Cytoscape.

Example dataset: pancreatic ductal adenocarcinoma vs. normal Downloaded
from Expression Atlas:
<https://www.ebi.ac.uk/gxa/experiments/E-GEOD-15471/>

Ensembl gene id, gene name, log2FC, p-value

``` r
data <- read.delim("data/E-GEOD-15471.tsv")

# overview table > pathway, genes, measured genes, DEGs+, DEGs-
pwy.data <- merge(pwy.genes, data, by.x="V2", by.y="Gene.ID")

print(paste0(length(pwy.data$Gene.Name)," out of ", length(pwy.genes$V1), " genes in PancCanNet pathways are measured in the dataset."))

# pathway-gene network
suid.viz <- RCy3::cloneNetwork(suid)
RCy3::setCurrentNetwork(suid.viz)
RCy3::renameNetwork("PancCanNet pathway-gene network: Data visualization")

RCy3::loadTableData(pwy.data, data.key.column = "V2", table.key.column = "name")

# adapt visual style
RCy3::copyVisualStyle("crosstalk","data")
RCy3::setVisualStyle("data")

RCy3::toggleGraphicsDetails()
node.colors <- c(rev(brewer.pal(3, "RdBu")))
setNodeColorMapping("log2foldchange", c(-1,0,1), node.colors, default.color = "#99FF99", style.name = "data")

RCy3::setNodeLabelMapping(table.column = "Gene.Name", style.name = "data")
```

You can now explore the gene expression changes in the PancCanNet
pathways. (blue = down-regulated, white = not changed, red =
up-regulated)
