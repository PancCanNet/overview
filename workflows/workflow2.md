# Pathway analysis of subtypes of pancreatic cancer (basal / classical vs. healthy)

* Download .Rmd file for the workflow [here](PancCanNet-workflow2.Rmd). 
* Open the .Rmd file in RStudio (see installation guide here) and walk through one section after the other. 

## Introduction

This workflow performs pathway enrichment analysis for differential gene
expression datasets. After a general enrichment anlaysis, the workflow
will put the results in the contect of pancreatic cancer using the
information from the PancCanNet pathway resource.

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
if(!"rWikiPathways" %in% installed.packages()) BiocManager::install
if(!"RCy3" %in% installed.packages()) BiocManager::install("RCy3")
if(!"EnhancedVolcano" %in% installed.packages()) BiocManager::install("EnhancedVolcano")
if(!"VennDiagram" %in% installed.packages()) BiocManager::install("VennDiagram")
if(!"clusterProfiler" %in% installed.packages()) BiocManager::install("clusterProfiler")
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr")
if(!"RColorBrewer" %in% installed.packages()) BiocManager::install("RColorBrewer")


# load required libraries
library(rstudioapi)
library(rWikiPathways)
library(RCy3)
library(EnhancedVolcano)
library(VennDiagram)
library(clusterProfiler)
library(dplyr)
library(RColorBrewer)

# set working environment to source file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Load and explore differential gene expression dataset

In the following section, you will load a differential gene expression
dataset. You can easily replace the example dataset with your own
dataset (copy it in the “data” folder and change the file name below).
In this workflow, we will identify affected pathways in two pancreatic
cancer subtypes, and visualize the data on the pathway models.

``` r
# load dataset
dataset <- read.delim("data/GSE71729-dataset.txt")

# filter genes without Entrez Gene identifier
data.panc <- dataset %>% tidyr::drop_na(Entrez.Gene)
colnames(data.panc)[2] <- "GeneName"
colnames(data.panc)[1] <- "GeneId"

EnhancedVolcano(subset(data.panc, select=c(1:7)), title = "Classic subtype", lab = data.panc$GeneName, labSize = 3, x = 'C_logFC', y = 'C_P.Value', pCutoff = 0.05, FCcutoff = 0.585)

png('output/w2-fig1.png')
EnhancedVolcano(subset(data.panc, select=c(1:7)), title = "Classic subtype", lab = data.panc$GeneName, labSize = 3, x = 'C_logFC', y = 'C_P.Value', pCutoff = 0.05, FCcutoff = 0.585)

EnhancedVolcano(subset(data.panc, select=c(1,2,8:12)), title = "Basal subtype", lab = data.panc$GeneName, labSize = 3, x = 'B_logFC', y = 'B_P.Value', pCutoff = 0.05, FCcutoff = 0.585)

png('output/w2-fig2.png')
EnhancedVolcano(subset(data.panc, select=c(1,2,8:12)), title = "Basal subtype", lab = data.panc$GeneName, labSize = 3, x = 'B_logFC', y = 'B_P.Value', pCutoff = 0.05, FCcutoff = 0.585)

deg.basal <- unique(data.panc[!is.na(data.panc$B_P.Value) & data.panc$B_P.Value < 0.05 & abs(data.panc$B_logFC) > 0.58,c(1,2)])
basal.up <- unique(data.panc[!is.na(data.panc$B_P.Value) & data.panc$B_P.Value < 0.05 & data.panc$B_logFC > 0.58,c(1,2)])
basal.down <- unique(data.panc[!is.na(data.panc$B_P.Value) & data.panc$B_P.Value < 0.05 & data.panc$B_logFC < -0.58,c(1,2)])

deg.classical <- unique(data.panc[!is.na(data.panc$C_P.Value) & data.panc$C_P.Value < 0.05 & abs(data.panc$C_logFC) > 0.58,c(1,2)])
classical.up <- unique(data.panc[!is.na(data.panc$C_P.Value) & data.panc$C_P.Value < 0.05 & data.panc$C_logFC > 0.58,c(1,2)])
classical.down <- unique(data.panc[!is.na(data.panc$C_P.Value) & data.panc$C_P.Value < 0.05 & data.panc$C_logFC < -0.58,c(1,2)])

venn.diagram(x = list(basal.up$GeneId, basal.down$GeneId, classical.up$GeneId, classical.down$GeneId),
  category.names = c("Basal up", "Basal down" ,"Classical up", "Classical down"),
  filename = 'output/w2-fig3.png',
  output=FALSE,
  col=c("#440154ff","#440154ff", '#21908dff','#21908dff'),
  fill = c(alpha("#440154ff",0.3),alpha("#440154ff",0.3), alpha('#21908dff',0.3),alpha('#21908dff',0.3)),
  cex = 1.5,
)

bkgd.genes <- unique(data.panc[,c(1,2)])
```

You can use the Volcano plots, venn diagrams and gene lists to see which
genes are up- and down-regulated in the two subtypes.

## Load the pathway collection of interest

There are different online pathway databases that can be used or even
combined in this step (e.g. WikiPathways, Reactome, KEGG). In the
enrichment analysis, it is easy to exchange/combine the resources and
some of the visualizations work as well. The actual pathway
visualization and PancCanNet analysis will only work with WikiPathways.

``` r
# load current GMT + merge with all PancCanNet pathways
gmt <- rWikiPathways::downloadPathwayArchive(organism = "Homo sapiens",format = "gmt", date = "20211210")
wp2gene <- readPathwayGMT(gmt)
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE

pwys <- rWikiPathways::getPathwayIdsByCurationTag("Curation:PancCanNet")

# retrieve gene list per pathway
mylist <- list()
x <- 1
for(p in pwys) {
  genes <- rWikiPathways::getXrefList(p,"L")
  for(g in genes) {
    mylist[[x]] <- c(p,g)
    x <- x+1
  }
}
pwy.genes <- as.data.frame(do.call("rbind",mylist))
colnames(pwy.genes) <- c("wpid","gene")
pwy.list <- unique(rbind(pwy.genes, wpid2gene))

rm(mylist,x,genes,p,g,wp2gene)

# output number of pathways
print(paste0("Number of pathways included in the collection: ", length(unique(pwy.list$wpid))))
```

## Perform enrichment analysis

We will now perform pathway enrichment with the gene sets loaded before.
The clusterProfiler R-package is used to perform overrepresentation
analysis (ORA). The function can be easily replaced to use other
enrichment methods (GSEA / rSEA / etc). We will run the analysis
separately for basal and classical subtype.

``` r
## BASAL SUBTYPE

ewp.basal <- clusterProfiler::enricher(
  deg.basal$GeneId,
  universe = bkgd.genes$GeneId,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.02,
  TERM2GENE = wpid2gene)

ewp.basal.res <- as.data.frame(ewp.basal) 

# number of genes measured in pathways
length(ewp.basal@universe)

# number of DEG in pathways
length(deg.basal$GeneId[deg.basal$GeneId %in% unique(wpid2gene$gene)])

num.pathways.basal <- dim(ewp.basal.res)[1]

# export enrichment result
png('output/w2-fig4.png', width = 1200, height=1000)
ggplot(ewp.basal[1:num.pathways.basal], aes(x=reorder(Description, -pvalue), y=Count)) +
  geom_bar(stat ="identity", fill="#BA8CD7") +
  coord_flip() +
  labs(x="", y="Basal DEG gene count", fill="") +
  theme(axis.text=element_text(size=15)) + 
  theme(legend.position="none")
dev.off()

ggplot(ewp.basal[1:num.pathways.basal], aes(x=reorder(Description, -pvalue), y=Count)) +
  geom_bar(stat ="identity", fill="#BA8CD7") +
  coord_flip() +
  labs(x="", y="Basal DEG gene count", fill="") +
  theme(axis.text=element_text(size=15)) + 
  theme(legend.position="none")

# TODO: highlight panccannet pathways!

## CLASSICAL SUBTYPE

ewp.classical <- clusterProfiler::enricher(
  deg.classical$GeneId,
  universe = bkgd.genes$GeneId,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  TERM2GENE = wpid2gene)

ewp.classical.res <- as.data.frame(ewp.classical) 

# number of genes measured in pathways
length(ewp.classical@universe)

# number of DEG in pathways
length(deg.classical$GeneId[deg.classical$GeneId %in% unique(wpid2gene$gene)])

num.pathways.classical <- dim(ewp.classical.res)[1]

# export enrichment result
ggplot(ewp.classical[1:num.pathways.classical], aes(x=reorder(Description, -pvalue), y=Count)) +
  geom_bar(stat ="identity", fill="#BA8CD7") +
  coord_flip() +
  labs(x="", y="Classical DEG gene count", fill="") +
  theme(axis.text=element_text(size=15)) + 
  theme(legend.position="none")

png('output/w2-fig5.png', width = 1200, height=1000)
ggplot(ewp.classical[1:num.pathways.classical], aes(x=reorder(Description, -pvalue), y=Count)) +
  geom_bar(stat ="identity", fill="#BA8CD7") +
  coord_flip() +
  labs(x="", y="Classical DEG gene count", fill="") +
  theme(axis.text=element_text(size=15)) + 
  theme(legend.position="none")
dev.off()

venn.diagram(x = list(ewp.basal.res$ID, ewp.classical.res$ID),
  category.names = c("Basal" , "Classical"),
  filename = 'output/w2-fig6.png',
  output=TRUE,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 1.5,
)
```

## Pathway overlap visualization

There is often crosstalk and overlap between pathways enriched in gene
expression analyses. The following step visualizes the overlap between
the enriched pathways in a pathway-gene network.

The genes not present in any pathway are included in the visualization
but can be removed if that is preferred.

Please make sure you have Cytoscape opened before you run this code!

``` r
pwy <- unique(ewp.basal.res[,c(1,2)])
colnames(pwy) <- c("id","label")
pwy$type <- 'pathway'

edges <- wpid2gene[wpid2gene$wpid %in% pwy$id,]
colnames(edges) <- c("source", "target")

genes <- unique(deg.basal)
colnames(genes) <- c("id","label")
genes$type <- 'gene'

edges <- unique(edges[edges$target %in% genes$id,])

nodes <- rbind(pwy, genes)
rownames(nodes) <- NULL

RCy3::createNetworkFromDataFrames(nodes=nodes,edges=edges,title="Basal Subtype", collection="Pathway-Gene-Associations Basal")

loadTableData(data.panc, data.key.column = "GeneId", table.key.column = "id")

# Visual style
RCy3::copyVisualStyle("default","wp.vis")

RCy3::setNodeLabelMapping("label", style.name="wp.vis")

RCy3::lockNodeDimensions(TRUE, style.name="wp.vis")
RCy3::setNodeShapeMapping('type', c('gene','pathway'), c("ellipse","hexagon"), style.name="wp.vis")

RCy3::setNodeSizeMapping('type', c('gene','pathway'), c(40,25), mapping.type = "d", style.name = "wp.vis")

data.values<-c(-1,0,1) 
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("B_logFC", data.values, node.colors, default.color = "#99FF99", style.name = "wp.vis")

RCy3::setVisualStyle("wp.vis")

RCy3::toggleGraphicsDetails()

# Saving output
png.file <- file.path(getwd(), "output/w2-fig7.png")
exportImage(png.file,'PNG', zoom = 500)

## CLASSICAL SUBTYPE

pwy <- unique(ewp.classical.res[,c(1,2)])
colnames(pwy) <- c("id","label")
pwy$type <- 'pathway'

edges <- wpid2gene[wpid2gene$wpid %in% pwy$id,]
colnames(edges) <- c("source", "target")

genes <- unique(deg.classical)
colnames(genes) <- c("id","label")
genes$type <- 'gene'

edges <- unique(edges[edges$target %in% genes$id,])

nodes <- rbind(pwy, genes)
rownames(nodes) <- NULL

RCy3::createNetworkFromDataFrames(nodes=nodes,edges=edges,title="Classical Subtype", collection="Pathway-Gene-Associations Classical")

loadTableData(data.panc, data.key.column = "GeneId", table.key.column = "id")

# Visual style
RCy3::copyVisualStyle("default","wp.vis.classical")

RCy3::setNodeLabelMapping("label", style.name="wp.vis.classical")

RCy3::lockNodeDimensions(TRUE, style.name="wp.vis.classical")
RCy3::setNodeShapeMapping('type', c('gene','pathway'), c("ellipse","hexagon"), style.name="wp.vis.classical")

RCy3::setNodeSizeMapping('type', c('gene','pathway'), c(40,25), mapping.type = "d", style.name = "wp.vis.classical")

data.values<-c(-1,0,1) 
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("C_logFC", data.values, node.colors, default.color = "#99FF99", style.name = "wp.vis.classical")

RCy3::setVisualStyle("wp.vis.classical")

RCy3::toggleGraphicsDetails()

# Saving output
png.file <- file.path(getwd(), "output/w2-fig8.png")
exportImage(png.file,'PNG', zoom = 500)
```

You can investigate the pathway-gene networks in Cytoscape and see which
genes are in multiple pathways and which ones are not present in any of
the altered pathways.

## Visualize your own data

To study the detailed molecular mechanisms, you can also open enriched
pathway models individually in Cytoscape and visualize the gene
expression data on the pathways.

``` r
# for all enriched pathways
enrich.pwy <- rbind(ewp.basal.res,ewp.classical.res)
for(p in unique(enrich.pwy$ID)) {
  RCy3::commandsRun(paste0("wikipathways import-as-pathway id=",p))
  toggleGraphicsDetails()
  loadTableData(table = "node", data = data.panc, data.key.column = "GeneName", table.key.column = "name")

  # apply visual style 
  RCy3::setNodeCustomHeatMapChart(c("B_logFC","C_logFC"), slot = 2, style.name = "WikiPathways")
  
  RCy3::setVisualStyle("WikiPathways")

  png.file <- file.path(getwd(), paste0("output/w2-",p,".png"))
  exportImage(png.file,'PNG', zoom = 500)
}
```

## Save session

```{r}
cys.file <- file.path(getwd(), "output/workflow2.cys")
saveSession(cys.file) 
```
