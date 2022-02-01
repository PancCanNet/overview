# Drug-target identification

* Download .Rmd file for the workflow [here](PancCanNet-workflow4.Rmd). 
* Open the .Rmd file in RStudio (see installation guide here) and walk through one section after the other. 

## Introduction

In this workflow, we will visualize transcriptomics data on the Pancreatic adenocarcinoma pathway (WP4263) network. We will extend this network with drug-target interaction information from DrugBank and visualize the results. We will then build a protein-protein interaction network using String Disease query (“pancreatic cancer”) and visualize transcriptomics data. We will again extend the PPI network with drug-target interaction information from DrugBank and visualize the results.

The PancCanNet project is a collaboration between the Erasmus MC, Maastricht University and Omnigen. The pathway portal is integrated in the collaborative pathway database WikiPathways (http://panccannet.wikipathways.org) and all analysis workflows are available on Github (https://panccannet.github.io). 

Pathway curation is an ongoing effort and we are continuously improving current pathway models and adding new pathway models. The workflow will always use the most-up-to-date information and results can slightly differ to the information on the project website. 

## R environment setup

First, we need to make sure all required R-packages are installed. Run the following code snippet to install and load the libraries. 
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr")
if(!"rWikiPathways" %in% installed.packages()) BiocManager::install("rWikiPathways")
if(!"RCy3" %in% installed.packages()) BiocManager::install("RCy3")
if(!"RColorBrewer" %in% installed.packages()) BiocManager::install("RColorBrewer")
if(!"rstudioapi" %in% installed.packages()) BiocManager::install("rstudioapi")

library(dplyr)
library(rWikiPathways)
library(RCy3)
library(RColorBrewer)
library(rstudioapi) 

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
```

## Load differential gene expression dataset
In the following section, you will load a differential gene expression dataset. You can easily replace the example dataset with your own dataset (copy it in the "data" folder and change the file name below). 

We will visualize the data on the pathway models. Note: the same dataset is used in workflow 2 and 3.
```{r}
# load dataset 
dataset <- read.delim("data/GSE71729-dataset.txt")

# filter genes without Entrez Gene identifier
data.panc <- dataset %>% tidyr::drop_na(Entrez.Gene)
colnames(data.panc)[2] <- "GeneName"
colnames(data.panc)[1] <- "GeneId"

```

For the rest of this workflow, Cytoscape must be open and running. 

```{r}
RCy3::cytoscapePing()

```

## Visualization of Pancreatic adenocarcinoma pathway (WP4263).

Pancreatic ductal adenocarcinoma (PDA) is the most common malignancy of the pancreas. We will load the Pancreatic adenocarcinoma pathway (WP4263) as a network and visualize transcriptomics data on it.

```{r}

## Open Pancreatic adenocarcinoma pathway as a network (WP4263)
RCy3::commandsRun('wikipathways import-as-network id=WP4263') 
toggleGraphicsDetails()
  
## Visualize of transcriptomics data on pathway 

# load data 
loadTableData(table = "node", data = data.panc, data.key.column = "GeneName", table.key.column = "name")

# node colors and shapes
node.types.a = c("GeneProduct","Metabolite","Pathway","Group","Label") 
node.colors.a <- c("#CCCCCC","#CCCCCC","#99FF99","#CCCCCC","#FFFFFF")
node.shapes.a <- c("rectangle","ellipse","rectangle","octagon","rectangle")

# apply visual style 
RCy3::copyVisualStyle("default","wp.vis1")
RCy3::setNodeLabelMapping("name", style.name="wp.vis1")
RCy3::lockNodeDimensions(TRUE, style.name="wp.vis1")
RCy3::setNodeShapeMapping('Type', node.types.a, node.shapes.a, style.name="wp.vis1")
RCy3::setNodeColorMapping("Type", node.types.a, node.colors.a, mapping.type = 'd', default.color = "#CCCCCC", style.name = "wp.vis1")
RCy3::setNodeBorderColorMapping("Type", c("GeneProduct","Metabolite"), c("#000000","#0000FF"), mapping.type = 'd', default.color = "#CCCCCC", style.name = "wp.vis1")
RCy3::setNodeBorderWidthMapping("Type", c("GeneProduct","Metabolite"), c(5,5), mapping.type = 'd', style.name = "wp.vis1")
RCy3::setNodeCustomHeatMapChart(c("B_logFC","C_logFC"), slot = 2, style.name = "wp.vis1",  colors = c("#CC3300","#FFFFFF","#6699FF","#CCCCCC00"))
RCy3::setVisualStyle("wp.vis1")
```

You can now explore the gene expression changes in the WP4263 network. (blue = down-regulated, white = not changed, red = up-regulated, gray = no data).

The following visual style is applied for the network. 
- GeneProduct: rectangle with thick black border (blue = down-regulated, white = not changed, red = up-regulated, gray = no data).
- Metabolite: gray ellipse with thick blue border
- Pathway: green rectangle
- Group: gray octagon
- Label: white rectangle without border


## Drug-target extension and visualization of Pancreatic adenocarcinoma pathway (WP4263).

We then extend the Pancreatic adenocarcinoma pathway (WP4263) network with Drugbank information using CyTargetLinker and visualize it. 

```{r}
## Extend pathway with drug-target interaction information from DrugBank
drugbank <- file.path(getwd(), "data/drugbank-5.1.7.xgmml")

# run CyTargetLinker
commandsRun(paste0('cytargetlinker extend idAttribute="name" linkSetFiles="', drugbank, '"'))
commandsRun('cytargetlinker applyLayout network="current"')

drug <- read.delim2("data/drugNodes-5.1.7.txt", stringsAsFactors=TRUE)
loadTableData(drug, data.key.column="Identifier", table.key.column="name")

## Improve visualization

# apply visual style
RCy3::setVisualStyle("wp.vis1")
my.drugs <- selectNodes("drug", by.col = "CTL.Type", preserve = FALSE)$nodes 
clearSelection()
setNodeColorBypass(my.drugs, "#DD99FF")
setNodeShapeBypass(my.drugs, "hexagon")
# set drugs label 
drug.labels <- getTableColumns(columns=c("SUID","label"))
drug.labels <- na.omit(drug.labels)
invisible(mapply(function(x,y) setNodeLabelBypass(x,y), drug.labels$SUID, drug.labels$label))

```

Drugs are visualized as purple hexagons. 


## Visualization of PPI network from String Disease query (“pancreatic cancer”)

Next, we will create and visualize a protein-protein interaction (PPI) network from String Disease query "pancreatic cancer" using RCy3 and StringApp. 

```{r}
## PPI network from String Disease query (“pancreatic cancer”) using RCy3 + stringApp
installApp('stringApp') 

# network will be opened in Cytoscape (this might take a while)
commandsRun(paste0('string disease query cutoff=0.9 disease=pancreatic cancer'))

# load data 
loadTableData(table = "node", data = data.panc, data.key.column = "GeneName", table.key.column = "display name")

## Improve visualization

## apply visual style 
RCy3::copyVisualStyle("default","wp.vis2")
RCy3::setNodeLabelMapping("display name", style.name="wp.vis2")
RCy3::lockNodeDimensions(TRUE, style.name="wp.vis2")
RCy3::setNodeColorDefault("#CCCCCC", style.name = "wp.vis2")
RCy3::setNodeCustomHeatMapChart(c("B_logFC","C_logFC"), slot = 2, style.name = "wp.vis2",  colors = c("#CC3300","#FFFFFF","#6699FF","#CCCCCC00"))
RCy3::setVisualStyle("wp.vis2")

```

You can now explore the gene expression changes in the PPI network. (blue = down-regulated, white = not changed, red = up-regulated).

The following visual style is applied for the network. 
- Protein: rectangle, fill color (blue = down-regulated, white = not changed, red = up-regulated, gray = no data).


## Drug-target extension and visualization of PPI network from String Disease query (“pancreatic cancer”)

We extend the PPI network with Drugbank information using CyTargetLinker and visualize it. 

```{r}

## Extend pathway with drug-target interaction information from DrugBank
commandsRun(paste0('cytargetlinker extend idAttribute="display name" linkSetFiles="', drugbank, '"'))
commandsRun('cytargetlinker applyLayout network="current"')

loadTableData(drug, data.key.column="Identifier", table.key.column="display name")
loadTableData(data.panc, data.key.column = "GeneId", table.key.column = "GeneId") #needed? 

## Improve visualization

## apply visual style
RCy3::setVisualStyle("wp.vis2")
my.drugs <- selectNodes("drug", by.col = "CTL.Type", preserve = FALSE)$nodes 
clearSelection()
setNodeColorBypass(my.drugs, "#DD99FF")
setNodeShapeBypass(my.drugs, "hexagon")
# set drugs label 
drug.labels <- getTableColumns(columns=c("SUID","label"))
drug.labels <- na.omit(drug.labels)
invisible(mapply(function(x,y) setNodeLabelBypass(x,y), drug.labels$SUID, drug.labels$label))

```

Drugs are visualized as purple hexagons. 


## Save session

```{r session}
cys.file <- file.path(getwd(), "output/workflow4.cys")
saveSession(cys.file) 
```

## Clear environment

```{r}
rm(list=ls())
```
