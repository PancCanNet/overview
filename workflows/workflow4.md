# Drug-target identification

* Download .Rmd file for the workflow [here](PancCanNet-workflow4.Rmd). 
* Open the .Rmd file in RStudio (see installation guide here) and walk through one section after the other. 

## Introduction

This workflow performs ... After ..., the workflow
will put the results in the context of pancreatic cancer using the
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

```

## Load and explore differential gene expression dataset

In the following section, you will load a differential gene expression
dataset. You can easily replace the example dataset with your own
dataset (copy it in the “data” folder and change the file name below).
In this workflow, we will identify affected pathways in two pancreatic
cancer subtypes, and visualize the data on the pathway models.

``` r
# load dataset

```

Interpretation


## Perform 

We will now perform ....

``` r

```

## Pathway  visualization

Please make sure you have Cytoscape opened before you run this code!

``` r

```

You can investigate the pathway-gene networks in Cytoscape and see which
genes are in multiple pathways and which ones are not present in any of
the altered pathways.

## Visualize your own data

To study the detailed molecular mechanisms, you can also open enriched
pathway models individually in Cytoscape and visualize the gene
expression data on the pathways.

``` r

```

## Save session

```{r}
cys.file <- file.path(getwd(), "output/workflow2.cys")
saveSession(cys.file) 
```
