# Analysing a pancreatic cancer immune panel dataset (PPI network + pathway information)

* Download .Rmd file for the workflow [here](PancCanNet-workflow3.Rmd). 
* Open the .Rmd file in RStudio (see installation guide here) and walk through one section after the other. 

## Introduction

This workflow builds a protein-protein interaction network for all altered genes in the immune panel dataset. Afterwards, the PPI network is extended with pathway information from WikiPathway to identify relevant pathways and missing pathway information (affected genes not present in any pathway).

The PancCanNet project is a collaboration between the Erasmus MC,
Maastricht University and Omnigen. The pathway portal is integrated in
the collaborative pathway database WikiPathways
(<http://panccannet.wikipathways.org>) and all analysis workflows are
available on Github (<https://panccannet.github.io>).

Pathway curation is an ongoing effort and we are continuously improving
current pathway models and adding new pathway models. The workflow will
always use the most-up-to-date information and results can slightly
differ to the information on the project website.
