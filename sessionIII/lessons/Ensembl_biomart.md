---
layout: topic
title: BioMart (Ensembl/R)
author: Mary Piper, Radhika Khetani
date: "Tuesday, May 16th, 2017"
---

Approximate time: 45 minutes

## Learning Objectives

* learn how to use features of the Ensembl biological database and genome browser to access information and data during an NGS analysis

## Ensembl

[*Ensembl*](http://useast.ensembl.org/index.html) provides a website that acts as a **single point of access to annotated genomes** for vertebrate species. 

![ensembl_homepage](../img/ensembl_interface.png)

Ensembl contains extensive genomic information and we can mine this information readily using the "BioMart" tool. You can access it directly using the web interface; alternatively, there is an R package (Bioconductor), "biomaRt", available for mining Ensembl data. 

### Biomart web interface

The Web interface can be accessed from the Ensembl home page.

![ensembl_biomart](../img/ensembl_biomart.png). 

Briefly, BioMart requires three pieces of information that can be selected from drop-down menus and using check boxes on the web interface:

1. **Choose a database to mine.** Options are Ensembl Gene, Ensembl Variation, Ensembl Regulation, and [Vega](https://en.wikipedia.org/wiki/Vertebrate_and_Genome_Annotation_Project) databases. You will be able to choose your species of interest within these databases.
2. **Select the type and content of input/query.** Your query can be genomic region(s), specific gene(s), known variant(s), etc.
3. **Choose the attributes/content to output.** Depending on your query, this can be almost any related genomic information.

If you are interested in exploring how to use the web interface, Ensembl has a [video tutorial](http://www.ensembl.org/Multi/Help/Movie?db=core;id=189) that goes over the various aspects. 

In this section we will focus on the R package to mine genomic information from Ensembl. Please note that the functions within the R package also require the 3 pieces of information listed above.

### biomaRt R package

This is a very convenient way to access and mine the database, since you can easily add any steps with biomart to your DGE workflow within R. 

Let's explore bioMart functionality in R using the row names or the Ensembl IDs from the mm_counts data frame. Out goal here is to **obtain the gene names for a list of Ensembl mouse IDs**, which is something you will find yourself doing relatively frequently as you start working on large genomic datasets. 

Click on the link to the [counts file](https://raw.githubusercontent.com/hbc/NGS_Data_Analysis_Course/master/sessionIV/results/counts.txt) and save it to your `data` folder.

Read in the counts file:

```r
# Read in counts file

example_counts <- read.table("data/counts.txt")

mm_counts <- head(example_counts, n=50)
```

Load the biomaRt library:

```r
# Load library

library("biomaRt")
```

As mentioned before, the same three components required by the web interface are required by the R package. To reiterate, the following steps will demonstrate how to ***obtain the gene names for a list of Ensembl mouse IDs***.

#### **Step 1: Choose a dataset**

Similar to the web interface, we will choose a database (`Ensembl Genes ##`) and a species dataset (`Mus musculus genes (GRCm38.p#)`). 

Choose a BioMart database - we will choose the `Ensembl Genes ##`:

```r
# To connect to a BioMart database - useMart()

listMarts(host =  'www.ensembl.org')

ensembl_genes <- useMart('ENSEMBL_MART_ENSEMBL',
                        host =  'www.ensembl.org')
```

Choose the *Mus musculus* genes (GRCm38.p4) dataset within the `Ensembl Genes ##` database:

```r
# To query the BioMart database for a specific species - useDataset() within the query command

datasets <- listDatasets(ensembl_genes)

View(datasets)

mouse <- useDataset("mmusculus_gene_ensembl",
                   useMart('ENSEMBL_MART_ENSEMBL',
                           host =  'www.ensembl.org'))
```

#### **Step 2: Select your filters or inputs**

We can build a query of our dataset using the `getBM()` function and specifying the filters, filter values, attributes, and mart object database to search: `getBM(filters, values, attributes, mart)`

First we can specify our input using the `filters`argument. 

**What is our input?** We want to return gene names for a list of Ensembl mouse IDs from within our `mm_counts` dataframe; therefore our input will be Ensembl IDs and their values will be the row names of our mm_counts dataframe.

```r
# To build a query - getBM(filters, values, ...)

## "Filters" is a vector for the type of input; in our case, our input is Ensembl IDs

filters <- listFilters(mouse)

View(filters)

getBM(filters= "ensembl_gene_id", ...)  # The "..." represents that the getBM() function is not complete

                    
## "Values" is a vector of values for the filter; in our case, our Ensembl IDs are the row names of the mm_counts dataset

getBM(filters= "ensembl_gene_id", 
		values= row.names(mm_counts), ...)
```

#### **Step 3: Choose the attributes to output**

We can continue building our `getBM()` function by specifying what we want output for each of our Ensembl IDs using the `attributes` argument. We would like to output the Ensembl ID and the gene name. 

```r
# To build a query - getBM(filters, values, attributes, ...)

## "Attributes" is a vector of attributes for the output we want to generate

attributes <- listAttributes(mouse)

View(attributes)

## Use BioMart to return gene names for a list of Ensembl IDs:

gene_names <- getBM(filters= "ensembl_gene_id",
                    values= row.names(mm_counts), 
                    attributes= c("ensembl_gene_id", "external_gene_name"), ...)
```

Finally, to complete the `getBM()` function, we need to specify which dataset to query.

```r
# To build a query - getBM(filters, values, attributes, mart)

gene_names <- getBM(filters= "ensembl_gene_id",
                    values= row.names(mm_counts), 
                    attributes= c("ensembl_gene_id", "external_gene_name"), 
                    mart= mouse)

## Now we can run the query. BioMart queries can take a bit of time depending on the size of your dataset and the attributes you are asking for.                    

View(gene_names)
                    
```

Now that we have our gene names, we need to match them to the Ensembl IDs in our mm_counts dataset. If the columns from two dataframes have the same name, we can merge the dataframes using those columns:

```r
# Merge the two dataframes by ensembl_gene_id

ensembl_results <- merge(mm_counts, gene_names, by.x="row.names", by.y = "ensembl_gene_id")

write.csv(ensembl_results, "results/annotated_mm_counts.csv", quote=F)
```
#### What if you are using an older genome? 

Check the archived BioMart sites to determine the archived database desired. 

If we want to use the archived databases in R, we need to change our query a bit to specify an older genome build:

```r
# Using an older genome build

mouse_mm9 <- useDataset("mmusculus_gene_ensembl",
                   useMart('ENSEMBL_MART_ENSEMBL',
                           host =  'may2012.archive.ensembl.org'))
```

The filters and attributes change for different builds of the genome, so you might find yourself looking them up if you change builds. Using our previous filters and attributes we would have run the following:                    

```r
# DO NOT RUN			   
gene.names_mm9 <- getBM(filters= "ensembl_gene_id", 
                    attributes= c("ensembl_gene_id", "external_gene_name"),
                    values= row.names(mm_counts),
                    mart= mouse_mm9)
```

However, for `attributes`, specifying gene names using `external_gene_name` no longer works. To find the correct way to specify gene names for this mouse build, you will need to look up the attributes again:

```r
attributes_mm9 <- listAttributes(mouse_mm9)

View(attributes_mm9)

gene.names_mm9 <- getBM(filters= "ensembl_gene_id", 
                    attributes= c("ensembl_gene_id", "external_gene_id"),
                    values= row.names(mm_counts),
                    mart= mouse_mm9)

# human archive for GRCH37 genome build: host = "grch37.ensembl.org"
```

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
