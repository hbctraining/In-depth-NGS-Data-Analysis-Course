---
title: "Data manipulation"
authors: Meeta Mistry, Mary Piper
date: "Tuesday, May 9, 2017"
---
Approximate time: 60 min

## Learning Objectives

* Demonstrate how to subset, merge, and create new datasets from existing data structures in R.
* Export data tables and plots for use outside of the R environment.

### Dataframes

Dataframes (and matrices) have 2 dimensions (rows and columns), so if we want to select some specific data from it we need to specify the "coordinates" we want from it. We use the same square bracket notation but rather than providing a single index, there are *two indexes required*. Within the square bracket, **row numbers come first followed by column numbers (and the two are separated by a comma)**. Let's explore the `metadata` dataframe, shown below are the first six samples:

![metadata](../img/metadata.png)

For example:

```r
metadata[1, 1]   # element from the first row in the first column of the data frame
metadata[1, 3]   # element from the first row in the 3rd column
```

Now if you only wanted to select based on rows, you would provide the index for the rows and leave the columns index blank. The key here is to include the comma, to let R know that you are accessing a 2-dimensional data structure:

```r
metadata[3, ]    # vector containing all elements in the 3rd row
```

If you were selecting specific columns from the data frame - the rows are left blank:

```r
metadata[ , 3]    # vector containing all elements in the 3rd column
```

Just like with vectors, you can select multiple rows and columns at a time. Within the square brackets, you need to provide a vector of the desired values:	

```r
metadata[ , 1:2] # dataframe containing first two columns
metadata[c(1,3,6), ] # dataframe containing first, third and sixth rows
```

For larger datasets, it can be tricky to remember the column number that corresponds to a particular variable. (Is celltype in column 1
or 2? oh, right... they are in column 1). In some cases, the column number for a variable can change if the script you are using adds or removes columns. It's therefore often better to use column names to refer to a particular variable, and it makes your code easier to read and your intentions clearer.

```r
metadata[1:3 , "celltype"] # elements of the celltype column corresponding to the first three samples
```


You can do operations on a particular column, by selecting it using the `$` sign. In this case, the entire column is a vector. For instance, to extract all the genotypes from our dataset, we can use: 

```r
metadata$genotype 
```
You can use `colnames(metadata)` or `names(metadata)` to remind yourself of the column names. We can then supply index values to select specific values from that vector. For example, if we wanted the genotype information for the first five samples in `metadata`:

```r
colnames(metadata)

metadata$genotype[1:5]
```

The `$` allows you to select a single column by name. To select multiple columns by name, you need to  concatenate a vector of strings that correspond to column names: 

```r
metadata[, c("genotype", "celltype")]
```

```r
          genotype celltype
sample1        Wt    typeA
sample2        Wt    typeA
sample3        Wt    typeA
sample4        KO    typeA
sample5        KO    typeA
sample6        KO    typeA
sample7        Wt    typeB
sample8        Wt    typeB
sample9        Wt    typeB
sample10       KO    typeB
sample11       KO    typeB
sample12       KO    typeB
```

While there is no equivalent `$` syntax to select a row by name, you can select specific rows using the row names. To remember the names of the rows, you can use the `rownames()` function:

```r
rownames(metadata)

metadata[c("sample10", "sample12"),]
```

#### Selecting using indexes with logical operators

With dataframes, similar to vectors, we can use logical vectors for specific columns in the dataframe to select only the rows in a dataframe with TRUE values at the same position or index as in the logical vector. We can then use the logical vector to return all of the rows in a dataframe where those values are TRUE.

```r
idx <- metadata$celltype == "typeA"
	
metadata[idx, ]
```

##### Selecting indexes with logical operators using the `which()` function
As you might have guessed, we can also use the `which()` function to return the indexes for which the logical expression is TRUE. For example, we can find the indexes where the `celltype` is `typeA` within the `metadata` dataframe:

```r
idx <- which(metadata$celltype == "typeA")
	
metadata[idx, ]
```

Or we could find the indexes for the metadata replicates 2 and 3:

```r
idx <- which(metadata$replicate > 1)
	
metadata[idx, ]
```

***

**Exercise**  

Subset the `metadata` dataframe to return only the rows of data with a genotype of `KO`.
	
***

#### Subsetting dataframes using logical operators and the `subset()` function

Another way of partitioning **dataframes** is using the `subset()` function to return the rows of the dataframe for which the logical expression is TRUE. Allowing us to the subset the data in a single step. The syntax for the `subset()` function is:

```r
subset(dataframe, column_name == "value") # Any logical expression could replace the `== "value"`
```

For example, we can look at the samples of a specific celltype "typeA":

```r
subset(metadata, celltype == "typeA")
```

```r
         genotype celltype replicate
sample1       Wt    typeA         1
sample2       Wt    typeA         2
sample3       Wt    typeA         3
sample4       KO    typeA         1
sample5       KO    typeA         2
sample6       KO    typeA         3
```

We can also use the `subset` function with the other logical operators in R. For example, suppose we wanted to subset to keep only the **Wt samples** from the **typeA** celltype.

```r
subset(metadata, celltype == "typeA" & genotype == "Wt")
```

```r
        genotype celltype replicate
sample1       Wt    typeA         1
sample2       Wt    typeA         2
sample3       Wt    typeA         3
```

Alternatively, we could try looking at only the first two replicates of each sample set. Here, we can use the less than operator since replicate is currently a numeric vector. Adding in the argument `select` allows us to specify which columns to keep, with the syntax:

```r
subset(dataframe, column_name == "value", select = "name of column(s) to return")
```

Which columns are left?

```r
sub_meta <- subset(metadata, replicate < 3, select = c('genotype', 'celltype'))
```

***

**Exercise** 

1. Return only the rows of data from the `metadata` dataframe with `genotype` of `Wt`. 

2. Return only the `celltype` information for those samples from `metadata` dataframe with genotype `KO`.

***

### Writing to file 

Everything we have done so far has only modified the data in R; the files have remained unchanged. Whenever we want to save our datasets to file, we need to use a `write` function in R. 

To write our matrix to file in comma separated format (.csv), we can use the `write.csv` function. There are two required arguments: the variable name of the data structure you are exporting, and the path and filename that you are exporting to. By default the delimiter is set, and columns will be separated by a comma:

```r
write.csv(sub_meta, file="data/subset_meta.csv")
```

Similar to reading in data, there are a wide variety of functions available allowing you to export data in specific formats. Another commonly used function is `write.table`, which allows you to specify the delimiter you wish to use. This function is commonly used to create tab-delimited files.

> **NOTE**:
>  
> Sometimes when writing a dataframe with row names to file, the column names will align starting with the row names column. To avoid this, you can include the argument `col.names = NA` when writing to file to ensure all of the column names line up with the correct column values.

> Writing a vector of values to file requires a different function than the functions available for writing dataframes. You can use `write()` to save a vector of values to file. For example:
>

```r
write(glengths, file="genome_lengths.txt", ncolumns=1)
```
>

***

> ### An R package for data manipulation
> The methods presented above are using base R functions for data manipulation. For more advanced R users, 
> the package `dplyr` is a fairly new (2014) package that tries to provide easy
> tools for the most common data manipulation tasks. It is built to work directly
> with data frames. The thinking behind it was largely inspired by the package
> `plyr` which has been in use for some time but suffered from being slow in some
> cases.` dplyr` addresses this by porting much of the computation to C++. An
> additional feature is the ability to work with data stored directly in an
> external database. The benefits of doing this are that the data can be managed
> natively in a relational database, queries can be conducted on that database, and only the results of the query returned.
>
> `dplyr` is part of the [tidyverse suite of packages](http://tidyverse.org), which were designed to work together to quickly and easily perform the most common operations needed for data science.









---

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson is adapted from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
