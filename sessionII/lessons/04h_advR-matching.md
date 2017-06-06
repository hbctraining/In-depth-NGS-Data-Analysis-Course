---
title: "Advanced R, matching"
authors: Meeta Mistry, Mary Piper
date: "Tuesday, May 9, 2017"
---
Approximate time: 110 min


## Learning Objectives

* Implement matching and re-ordering data within data structures.

## Matching data 

Often when working with genomic data, we have a data file that corresponds with our metadata file. The data file contains measurements from the biological assay for each individual sample. In our case, the biological assay is gene expression and data was generated using RNA-Seq. 

[Download counts data here by right clicking and saving to the `data` folder.](https://raw.githubusercontent.com/hbctraining/Intro-to-R/master/data/counts.rpkm)

Let's read in our expression data (RPKM matrix) that we downloaded previously:

```r
rpkm_data <- read.csv("data/counts.rpkm.txt")
```

Take a look at the first few lines of the data matrix to see what's in there.

```r
head(rpkm_data)
```

It looks as if the sample names (header) in our data matrix are similar to the row names of our metadata file, but it's hard to tell since they are not in the same order. We can do a quick check of the number of columns in the count data and the rows in the metadata and at least see if the numbers match up. 

```r
ncol(rpkm_data)
nrow(metadata)
```

What we want to know is, **do we have data for every sample that we have metadata?** 

## The `%in%` operator
 
Although lacking in [documentation](http://dr-k-lo.blogspot.com/2013/11/) this operator is well-used and convenient once you get the hang of it. The operator is used with the following syntax: 

```r
vector1_of_values %in% vector2_of_values
```

It will take a vector as input to the left and will **evaluate each element to see if there is a match in the vector that follows on the right of the operator.** *The two vectors do not have to be the same size.* This operation will return a vector of the same length as vector1 containing logical values to indicate whether or not there was a match. Take a look at the example below:

```r
A <- c(1,3,5,7,9,11)   # odd numbers
B <- c(2,4,6,8,10,12)  # even numbers

# test to see if each of the elements of A is in B	
A %in% B
```

```r
## [1] FALSE FALSE FALSE FALSE FALSE FALSE
```

Since vector A contains only odd numbers and vector B contains only even numbers, there is no overlap and so the vector returned contains a `FALSE` for each element. Let's change a couple of numbers inside vector B to match vector A:


```r
A <- c(1,3,5,7,9,11)   # odd numbers
B <- c(2,4,6,8,1,5)  # add some odd numbers in 
```

```r
# test to see if each of the elements of A is in B
A %in% B
```

```r
## [1]  TRUE FALSE  TRUE FALSE FALSE FALSE
```

The logical vector returned denotes which elements in `A` are also in `B` and which are not.  

We saw previously that we could use the output from a logical expression to subset data by returning only the values corresponding to `TRUE`. Therefore, we can use the output logical vector to subset our data, and return only those elements in `A`, which are also in `B` by returning only the TRUE values:

![matching1](../img/in-operator1.png)

```r
intersection <- A %in% B
intersection
```

![matching2](../img/in-operator2.png)

```r
A[intersection]
```

![matching3](../img/in-operator3.png)

In these previous examples, the vectors were small and so it's easy to count by eye; but when we work with large datasets this is not practical. A quick way to assess whether or not we had any matches would be to use the `any` function to see if **any of the values contained in vector A are also in vector B**:

```r
any(A %in% B)
```

The `all` function is also useful. Given a logical vector, it will tell you whether all values returned are `TRUE`. If there is at least one `FALSE` value, the `all` function will return a `FALSE` and you know that all of A are not contained in B.

```r
all(A %in% B)
```
***
[**Exercise 1**](https://github.com/hbctraining/Intro-to-R/blob/master/results/answer_keys/07_matching_answer_key.md)

1. Using the `A` and `B` vectors created above, evaluate each element in `B` to see if there is a match in `A`

2. Subset the `B` vector to only return those values that are also in `A`.

***
Suppose we had **two vectors that had the same values but just not in the same order**. We could also use `all` to test for that. Rather than using the `%in%` operator we would use `==` and compare each element to the same position in the other vector. Unlike the `%in%` operator, **for this to work you must have two vectors that are of equal length**.

```r
A <- c(10,20,30,40,50)
B <- c(50,40,30,20,10)  # same numbers but backwards 

# test to see if each element of A is in B
A %in% B

# test to see if each element of A is in the same position in B
A == B

# use all() to check if they are a perfect match
all(A == B)

```

Let's try this on our data and see whether we have metadata information for all samples in our expression data. We'll start by creating two vectors; one with the `rownames` of the metadata and `colnames` of the RPKM data. These are base functions in R which allow you to extract the row and column names as a vector:

```r
x <- rownames(metadata)
y <- colnames(rpkm_data)
```

Now check to see that all of `x` are in `y`:

```r
all(x %in% y)
```

*Note that we can use nested functions in place of `x` and `y`:*

```r
all(rownames(metadata) %in% colnames(rpkm_data))
```

We know that all samples are present, but are they in the same order:

```r
all(rownames(metadata) == colnames(rpkm_data))
```

**Looks like all of the samples are there, but will need to be reordered. To reorder our genomic samples, we need to first learn different ways to reorder data. Therefore, we will step away from our genomic data briefly to learn about reordering, then return to it at the end of this lesson.**

***
[**Exercise 2**](https://github.com/hbctraining/Intro-to-R/blob/master/results/answer_keys/07_matching_answer_key.md)

We have a list of IDs for marker genes of particular interest. We want to extract count information associated with each of these genes, without having to scroll through our matrix of count data. We can do this using the `%in%` operator to extract the information for those genes from `rpkm_data`.

1. Create a vector for your important gene IDs, and use the `%in%`operator to determine whether these genes are contained in the row names of our `rpkm_data` dataset.

	```r
	important_genes <- c("ENSMUSG00000083700", "ENSMUSG00000080990", "ENSMUSG00000065619", "ENSMUSG00000047945", "ENSMUSG00000081010", 	"ENSMUSG00000030970")
	```
	
2. Extract the rows containing the important genes from your `rpkm_data` dataset using the `%in%` operator.

3. **Extra Credit:** Using the `important_genes` vector, extract the rows containing the important genes from your `rpkm_data` dataset without using the `%in%` operator.

***

## Reordering data using indexes
Indexing `[ ]` can be used to extract values from a dataset as we saw earlier, but we can also use it to rearrange our data values. 

```r
teaching_team <- c("Mary", "Meeta", "Radhika")
```
![reordering](../img/teachin-team.png)

Remember that we can return values in a vector by specifying it's position or index:

```r
teaching_team[c(2, 3)] # Extracting values from a vector
teaching_team
```

We can also extract the values and reorder them:

```r
teaching_team[c(3, 2)] # Extracting values and reordering them
```

Similarly, we can extract all of the values and reorder them:

```r
teaching_team[c(3, 1, 2)]
```

If we want to save our results, we need to assign to a variable:

```r
reorder_teach <- teaching_team[c(3, 1, 2)] # Saving the results to a variable
```

***
[**Exercise 3**](https://github.com/hbctraining/Intro-to-R/blob/master/results/answer_keys/07_matching_answer_key.md)

For a research project, we asked healthy volunteers and cancer patients questions about their diet and exercise. We also collected blood work for each individual, and each person was given a unique ID. Create the following dataframes, `behavior` and `blood` by copy/pasting the code below:

```r
# Creating behavior dataframe
	
ID <- c(546, 983, 042, 952, 853, 061)
diet <- c("veg", "pes", "omni", "omni", "omni", "omni")
exercise <- c("high", "low", "low", "low", "med", "high")
behavior <- data.frame(ID, diet, exercise)
	
# Creating blood dataframe
	
ID <- c(983, 952, 704, 555, 853, 061, 042, 237, 145, 581, 249, 467, 841, 546)
blood_levels <- c(43543, 465, 4634, 94568, 134, 347, 2345, 5439, 850, 6840, 5483, 66452, 54371, 1347)
blood <- data.frame(ID, blood_levels)
```

1. We would like to see if we have diet and exercise information for all of our blood samples. Using the `ID` information, determine whether all individuals with blood samples have associated behavioral information. Which individuals do not have behavioral information?

2. The samples lacking behavioral information correspond to individuals who opted out of the study after having their blood drawn. Subset the blood data to keep only samples that have behavioral information and save the dataframe back to the `blood` variable.

3. We would like to combine the `blood` and `behavior` dataframes together, but first we need to make sure the data is in the same order. 
	
	**a.** Take a look at each of the dataframes and manually identify the correct order for the blood dataframe such that it matches the order of IDs in the behavior dataframe. 
		
	**b.** Reorder the blood data to match the order of the IDs in the behavior dataframe and save the reordered blood dataframe as `blood_reordered`. *Hint: you will need to have a vector of index values from a. to reorder.* Once you have created `blood_reordered` you can use the `all()` function as a sanity check to make sure it was done correctly.
	
	**c.** Combine the dataframes blood_reordered and behavior using the data.frame() function and save this to a new dataframe called `blood_behavior`. *Note: you will find that there are now two "ID" columns, this will help verify that you have reordered correctly.*
	
***

## The `match` function

Now that we know how to reorder using indices, we can use the `match()` function to match the values in two vectors. We'll be using it to evaluate which samples are present in both our counts and metadata dataframes, and then to re-order the columns in the counts matrix to match the row names in the metadata matrix. 

`match()` takes at least 2 arguments: 

1. a vector of values in the order you want
2. a vector of values to be reordered

The function returns the position of the matches (indices) with respect to the second vector, which can be used to re-order it so that it matches the order in the first vector.  Let's create vectors `first` and `second` to demonstrate how it works:

```r
first <- c("A","B","C","D","E")
second <- c("B","D","E","A","C")  # same letters but different order
```
![matching4](../img/match1.png)

***How would you reorder `second` vector to match `first` using indexes?***

If we had large datasets, it would be difficult to reorder them by searching for the indices of the matching elements. This is where the `match` function comes in really handy:
	
```r
match(first,second)
[1] 4 1 5 2 3
```

The function should return a vector of size `length(first)`. Each number that is returned represents the index of the `second` vector where the matching value was observed. 

Let's change vector `second` so that only a subset are retained:

```r	
first <- c("A","B","C","D","E")
second <- c("D","B","A")  # remove values
```

![matching5](../img/match2.png)

And try to `match` again:

```r
match(first,second)

[1]  3  2 NA  1 NA
```

>**NOTE:** For values that don't match by default return an `NA` value. You can specify what values you would have it assigned using `nomatch` argument. Also, if there is more than one matching value found only the first is reported.

We can also reorder data using the output of the `match` function. We can reorder the `second` vector using the indexes from the `match` function of where the elements of the `first` vector occur in the `second` vector. First,  we save the match indexes to a variable:

```r
first <- c("A","B","C","D","E")
second <- c("B","D","E","A","C") 
```

![matching4](../img/match1.png)

```
reorder_idx <- match(first,second)
reorder_idx
[1] 4 1 5 2 3
```

![matching5](../img/match4-idx.png)

Now, we can just use the indexes to reorder the elements of the `second` vector to be in the same positions as the matching elements in the `first` vector:

```r
second[reorder_idx]  # Reordering the second vector to match the order of the first vector
second_reordered <- second[reorder_idx]  # Reordering and saving the output to a variable
```

![matching7](../img/match3-reordered.png)

***
[**Exercise 4**](https://github.com/hbctraining/Intro-to-R/blob/master/results/answer_keys/07_matching_answer_key.md)

Similar to the previous exercise, perform the reordering of the `blood` data to match the order of the IDs in the `behavior` dataframe, but this time use the `match()` function. Save the reordered blood dataframe as `blood_reordered_match`. 
***

### Reordering genomic data using `match` function
Using the `match` function, we now would like to *match the row names of our metadata to the column names of our expression data*, so these will be the arguments for `match`. Using these two arguments we will retrieve a vector of match indexes. The resulting vector represents the re-ordering of the column names in our data matrix to be identical to the rows in metadata:
 
 ```r
rownames(metadata)
	
colnames(rpkm_data)
	
genomic_idx <- match(rownames(metadata), colnames(rpkm_data))
genomic_idx
```

Now we can create a new data matrix in which columns are re-ordered based on the match indices:

```r
rpkm_ordered  <- rpkm_data[,genomic_idx]
```

Check and see what happened by using `head`. You can also verify that column names of this new data matrix matches the metadata row names by using the `all` function:

```r
head(rpkm_ordered)
all(rownames(metadata) == colnames(rpkm_ordered))
```

Now that our samples are ordered the same in our metadata and counts data, we could proceed to perform differential expression analysis with this dataset.


---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
