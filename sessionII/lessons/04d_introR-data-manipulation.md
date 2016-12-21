---
title: "Data manipulation"
authors: Meeta Mistry, Mary Piper, Radhika Khetani
date: "Wednesday, September 28, 2016"
---
Approximate time: 60 min

## Learning Objectives
* Reading data into R
* Inspecting data structures
* Using indices and sequences to select data from vectors


## Reading data into R
Regardless of the specific analysis in R we are performing, we usually need to bring data in for the analysis. The function in R we use will depend on the type of data file we are bringing in (e.g. text, Stata, SPSS, SAS, Excel, etc.) and how the data in that file are separated, or delimited. The table below lists functions that can be used to import data from common file formats.

| Data Type  | Function | Package
| -----------:|:----------------:|:---------------:|
| comma separated (.csv)  | read.csv()	| utils (default) |
| other delimited formats (.txt) | read.table(); read.csv() | utils (default) |
| Stata version 7-12 (.dta) | read.dta() | foreign |
| Stata version 13-14 (.dta) | readdta() | haven |
| SPSS (.sav) |	read.spss() | foreign |
| SAS (.sas7bdat) | read.sas7bdat() | sas7bdat |
| Excel (.xls, .xlsx) | readWorksheetFromFile() | XLConnect |
 

For example, if we have text file separated by commas (comma-separated values), we could use the function `read.csv`. However, if the data are separated by a different delimiter in a text file, we could use the generic `read.table` function and specify the delimiter as an argument in the function. 

When working with genomic data, we often have a metadata file containing information on each sample in our dataset. Let's bring in the metadata file using the `read.csv` function. Check the arguments for the function to get an idea of the function options:

```r
?read.csv
```

The `read.csv` function has *one required argument* and several *options* that can be specified. The mandatory argument is a path to the file and filename, which in our case is `data/mouse_exp_design.csv`. We will put the function to the right of the assignment operator, meaning that **any output will be saved as the variable name provided on the left**.

```r
metadata <- read.csv(file="data/mouse_exp_design.csv")
```

> *Note: By default, `read.csv` converts (= coerces) columns that contain characters (i.e., text) into the `factor` data type. Depending on what you want to do with the data, you may want to keep these columns as `character`. To do so, `read.csv()` and `read.table()` have an argument called `stringsAsFactors` which can be set to `FALSE`.*
> 


## Inspecting data structures

There are a wide selection of base functions in R that are useful for inspecting your data and summarizing it. Let's use the `metadata` file that we created to test out data inspection functions. 

Take a look at the dataframe by typing out the variable name `metadata` and pressing return; the variable contains information describing the samples in our study. Each row holds information for a single sample, and the columns contain categorical information about the sample `genotype`(WT or KO),  `celltype` (typeA or typeB), and `replicate number` (1,2, or 3).


```r
metadata

          genotype celltype replicate
sample1        Wt    typeA		1
sample2        Wt    typeA		2
sample3        Wt    typeA		3
sample4        KO    typeA		1
sample5        KO    typeA		2
sample6        KO    typeA		3
sample7        Wt    typeB		1
sample8        Wt    typeB		2
sample9        Wt    typeB		3
sample10       KO    typeB		1
sample11       KO    typeB		2
sample12       KO    typeB		3

```

Suppose we had a larger file, we might not want to display all the contents in the console. Instead we could check the top (the first 6 lines) of this `data.frame` using the function `head()`:

```r
head(metadata)
```

Previously, we had mentioned that character values get converted to factors by default using `data.frame`. One way to assess this change would be to use the __`str`__ucture function. You will get specific details on each column:


```r
str(metadata)

'data.frame':	12 obs. of  3 variables:
 $ genotype : Factor w/ 2 levels "KO","Wt": 2 2 2 1 1 1 2 2 2 1 ...
 $ celltype : Factor w/ 2 levels "typeA","typeB": 1 1 1 1 1 1 2 2 2 2 ...
 $ replicate: num  1 2 3 1 2 3 1 2 3 1 ...
```

As you can see, the columns `genotype` and `celltype` are of the `factor` class, whereas the replicate column has been interpreted as integer data type.

__You can also get this information from the "Environment" tab in RStudio.__

### List of functions for data inspection

We already saw how the functions `head()` and `str()` can be useful to check the
content and the structure of a `data.frame`. Here is a non-exhaustive list of
functions to get a sense of the content/structure of data.

* All data structures - content display:
	- **`str()`:** compact display of data contents (env.)
	- **`class()`:** data type (e.g. character, numeric, etc.) of vectors and data structure of dataframes, matrices, and lists.
	- **`summary()`:** detailed display, including descriptive statistics, frequencies
	- **`head()`:** will print the beginning entries for the variable
	- **`tail()`:** will print the end entries for the variable
* Vector and factor variables: 
	- **`length()`:** returns the number of elements in the vector or factor
* Dataframe and matrix variables:
	- **`dim()`:** returns dimensions of the dataset
	- **`nrow()`:** returns the number of rows in the dataset
	- **`ncol()`:** returns the number of columns in the dataset
	- **`rownames()`:** returns the row names in the dataset  
	- **`colnames()`:** returns the column names in the dataset

## Selecting data using indices and sequences

When analyzing data, we often want to **partition the data so that we are only working with selected columns or rows.** A data frame or data matrix is simply a collection of vectors combined together. So let's begin with vectors and how to access different elements, and then extend those concepts to dataframes.

### Vectors

#### Selecting using indices

If we want to extract one or several values from a vector, we must provide one or several indices using square brackets `[ ]` syntax. The **index represents the element number within a vector** (or the compartment number, if you think of the bucket analogy). R indices start at 1. Programming languages like Fortran, MATLAB, and R start counting at 1, because that's what human beings typically do. Languages in the C family (including C++, Java, Perl, and Python) count from 0 because that's simpler for computers to do.

Let's start by creating a vector called age:

```r
age <- c(15, 22, 45, 52, 73, 81)
```

![vector indices](../img/vector-index.png)

Suppose we only wanted the fifth value of this vector, we would use the following syntax:

```r
age[5]
```

If we wanted all values except the fifth value of this vector, we would use the following:

```r
age[-5]
```

If we wanted to select more than one element we would still use the square bracket syntax, but rather than using a single value we would pass in a *vector of several index values*:

```r
idx <- c(3,5,6) # create vector of the elements of interest
age[idx]
```

To select a sequence of continuous values from a vector, we would use `:` which is a special function that creates numeric vectors of integer in increasing or decreasing order. Let's select the *first four values* from age:

```r
age[1:4]
```

Alternatively, if you wanted the reverse could try `4:1` for instance, and see what is returned. 

***
**Exercises** 

1. Create a vector called alphabets with the following alphabets, C, D, X, L, F.
2. Use the associated indices along with `[ ]` to do the following:
	* only display C, D and F
	* display all except X
	* (optional) display the alphabets in the opposite order (F, L, X, D, C)

***

#### Selecting using indices with logical operators

We can also use indices with logical operators. Logical operators include greater than (>), less than (<), and equal to (==). A full list of logical operators in R is displayed below:

| Operator | Description |
| :-----------:|:----------------|
| > | greater than |
| >= | greater than or equal to|
| < | less than |
| <= | less than or equal to |
| == | equal to |
| != | not equal to |
| & | and |
| \| |or |

We can use logical expressions to determine whether a particular condition is true or false. For example, let's use our age vector: 
	
```r
age
```

If we wanted to know if each element in our age vector is greater than 50, we could write the following expression:	

```r
age > 50
```

Returned is a vector of logical values the same length as age with TRUE and FALSE values indicating whether each element in the vector is greater than 50.

```r
[1] FALSE FALSE FALSE  TRUE  TRUE  TRUE
```

We can use these logical vectors to select only the elements in a vector with TRUE values at the same position or index as in the logical vector.

Create an index with logical operators to select all values in the `age` vector over 50 **or** `age` less than 18:

```r
idx <- age > 50 | age < 18
	
idx
	
age

age[idx]
```

##### Indexing with logical operators using the `which()` function

While logical expressions will return a vector of TRUE and FALSE  values of the same length, we could use the `which()` function to output the indices where the values are TRUE. Indexing with either method generates the same results, and personal preference determines which method you choose to use. For example:

```r
idx <- which(age > 50 | age < 18)

idx

age[idx]
```

Notice that we get the same results regardless of whether or not we use the `which()`. Also note that while `which()` works the same as the logical expressions for indexing, it can be used for multiple other operations, where it is not interchangeable with logical expressions.


### Factors

Since factors are special vectors, the same rules for selecting values using indices apply. The elements of the expression factor created previously had the following categories or levels: low, medium, and high. 

Let's extract the values of the factor with high expression:

First, we create a logical vector of TRUE and FALSE values:

```r
idx <- expression == "high"
```

Then, we use the brackets [ ] to extract the TRUE values from the dataset:

```r
expression[idx]
```
	
#### Releveling factors

We have briefly talked about factors, but this data type only becomes more intuitive once you've had a chance to work with it.  Let's take a slight detour and learn about how to **order and relevel categories within a factor**. As we learned earlier, the categories in the `expression` factor were assigned integers alphabetically, with high=1, low=2, medium=3. To view the integer assignments under the hood you can use str:

```r
str(expression)
Factor w/ 3 levels "high","low","medium": 2 1 3 1 2 3 1
```

The unique elements are referred to as "factor levels".

In the example above, the factor has levels but it is unordered, i.e. there is no notation to say that high is greater than medium etc. In fact, the high category is the middle category because of alphabetical order. 

To order factor levels, you can add an argument to the `factor()` function, ordered=TRUE:

```r
expression <- factor(expression, ordered=TRUE)    ## Note that the `factor()` function is used to create a factor, & to modify the characteristics of an existing factor

str(expression)
Ord.factor w/ 3 levels "low"<"high"<..: 1 3 2 3 1 2 3
```

Now the output of the `str()` function states that this is an "Ord.factor", and there are "<" signs to denote that low is the lowest category. 

However, the order of categories is still incorrect, because R is ordering them alphabetically. R does not consider the meaning of the words here, so we have to coerce the proper ordering (i.e. "low" < "medium" < "high") using the `factor()` function once again.

```r
expression <- factor(expression, levels=c("low", "medium", "high"), ordered=TRUE)    
	
str(expression)
```

***
**Exercise**

1. Extract only the elements in `samplegroup` that are not KO.
2. Use the `samplegroup` vector we created in a previous lesson, and change that to an ordered factor such that KO < CTL < OE. 

***

---

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson is adapted from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
