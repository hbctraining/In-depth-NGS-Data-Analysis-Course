---
title: "Nested Functions"
authors: Meeta Mistry, Mary Piper
date: "Wednesday, September 28, 2016"
---
Approximate time: 30 min

## Learning Objectives
* Understand and implement functions nested within other functions


## Nested functions

Thus far, to perform any specific task, we have executed every function separately; if we wanted to use the results of a function for downstream purposes, we saved the results to a variable. As you become more comfortable with R, you will find that it is more efficient to code using nested functions, or functions within other functions, which will allow you to execute multiple commands at the same time.

Let's start with an example from the previous lesson to demonstrate. To obtain the rows in `metadata` which correspond to celltype A, we used two lines of code: 

```r
idx <- metadata$celltype == "typeA"

metadata[idx, ]
```

However, we could have also done this in a single line of code and avoid having to create the variable `idx`:

```r
metadata[metadata$celltype == "typeA",]
```

This is a rather simple example combining only two lines of code, but you see how the code becomes lengthy and slightly more difficult to read. Even if you decide to avoid writing nested functions for the time being, you should still have experience reading and understanding them. The key to understanding nested functions is to **read from the inside out**.

Let's work through some examples of nested functions!

### Nested functions practice #1

Based on our metadata, how many samples are Wt *and* from celltype A?

**Step 1:** Classifying the samples using criteria based on logical operators:

```r
vec <- metadata$celltype == "typeA" & metadata$genotype == "Wt"
```

**Step 2:** Identifying which samples meet our criteria (i.e the `TRUE` samples):

```r
true_vec <- which(vec)
```

**Step 3:** Counting how many samples meet our criteria:

```r
length(true_vec)
```

**Nested code:**
Rather that assigning the output from each step to a separate variable, we could also just nest them into one another. **To create a nested function, simply replace the variable name with the code on the right hand side of the assignment operator.** *Don't forget to keep track of the sets of parentheses!*

```r
length(which(metadata$celltype == "typeA" & metadata$genotype == "Wt"))
```

### Nested functions practice #2
You realize that you forgot to include important metadata regarding sex of your samples in your `metadata` file. You would like to add this data in and create a new `metadata_new` dataframe. Using functions separately would require us to execute three separate steps:  

**Step 1:** Create the `sex` vector: 
	
```r
sex <- c("M","F","M","M","F","M","M","F","M","M","F","M")
```

**Step 2:** Turn the `sex` vector into a factor variable:
	 
```r
sex_fr <- factor(sex)
```

**Step 3:** Use the `cbind` function to add the column to the **end** of the `metadata` dataframe: 

```r
metadata_new <- cbind(metadata, sex=sex_fr)
```

Instead of performing all three steps, we would like to create a nested function. We could first replace the `sex_fr` variable with the assignment in **Step2**:

```r
metadata_new <- cbind(metadata, sex=factor(sex))
```

Or we can go a step further and combine all steps, making your code slightly more difficult to read (we don't recommend doing this):

```r
metadata_new <- cbind(metadata,
			sex=factor(c("M","F","M","M","F","M","M","F","M","M","F","M")))
```

### Nested functions practice #3			
Now, let's say that you are interested in finding out which samples (listed by sample name) in your dataset have  "Wt" genotype within our `metadata` file. For our small dataset, we can simply do this by eye but for larger datasets it is easier to write code to do so:

**Step 1:** Determine the row **locations** in `metadata` for those samples with `genotype` equal to "Wt":
	
```r
rloc <- which(metadata$genotype == "Wt")
```

**Step 2:** Obtain the vector of row names from `metadata`:
	
```r
rnames <- row.names(metadata)
```

**Step 3:** Identify the sample names by using the indexes determined in **Step 2**:

```r
wt_samples <- rnames[rloc]
```

Alternatively, we could combine the steps:

```r
wt_samples <- row.names(metadata)[which(metadata$genotype == "Wt")]
```

Learning to understand nested functions is a critical part of your mastery of R. Not only will their use improve your efficiency, but nested functions are frequently encountered in help forums and R package documentation, so understanding them is critical to your learning process. 

***
**Exercise**

Use nested functions and the `subset()` function to find out which samples (listed by sample name) in your dataset have  "typeB" cell type within our `metadata` file.

---

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson is adapted from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
