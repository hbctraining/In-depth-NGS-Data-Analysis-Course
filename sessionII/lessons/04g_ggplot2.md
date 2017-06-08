---
layout: topic
title: Plotting and data visualization in R
author: Mary Piper, Meeta Mistry, Radhika Khetani
date: "Tuesday, May 9, 2017"
---

Approximate time: 60 minutes

## Learning Objectives 

* Use the "apply" function for iterative tasks on data structures.
* Plot graphs using the external package "ggplot2".
* Export plots for use outside of the R environment.


## Calculating simple statistics

Let's take a closer look at our counts data by bringing it into a data frame in R. We are going to download a [new counts data file](https://raw.githubusercontent.com/hbctraining/In-depth-NGS-Data-Analysis-Course/may2017/sessionII/data/counts.rpkm_ordered.csv) into our data folder for plotting. 

```r
rpkm_data <- read.csv("data/counts.rpkm_ordered.csv", row.names = 1)

View(rpkm_data)
```

Each column represents a sample in our experiment, and each sample has ~38K values corresponding to the expression of different transcripts. Suppose we wanted to compute the average value of expression for a sample (across all transcripts), the R base package provides many built-in functions such as `mean`, `median`, `min`, `max`, and `range`, just to name a few. Try computing the mean for "sample1" (_Hint: apply what you have learned previously using indexes_)  

```r
mean(rpkm_data[,'sample1'])
```

> ### Missing values
> By default, all **R functions operating on vectors that contains missing data will return NA**. It's a way to make sure that users know they have missing data, and make a conscious decision on how to deal with it. When dealing with simple statistics like the mean, the easiest way to ignore `NA` (the missing data) is to use `na.rm=TRUE` (`rm` stands for remove). 
> In some cases, it might be useful to remove the missing data from the vector. For this purpose, R comes with the function `na.omit` to generate a vector that has NA's removed. For some applications, it's useful to keep all observations, for others, it might be best to remove all observations that contain missing data. The function
`complete.cases()` returns a logical vector indicating which rows have no missing values. 

## The `apply` Function
To obtain mean values for all samples we can use `mean` on each column individually, but there is also an easier way to go about it. Programming languages typically have a way to allow the execution of a single line of code or several lines of code multiple times, or in a "loop". By default R is not very good at looping, hence the `apply()` family of functions are used for this purpose. This family includes several functions, each differing slightly on the inputs or outputs. For example, we can use `apply()` to execute some task on every element in a vector, every row/column in a dataframe, and so on. 

```r
base::apply             Apply Functions Over Array Margins
base::by                Apply a Function to a Data Frame Split by Factors
base::eapply            Apply a Function Over Values in an Environment
base::lapply            Apply a Function over a List or Vector (returns list)
base::sapply            Apply a Function over a List or Vector (returns vector)
base::mapply            Apply a Function to Multiple List or Vector Arguments
base::rapply            Recursively Apply a Function to a List
base::tapply            Apply a Function Over a Ragged Array
```

We will be using `apply` in our examples today, but do take a moment on your own to explore the many options that are available. The `apply` function returns a vector or array or list of values obtained by applying a function to margins of an array or matrix. We know about vectors/arrays and functions, but what are these “margins”? Margins are referring to either the rows (denoted by 1), the columns (denoted by 2) or both (1:2). By “both”, we mean  apply the function to each individual value.

The syntax for the apply function is: 

```r
## DO NOT RUN
apply(dataframe/matrix, margin, function_to_apply)
```

Let's try this to obtain mean expression values for each sample in our RPKM matrix:

```r
samplemeans <- apply(rpkm_data, 2, mean) 
```

Now, add `samplemeans` to the end of the `metadata` dataframe:
	
```r
new_metadata <- cbind(metadata, samplemeans)
```

Before we start to plot, we also want to add an additional metadata column to `new_metadata`, this new column lists the age of each of the mouse samples in days.

```r
age_in_days <- c(40, 32, 38, 35, 41, 32, 34, 26, 28, 28, 30, 32)    
# Create a numeric vector with ages. Note that there are 12 elements here.
	
new_metadata <- cbind(new_metadata, age_in_days)    
# add the new vector as the last column to the new_metadata dataframe
```

We are now ready for plotting and data visualization!

## Data Visualization with `ggplot2`

When we are working with large sets of numbers it can be useful to display that information graphically to gain more insight. Visualization deserves an entire course of its own (there is that much to know!). If you are interested in learning about plotting with base R functions, we have a short lesson [available here](basic_plots_in_r.md). In this lesson we will be plotting with the popular Bioconductor package [`ggplot2`](http://docs.ggplot2.org/).

More recently, R users have moved away from base graphic options towards `ggplot2` since it offers a lot more functionality as compared to the base R plotting functions. The `ggplot2` syntax takes some getting used to, but once you get it, you will find it's extremely powerful and flexible. We will start with drawing a simple x-y scatterplot of `samplemeans` versus `age_in_days` from the `new_metadata` data frame. `ggplot2` assumes that the input is a  data frame.

Let's start by loading the `ggplot2` library.

```{r}
library(ggplot2)
```

The `ggplot()` function is used to **initialize the basic graph structure**, then we add to it. The basic idea is that you specify different parts of the plot, and add them together using the `+` operator. These parts are often referred to as layers.

Let's start: 

```{r, eval=FALSE}
ggplot(new_metadata) # what happens? 
```
You get an blank plot, because you need to **specify layers** using the `+` operator.

One type of layer is **geometric objects**. These are the actual marks we put on a plot. Examples include:

* points (`geom_point`, `geom_jitter` for scatter plots, dot plots, etc)
* lines (`geom_line`, for time series, trend lines, etc)
* boxplot (`geom_boxplot`, for, well, boxplots!)

For a more exhaustive list on all possible geometric objects and when to use them check out [Hadley Wickham's RPubs](http://rpubs.com/hadley/ggplot2-layers). 

A plot **must have at least one `geom`**; there is no upper limit. You can add a `geom` to a plot using the `+` operator

```{r, eval=FALSE}
ggplot(new_metadata) +
  geom_point() # note what happens here
```

You will find that even though we have added a layer by specifying `geom_point`, we get an error. This is because each type of geom usually has a **required set of aesthetics** to be set. Aesthetic mappings are set with the aes() function and can be set inside `geom_point()` to be specifically applied to that layer. If we supplied aesthetics within `ggplot()`, they will be used as defaults for every layer. Examples of aesthetics include:

* position (i.e., on the x and y axes)
* color ("outside" color)
* fill ("inside" color) 
* shape (of points)
* linetype
* size

To start, we will add position for the x- and y-axis since `geom_point` requires the most basic information about a scatterplot, i.e. what you want to plot on the x and y axes. All of the others mentioned above are optional.

```{r, fig.align='center'}
ggplot(new_metadata) +
     geom_point(aes(x = age_in_days, y= samplemeans))
```

 ![ggscatter1](../img/ggscatter-1.png) 


Now that we have the required aesthetics, let's add some extras like color to the plot. We can **`color` the points on the plot based on genotype**, by specifying the column header. You will notice that there are a default set of colors that will be used so we do not have to specify. Also, the **legend has been conveniently plotted for us!**


```{r, fig.align='center'}
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype)) 
```

 ![ggscatter1.1](../img/ggscatter-2.png) 


Alternatively, we could color based on celltype by changing it to `color =celltype`. Let's try something different and have both **celltype and genotype identified on the plot**. To do this we can assign the `shape` aesthetic the column header, so that each celltype is plotted with a different shaped data point. Add in `shape = celltype` to your aesthetic and see how it changes your plot:

```{r, fig.align='center'}
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
  			shape=celltype)) 
```

 ![ggscatter3](../img/ggscatter-3.png) 


The **size of the data points** are quite small. We can adjust that within the `geom_point()` layer, but does **not** need to be **included in `aes()`** since we are specifying how large we want the data points, rather than mapping it to a variable. Add in the `size` argument by specifying a number for the size of the data point:

```{r, fig.align='center'}
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
  			shape=celltype), size=3.0) 
```

 ![ggscatter4](../img/ggscatter-4.png)
  

The labels on the x- and y-axis are also quite small and hard to read. To change their size, we need to add an additional **theme layer**. The ggplot2 `theme` system handles non-data plot elements such as:

* Axis label aesthetics
* Plot background
* Facet label backround
* Legend appearance

There are built-in themes we can use (i.e. `theme_bw()`) that mostly change the background/foreground colours, by adding it as additional layer. Or we can adjust specific elements of the current default theme by adding the `theme()` layer and passing in arguments for the things we wish to change. Or we can use both.

Let's add a layer `theme_bw()`. Do the axis labels or the tick labels get any larger by changing themes?

```{r, fig.align='center'}
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
  			shape=celltype), size=3.0) +
  theme_bw() 
```

Not in this case. But we can add arguments using `theme()` to change it ourselves. Since we are adding this layer on top (i.e later in sequence), any features we change will override what is set in the `theme_bw()`. Here we'll **increase the size of the axes labels and axes tick labels to be 1.5 times the default size.** When modfying the size of text we often use the `rel()` function. In this way the size we specify is relative to the default (similar to `cex` for base plotting). We can also provide the number vaue as we did with the data point size, but can be cumbersome if you don't know what the default font size is to begin with. 

```{r, fig.align='center'}
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
  			shape=celltype), size=3.0) +
  theme_bw() +
  theme(axis.text = element_text(size=rel(1.5)),
  		axis.title = element_text(size=rel(1.5)))			
```
 
 ![ggscatter5](../img/ggscatter-5.png)
 

> *NOTE:* You can use the `example("geom_point")` function here to explore a multitude of different aesthetics and layers that can be added to your plot. As you scroll through the different plots, take note of how the code is modified. You can use this with any of the different geometric object layers available in ggplot2 to learn how you can easily modify your plots! 

> *NOTE:* RStudio provide this very [useful cheatsheet](https://www.rstudio.com/wp-content/uploads/2016/11/ggplot2-cheatsheet-2.1.pdf) for plotting using `ggplot2`. Different example plots are provided and the associated code (i.e which `geom` or `theme` to use in the appropriate situation.)


***

**Exercise**

1. The current axis label text defaults to what we gave as input to `geom_point` (i.e the column headers). We can change this by **adding additional layers** called `xlab()` and `ylab()` for the x- and y-axis, respectively. Add these layers to the current plot such that the x-axis is labeled "Age (days)" and the y-axis is labeled "Mean expression".
2. Use the `ggtitle` layer to add a title to your plot. *NOTE: Useful code to center your title over your plot can be done using `theme(plot.title=element_text(hjust=0.5))`.*

***



## Histogram

To plot a histogram we require another type of geometric object called `geom_histogram`, which requires a statistical transformation. Some plot types (such as scatterplots) do not require transformations, each point is plotted at x and y coordinates equal to the original value. Other plots, such as boxplots, histograms, prediction lines etc. need to be transformed. Usually these objects have has a default statistic for the transformation, but that can be changed via the `stat_bin` argument. 

Let's plot a histogram of sample mean expression in our data:

```r
ggplot(new_metadata) +
  geom_histogram(aes(x = samplemeans))
```


 ![gghist1](../img/gghist-1.png) 
 

You will notice that even though the histogram is plotted, R gives a warning message ``stat_bin()` using `bins = 30`. Pick better value with `binwidth`.` These are the transformations we discussed. Apparently the default is not good enough. 

Let's change the binwidth values. How does the plot differ?

```{r, fig.align='center'}
ggplot(new_metadata) +
  geom_histogram(aes(x = samplemeans), stat = "bin", binwidth=0.8)
```

 ![gghist2](../img/gghist-2.png) 

***

**Exercise**

1. Use what you learned previously about the `theme()` layer to make the text larger for x- and y-axis (not the tick labels).
2. Also, add an appropriate title to this plot.

***


## Boxplot

Now that we have all the required information for plotting with ggplot2 let's try plotting a boxplot. A boxplot provides a graphical view of the distribution of data based on a five number summary. The top and bottom of the box represent the (1) first and (2) third quartiles (25th and 75th percentiles, respectively). The line inside the box represents the (3) median (50th percentile). The whiskers extending above and below the box represent the (4) maximum, and (5) minimum of a data set. The whiskers of the plot reach the minimum and maximum values that are not outliers. 

Outliers are determined using the interquartile range (IQR), which is defined as: Q3 - Q1. Any values that exceeds 1.5 x IQR below Q1 or above Q3 are considered outliers and are represented as points above or below the whiskers. These outliers are useful to identify any unexpected observations.

1. Use the `geom_boxplot()` layer to plot the differences in sample means between the Wt and KO genotypes.
2. Add a title to your plot.
3. Add 'Genotype' as your x-axis label and 'Mean expression' as your y-axis labels.
4. Change the size of your axes labels to 1.5x larger than the default.
5. Change the size of your axes text (the labels on the tick marks) to 1.25x larger than the default.
6. Change the size of your plot title in the same way that you change the size of the axes text but use `plot.title`.

*BONUS: Use the `fill` aesthetic to look at differences in sample means between celltypes within each genotype.*


**Our final figure should look something like that provided below.**

 ![ggbox](../img/ggboxplot.png)


## Exporting figures to file

There are two ways in which figures and plots can be output to a file (rather than simply displaying on screen). The first (and easiest) is to export directly from the RStudio 'Plots' panel, by clicking on `Export` when the image is plotted. This will give you the option of `png` or `pdf` and selecting the directory to which you wish to save it to. It will also give you options to dictate the size and resolution of the output image.

The second option is to use R functions and have the write to file hard-coded in to your script. This would allow you to run the script from start to finish and automate the process (not requiring human point-and-click actions to save).  In R’s terminology, **output is directed to a particular output device and that dictates the output format that will be produced**.  A device must be created or “opened” in order to receive graphical output and, for devices that create a file
on disk, the device must also be closed in order to complete the output.

Let's print our scatterplot to a pdf file format. First you need to initialize a plot using a function which specifies the graphical format you intend on creating i.e.`pdf()`, `png()`, `tiff()` etc. Within the function you will need to specify a name for your image, and the with and height (optional). This will open up the device that you wish to write to:
```r
pdf("figures/scatterplot.pdf")
```

If you wish to modify the size and resolution of the image you will need to add in the appropriate parameters as arguments to the function when you initialize. Then we plot the image to the device, using the ggplot scatterplot that we just created. 

```r
ggplot(new_metadata) +
  geom_point(aes(x = age_in_days, y= samplemeans, color = genotype,
  			shape=celltype), size=rel(3.0)) 
```

Finally, close the "device", or file, using the `dev.off()` function. There are also `bmp`, `tiff`, and `jpeg` functions, though the jpeg function has proven less stable than the others. 
  			
```r    
dev.off()
```

***Note 1:*** *You will not be able to open and look at your file using standard methods (Adobe Acrobat or Preview etc.) until you execute the `dev.off()` function.*

***Note 2:*** *If you had made any additional plots before closing the device, they will all be stored in the same file; each plot usually gets its own page, unless you specify otherwise.*



---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
