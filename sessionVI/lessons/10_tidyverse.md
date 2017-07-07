The tidyverse
================
Michael J. Steinbaugh
2017-07-07

The [tidyverse](http://tidyverse.org/) is a suite of integrated packages designed to make common operations performed in [R](https://www.r-project.org/) more user friendly. This was initially developed by [Hadley Wickham](http://hadley.nz), Chief Scientist at [RStudio](https://www.rstudio.com/), but is now maintained by a number of talented developers moving the [R](https://www.r-project.org/) language forward.

![](../img/tidyverse_website.png)

------------------------------------------------------------------------

Installation
============

The core collection of tidyverse packages are managed by the [tidyverse CRAN package](https://cran.r-project.org/web/packages/tidyverse/index.html).

``` r
install.packages("tidyverse")
```

``` r
library(tidyverse)
```

    ## Loading tidyverse: ggplot2
    ## Loading tidyverse: tibble
    ## Loading tidyverse: tidyr
    ## Loading tidyverse: readr
    ## Loading tidyverse: purrr
    ## Loading tidyverse: dplyr

    ## Conflicts with tidy packages ----------------------------------------------

    ## filter(): dplyr, stats
    ## lag():    dplyr, stats

When you load the `tidyverse` library, these core packages will be loaded into your environment:

-   [ggplot2](http://ggplot2.tidyverse.org/)
-   [tibble](http://tibble.tidyverse.org/)
-   [tidyr](http://tidyr.tidyverse.org/)
-   [readr](http://readr.tidyverse.org/)
-   [purr](http://purrr.tidyverse.org/)
-   [dplyr](http://dplyr.tidyverse.org/)

Recommended optional packages
-----------------------------

There are a number of additional [tidyverse](http://tidyverse.org/) packages that we highly recommend for performing data analysis, including:

-   [magrittr](http://magrittr.tidyverse.org/): Defines the pipe operator (`%>%`), which is used to write left-to-right chain operations. We'll cover this below.
-   [stringr](http://stringr.tidyverse.org/): Enables easier manipulation of vectors ("strings").
-   [readxl](http://readxl.tidyverse.org/): Current recommended best practice for import of Excel workbooks.

Function name conflicts
-----------------------

When you load the [tidyverse](http://tidyverse.org/), you'll see messages about tidy package conflicts at the end. *This is normal.* Currently, [dplyr](http://dplyr.tidyverse.org/) masks `stats::filter()` and `stats::lag()` by default. What's happening here is that the [tidyverse](http://tidyverse.org/) has some functions with the same names as base [R](https://www.r-project.org/) packages.

**Note**: This remains a common issue when loading multiple libraries in a single script. For example, many [Bioconductor](https://bioconductor.org/) packages have generic functions with the same name as base [R](https://www.r-project.org/) and [tidyverse](http://tidyverse.org/) packages. For example, `biomaRt::select()` and `dplyr::select()` have the same function name but require different arguments. Whichever library you load last will define the function name (`select()`). If you need to use two packages with the same function name at the same time, you can reference them explicitly (e.g. `dplyr::select()`). Therefore, when starting a new analysis using [tidyverse](http://tidyverse.org/) packages, we highly recommend slotting `library(tidyverse)` at the end of your list of libraries.

tibbles
=======

A core component of the [tidyverse](http://tidyverse.org/) is the [tibble](http://tibble.tidyverse.org/). Tibbles are a modern rework of the standard `data.frame`, with some internal improvements to make code more reliable. Most notably, tibbles will return a more reasonable number of rows to the console.

Before we begin the tidy workflow, coerce your counts, metadata, and DESeqResults to data frames:

``` r
counts <- as.data.frame(data)
meta <- as.data.frame(meta)
results <- as.data.frame(res_tableOE)
```

Now let's create tibbles from our data frames:

``` r
counts_tbl <- counts %>%
    as_tibble %>%
    rownames_to_column("ensgene")
meta_tbl <- meta %>%
    as_tibble %>%
    rownames_to_column("sample_name") %>%
    rename(sample_type = sampletype,
           mov_expression = MOVexpr) %>%
    mutate(sample_name = tolower(sample_name))
results_tbl <- results %>%
    as_tibble %>%
    rownames_to_column("symbol")
```

First, try returning the `counts` data frame in your console.

``` r
counts
```

Next, try returning the counts [tibble](http://tibble.tidyverse.org/).

``` r
counts_tbl
```

    ## # A tibble: 38,828 x 13
    ##    ensgene    sample1   sample2   sample3   sample4    sample5    sample6
    ##      <chr>      <dbl>     <dbl>     <dbl>     <dbl>      <dbl>      <dbl>
    ##  1       1 19.7848000 19.265000 20.889500 24.076700 23.7222000 20.8198000
    ##  2       2  0.0000000  0.000000  0.000000  0.000000  0.0000000  0.0000000
    ##  3       3  0.9377920  1.032290  0.892183  0.827891  0.8269540  1.1686300
    ##  4       4  0.0359631  0.000000  0.000000  0.000000  0.0000000  0.0511932
    ##  5       5  0.1514170  0.056033  0.146196  0.180883  0.0473238  0.1438840
    ##  6       6  0.2567330  0.258134  0.421286  2.191720  1.0730200  1.6853800
    ##  7       7  5.9998100  6.047990  6.020250  5.620120  6.4116300  5.5177700
    ##  8       8  5.2784700  3.971810  6.161450  7.045910  5.2136600  7.3951100
    ##  9       9 37.2141000 32.303000 31.249600  4.260570  3.5634100  5.3828200
    ## 10      10 25.4044000 22.950700 21.415600 16.825400 17.9712000 17.7215000
    ## # ... with 38,818 more rows, and 6 more variables: sample7 <dbl>,
    ## #   sample8 <dbl>, sample9 <dbl>, sample10 <dbl>, sample11 <dbl>,
    ## #   sample12 <dbl>

See how [R](https://www.r-project.org/) only prints 10 rows instead of returning all 38k? This is much more user friendly.

Internally, a [tibble](http://tibble.tidyverse.org/) is essentially a class variant of `data.frame`, with some extra tibble (`tbl`) magic baked in:

``` r
class(meta_tbl)
```

    ## [1] "tbl_df"     "tbl"        "data.frame"

`glimpse()` is a modern rework of `str()`, optimized for tibbles:

``` r
glimpse(counts_tbl)
glimpse(meta_tbl)
glimpse(results_tbl)
```

Row names
---------

*Important*: [tidyverse](http://tidyverse.org/) is very opininationed about row names. These packages insist that all column data (e.g. `data.frame`) be treated equally, and that special designation of a column as `rownames` should be deprecated. [tibble](http://tibble.tidyverse.org/) provides simple utility functions to to handle rownames: `rownames_to_column()` and `column_to_rownames()`.

``` r
help("rownames", "tibble")
```

Code style
==========

One problem with [R](https://www.r-project.org/) is a lack of consistency across packages in how functions and arguments are named.

-   Base [R](https://www.r-project.org/) functions are formatted in dotted case: `read.csv()`.
-   [tidyverse](http://tidyverse.org/) functions are formatted in snake\_case: `read_csv()`.
-   [Bioconductor](https://bioconductor.org/) functions are generally formatted in lowerCamelCase (and sometimes UpperCamelCase).

The [tidyverse](http://tidyverse.org/) collection of packages are very opinionated in this regard and consistently use `snake_case` formatting for all function names and arguments. When using these functions, we recommend that you follow the [tidy style guide](http://style.tidyverse.org/).

Non-standard evaluation
=======================

[tidyverse](http://tidyverse.org/) packages improve code readability by changing how functions interpret object names. This is achieved through the use of "non-standard evaluation" instead of base [R](https://www.r-project.org/)'s "standard evaluation". This probably sounds confusing but is actually pretty simple. In fact, we've already used a function (`subset()`) in the class that works with non-standard evaluation.

``` r
subset(meta_tbl, mov_expression == "high")
```

    ## # A tibble: 3 x 3
    ##   sample_name          sample_type mov_expression
    ##         <chr>                <chr>          <chr>
    ## 1  mov10_oe_1 MOV10_overexpression           high
    ## 2  mov10_oe_2 MOV10_overexpression           high
    ## 3  mov10_oe_3 MOV10_overexpression           high

Here's the [tidyverse](http://tidyverse.org/) variant:

``` r
filter(meta_tbl, mov_expression == "high")
```

    ## # A tibble: 3 x 3
    ##   sample_name          sample_type mov_expression
    ##         <chr>                <chr>          <chr>
    ## 1  mov10_oe_1 MOV10_overexpression           high
    ## 2  mov10_oe_2 MOV10_overexpression           high
    ## 3  mov10_oe_3 MOV10_overexpression           high

See how both functions refer to the `mov_expression` column directly and not in quotations? That's *non-standard evaluation*. It makes code easier to read.

dplyr
=====

The most useful tool in the [tidyverse](http://tidyverse.org/) is [dplyr](http://dplyr.tidyverse.org/). It's a swiss-army knife for data manipulation. [dplyr](http://dplyr.tidyverse.org/) has 5 core functions that we recommend incorporating into your analysis:

-   `filter()` picks cases based on their values.
-   `arrange()` changes the ordering of the rows.
-   `select()` picks variables based on their names.
-   `mutate()` adds new variables that are functions of existing variables
-   `summarise()` reduces multiple values down to a single summary.

**Note:** [dplyr](http://dplyr.tidyverse.org/) underwent a massive revision this year, switching versions from 0.5 to 0.7. If you consult other [dplyr](http://dplyr.tidyverse.org/) tutorials online, note that many materials developed prior to 2017 are no longer correct. In particular, this applies to writing functions with [dplyr](http://dplyr.tidyverse.org/) (see Notes section below).

Let's make a report tibble of our `DESeqResults`.

`select()`
----------

First, we only need a few columns of interest from our results tibble. We can do this easily in [dplyr](http://dplyr.tidyverse.org/) with `select()`.

``` r
report <- results_tbl %>%
    select(symbol, baseMean, log2FoldChange, padj)
```

Conversely, you can remove columns you don't want with negative selection.

``` r
select(results_tbl, -c(lfcSE, stat, pvalue))
```

    ## # A tibble: 23,368 x 4
    ##         symbol    baseMean log2FoldChange         padj
    ##          <chr>       <dbl>          <dbl>        <dbl>
    ##  1 1/2-SBSRNA4  45.6520399    0.266586547 2.708964e-01
    ##  2        A1BG  61.0931017    0.208057615 3.638671e-01
    ##  3    A1BG-AS1 175.6658069   -0.051825739 7.837586e-01
    ##  4        A1CF   0.2376919    0.012557390           NA
    ##  5       A2LD1  89.6179845    0.343006364 7.652553e-02
    ##  6         A2M   5.8600841   -0.270449534 2.318666e-01
    ##  7       A2ML1   2.4240553    0.236041349           NA
    ##  8       A2MP1   1.3203237    0.079525469           NA
    ##  9      A4GALT  64.5409534    0.795049160 2.875565e-05
    ## 10       A4GNT   0.1912781    0.009458374           NA
    ## # ... with 23,358 more rows

`arrange()`
-----------

Note that the rows are sorted by the gene symbol. Let's fix that and sort them by adjusted P value instead with `arrange()`.

``` r
report <- arrange(report, padj)
```

`filter()`
----------

Let's keep only genes that are expressed (`baseMean` above 0) with an adjusted P value below 0.01. You can perform multiple `filter()` operations together in a single command.

``` r
report <- report %>%
    filter(baseMean > 0,
           padj < 0.01)
```

`mutate()`
----------

`mutate()` enables you to create a new column from an existing column. Let's generate log10 calculations of our baseMeans for each gene.

``` r
report %>%
    mutate(log10BaseMean = log10(baseMean)) %>%
    select(symbol, baseMean, log10BaseMean)
```

    ## # A tibble: 4,909 x 3
    ##      symbol   baseMean log10BaseMean
    ##       <chr>      <dbl>         <dbl>
    ##  1    MOV10 21681.7998      4.336095
    ##  2     H1F0  7881.0811      3.896586
    ##  3    HSPA6   168.2522      2.225961
    ##  4 HIST1H1C  1741.3830      3.240894
    ##  5    TXNIP  5133.7486      3.710435
    ##  6    NEAT1 21973.7061      4.341903
    ##  7    KLF10  1694.2109      3.228967
    ##  8   INSIG1 11872.5106      4.074543
    ##  9    NR1D1   969.9119      2.986732
    ## 10    WDFY1  1422.7361      3.153124
    ## # ... with 4,899 more rows

`rename()`
----------

You can quickly rename an existing column with `rename()`. The syntax is `new_name` = `old_name`.

``` r
report %>%
    rename(gene = symbol)
```

    ## # A tibble: 4,909 x 4
    ##        gene   baseMean log2FoldChange          padj
    ##       <chr>      <dbl>          <dbl>         <dbl>
    ##  1    MOV10 21681.7998      4.7695983  0.000000e+00
    ##  2     H1F0  7881.0811      1.5250811 2.007733e-162
    ##  3    HSPA6   168.2522      4.4993734 1.969313e-134
    ##  4 HIST1H1C  1741.3830      1.4868361 5.116720e-101
    ##  5    TXNIP  5133.7486      1.3868320  4.882246e-90
    ##  6    NEAT1 21973.7061      0.9087853  2.269464e-83
    ##  7    KLF10  1694.2109      1.2093969  9.257431e-78
    ##  8   INSIG1 11872.5106      1.2260848  8.853278e-70
    ##  9    NR1D1   969.9119      1.5236259  1.376753e-64
    ## 10    WDFY1  1422.7361      1.0629160  1.298076e-61
    ## # ... with 4,899 more rows

`summarise()`
-------------

You can perform column summarization operations with `summarise()`.

``` r
report %>%
    summarise(avgBaseMean = mean(baseMean))
```

    ## # A tibble: 1 x 1
    ##   avgBaseMean
    ##         <dbl>
    ## 1      1911.6

*Advanced:* `summarise()` is particularly powerful in combination with the `group_by()` function, which allows you to group related rows together.

*Note*: `summarize()` also works if you prefer to use American English. This applies across the board to any tidy functions, including in [ggplot2](http://ggplot2.tidyverse.org/) (e.g. `color` in place of `colour`).

`pull()`
--------

In the recent [dplyr](http://dplyr.tidyverse.org/) 0.7 update, `pull()` was added as a quick way to access column data as a vector. This is very handy in chain operations with the pipe operator.

``` r
pull(report, symbol) %>% .[1:10]
```

Joins
-----

To demonstrate [dplyr](http://dplyr.tidyverse.org/)'s powerful suite of join operations, let's import Ensembl gene annotations from the \[annotables\]\[\] package and add them to our report.

``` r
install.packages("devtools")
devtools::install_github("stephen_turner/annotables")
```

``` r
library(annotables)
annotable <- grch37 %>%
    select(symbol, biotype, description) %>%
    distinct
```

``` r
report <- left_join(report, annotable, by = "symbol")
```

------------------------------------------------------------------------

Notes
=====

Programming
-----------

Underneath the hood, [tidyverse](http://tidyverse.org/) packages build upon the base [R](https://www.r-project.org/) language using [rlang](https://github.com/tidyverse/rlang/), which is a **complete rework** of how functions handle variable names and evaluate arguments. This is achieved through the `tidyeval` framework, which interprates command operations using `tidy evaluation`. This is outside of the scope of the course, but explained in detail in the [Programming with dplyr](http://dplyr.tidyverse.org/articles/programming.html) vignette, in case you'd like to understand how these new tools behave differently from base [R](https://www.r-project.org/).

Additional resources
====================

-   [R for Data Science](http://r4ds.had.co.nz)
-   [teach the tidyverse](http://varianceexplained.org/r/teach-tidyverse/)
-   [tidy style guide](http://style.tidyverse.org/)
