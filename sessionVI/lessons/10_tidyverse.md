The tidyverse
================
Michael J. Steinbaugh
2017-07-07

The [tidyverse](http://tidyverse.org/) is a suite of integrated packages designed to make common operations performed in [R](https://www.r-project.org/) more user friendly.

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

Load the example data
=====================

``` r
results_dir <- file.path(
    "https://raw.githubusercontent.com",
    "hbctraining",
    "In-depth-NGS-Data-Analysis-Course",
    "may2017",
    "sessionVI",
    "results")
counts_mat <- file.path(results_dir, "counts.txt") %>%
    read.table %>%
    as.matrix
counts_tbl <- counts_mat %>%
    as.data.frame %>%
    as_tibble %>%
    rownames_to_column("ensgene")
meta_tbl <- file.path(results_dir, "meta.txt") %>%
    read.table %>%
    as_tibble %>%
    rownames_to_column("sample_name") %>%
    rename(sample_type = sampletype,
           mov_expression = MOVexpr) %>%
    mutate(sample_name = tolower(sample_name))
results_tbl <- file.path(results_dir, "results_Mov10_oe.txt") %>%
    read.table %>%
    as_tibble %>%
    rownames_to_column("symbol")
```

tibbles
=======

A core component of the [tidyverse](http://tidyverse.org/) is the [tibble](http://tibble.tidyverse.org/). Tibbles are a modern rework of the standard `data.frame`, with some internal improvements to make code more reliable. Most notably, tibbles will return a more reasonable number of rows to the console.

Try executing `counts` in your console:

``` r
counts_tbl
```

    ## # A tibble: 38,828 x 13
    ##               ensgene    sample1   sample2   sample3   sample4    sample5
    ##                 <chr>      <dbl>     <dbl>     <dbl>     <dbl>      <dbl>
    ##  1 ENSMUSG00000000001 19.7848000 19.265000 20.889500 24.076700 23.7222000
    ##  2 ENSMUSG00000000003  0.0000000  0.000000  0.000000  0.000000  0.0000000
    ##  3 ENSMUSG00000000028  0.9377920  1.032290  0.892183  0.827891  0.8269540
    ##  4 ENSMUSG00000000031  0.0359631  0.000000  0.000000  0.000000  0.0000000
    ##  5 ENSMUSG00000000037  0.1514170  0.056033  0.146196  0.180883  0.0473238
    ##  6 ENSMUSG00000000049  0.2567330  0.258134  0.421286  2.191720  1.0730200
    ##  7 ENSMUSG00000000056  5.9998100  6.047990  6.020250  5.620120  6.4116300
    ##  8 ENSMUSG00000000058  5.2784700  3.971810  6.161450  7.045910  5.2136600
    ##  9 ENSMUSG00000000078 37.2141000 32.303000 31.249600  4.260570  3.5634100
    ## 10 ENSMUSG00000000085 25.4044000 22.950700 21.415600 16.825400 17.9712000
    ## # ... with 38,818 more rows, and 7 more variables: sample6 <dbl>,
    ## #   sample7 <dbl>, sample8 <dbl>, sample9 <dbl>, sample10 <dbl>,
    ## #   sample11 <dbl>, sample12 <dbl>

See how [R](https://www.r-project.org/) only prints 10 rows instead of returning all 38k?

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

See how both functions refer to `mov_expression` directly and not in quotations? That's *non-standard evaluation*. It makes code easier to read.

dplyr
=====

The most useful tool in the [tidyverse](http://tidyverse.org/) is [dplyr](http://dplyr.tidyverse.org/). It's a swiss-army knife for data manipulation. [dplyr](http://dplyr.tidyverse.org/) has 5 core functions that we recommend incorporating into your analysis:

-   `mutate()` adds new variables that are functions of existing variables
-   `select()` picks variables based on their names.
-   `filter()` picks cases based on their values.
-   `summarise()` reduces multiple values down to a single summary.
-   `arrange()` changes the ordering of the rows.

------------------------------------------------------------------------

Notes
=====

Programming
-----------

Underneath the hood, [tidyverse](http://tidyverse.org/) packages build upon the base [R](https://www.r-project.org/) language using [rlang](https://github.com/tidyverse/rlang/), which is a **complete rework** of how functions handle variable names and evaluate arguments. This is achieved through the `tidyeval` framework, which interprates command operations using `tidy evaluation`. This is outside of the scope of the course, but explained in detail in the [Programming with dplyr](http://dplyr.tidyverse.org/articles/programming.html) vignette, in case you'd like to understand how these new tools behave differently from base [R](https://www.r-project.org/).

Row names
---------

*Important*: [tidyverse](http://tidyverse.org/) is very opininationed about row names. These packages insist that all column data (e.g. `data.frame`) be treated equally, and that special designation of a column as `rownames` should be deprecated. [tibble](http://tibble.tidyverse.org/) provides simple utility functions to to handle rownames: `rownames_to_column()` and `columns_to_rowname()`.

Additional resources
====================

-   [R for Data Science](http://r4ds.had.co.nz)
-   [teach the tidyverse](http://varianceexplained.org/r/teach-tidyverse/)
-   [tidy style guide](http://style.tidyverse.org/)
