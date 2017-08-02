Introduction to knitr
================
Michael J. Steinbaugh
2017-07-07

``` r
library(knitr)
```

[knitr](https://yihui.name/knitr/), developed by [Yihui Xie](https://yihui.name), is an [R](https://www.r-project.org/) package designed for report generation within [RStudio](https://www.rstudio.com/). It enables dynamic generation of multiple file formats from an [RMarkdown](http://rmarkdown.rstudio.com/) file, including HTML and PDF documents. As [RMarkdown](http://rmarkdown.rstudio.com/) grows as an acceptable [reproducible manuscript](https://elifesciences.org/labs/cad57bcf/composing-reproducible-manuscripts-using-r-markdown) format, using [knitr](https://yihui.name/knitr/) to generate a report summary is becoming common practice. Knit report generation is now integrated into [RStudio](https://www.rstudio.com/), and can be accessed using the GUI or console.

Chunks
======

The basic idea of [knitr](https://yihui.name/knitr/) (along with [RMarkdown](http://rmarkdown.rstudio.com/)) is that you can write your analysis workflow in plain text [Markdown](https://daringfireball.net/projects/markdown/syntax), and intersperse chunks of code delimited with a special marker (\`\`\`). Backticks (\`) commonly indicate code and are also used on [GitHub](https://github.com). Each chunk should be given a **unique** name. [knitr](https://yihui.name/knitr/) isn't very picky how you name the code chunks, but I recommend using `snake_case` for the names whenever possible. The output of each chunk can be cusotmized, using a number of options described below.

Additionally, you can write inline [R](https://www.r-project.org/) code enclosed by single backticks (\`) containing a lowercase `r` (like \`\`\` code chunks). This allows for variable returns outside of code chunks, and is extremely useful for making report text more dynamic. For example, you can print the current date inline with this syntax: `` ` r Sys.Date() ` `` (no spaces).

Code languages
==============

Recently [RStudio](https://www.rstudio.com/) has added support for additional [code languages](http://rmarkdown.rstudio.com/lesson-5.html) besides [R](https://www.r-project.org/):

-   Python
-   SQL
-   Bash
-   Rcpp
-   Stan
-   JavaScript
-   CSS

This is still in development but allows for integration of multiple programming languages into a single workflow, and is shaping up to be a very powerful tool for NGS analysis in the near future.

Global options
==============

[knitr](https://yihui.name/knitr/) allows for global options to be set on all chunks in an [RMarkdown](http://rmarkdown.rstudio.com/) file, and I recommend these be placed inside your `setup` chunk.

``` r
opts_chunk$set(
    autodep = TRUE,
    cache = TRUE,
    cache.lazy = TRUE,
    dev = c("png", "pdf", "svg"),
    error = TRUE,
    fig.height = 6,
    fig.retina = 2,
    fig.width = 6,
    highlight = TRUE,
    message = FALSE,
    prompt = TRUE,
    tidy = TRUE,
    warning = FALSE)
```

An additional cool trick is that you can save `opts_chunk$set` settings in `~/.Rprofile` and these [knitr](https://yihui.name/knitr/) options will apply to all of your [RMarkdown](http://rmarkdown.rstudio.com/) documents.

The setup chunk
===============

The `setup` chunk is a special knitr chunk that should be placed at the start of the document. We recommend storing all `library()` loads required for the script and other `load()` requests for external files here. In our RMarkdown templates, such as the bcbioRnaseq [differential expression template](https://github.com/hbc/bcbioRnaseq/blob/master/inst/rmarkdown/templates/differential_expression/skeleton/skeleton.Rmd), we store all the user-defined parameters in the `setup` chunk that are required for successful knitting.

    {r setup, include=FALSE}
    knitr::opts_chunk$set(echo = TRUE)

Per chunk options
=================

    {r chunk_name, options}

[knitr](https://yihui.name/knitr/) provides a lot of customization options, and this can be overwhelming at first. Here's a short list of some options I commonly use in chunks:

-   `echo = FALSE`
-   `eval = FALSE`
-   `include = FALSE`
-   `message = FALSE`
-   `results = "asis"`
-   `warning = FALSE`

In particular, you can easily resize images:

-   `fig.height = 6`
-   `fig.width = 4`

Figures
=======

My favorite feature of [knitr](https://yihui.name/knitr/) is how much simpler it makes generating figures. You can simply return a plot in a chunk, and [knitr](https://yihui.name/knitr/) will automatically write the files to disk, in an organized subfolder. By specifying options in a `setup` chunk (below), you can have R automatically save your plots in multiple file formats at once, including PNG, PDF, and SVG. A single chunk can support multiple plots, and they will be arranged in squares below the chunk in [RStudio](https://www.rstudio.com/).

Tables
======

[knitr](https://yihui.name/knitr/) includes a simple but powerful function for generating stylish tables in a knit report named `kable()`. Here's an example using [R](https://www.r-project.org/)'s built-in `mtcars` dataset:

``` r
help("kable", "knitr")
mtcars %>%
    head %>%
    kable
```

|                   |   mpg|  cyl|  disp|   hp|  drat|     wt|   qsec|   vs|   am|  gear|  carb|
|-------------------|-----:|----:|-----:|----:|-----:|------:|------:|----:|----:|-----:|-----:|
| Mazda RX4         |  21.0|    6|   160|  110|  3.90|  2.620|  16.46|    0|    1|     4|     4|
| Mazda RX4 Wag     |  21.0|    6|   160|  110|  3.90|  2.875|  17.02|    0|    1|     4|     4|
| Datsun 710        |  22.8|    4|   108|   93|  3.85|  2.320|  18.61|    1|    1|     4|     1|
| Hornet 4 Drive    |  21.4|    6|   258|  110|  3.08|  3.215|  19.44|    1|    0|     3|     1|
| Hornet Sportabout |  18.7|    8|   360|  175|  3.15|  3.440|  17.02|    0|    0|     3|     2|
| Valiant           |  18.1|    6|   225|  105|  2.76|  3.460|  20.22|    1|    0|     3|     1|

There are some other functions that allow for more powerful customization of tables, including `pander::pander()` and `xtable::xtable()`, but I generally perfer the simplicity and cross-platform reliability of `knitr::kable()`.

Generating the report
=====================

`knit()` (recommended)
----------------------

``` r
help("knit", "knitr")
```

Once we've finished creating an [RMarkdown](http://rmarkdown.rstudio.com/) file containing code chunks, we finally need to knit the report. When executing `knit()` on a document, by default this will generate an HTML report. If you would prefer a different document format, this can be specified in the YAML header with the `output:` parameter:

-   html\_document
-   pdf\_document
-   github\_document

[RStudio](https://www.rstudio.com/) now supports a [number of formats](http://rmarkdown.rstudio.com/formats.html), each with their own customization options. Consult their webiste for more details.

`render()` (advanced)
---------------------

``` r
help("render", "rmarkdown")
```

The `knit()` command works great if you only need to generate a single document format. [RMarkdown](http://rmarkdown.rstudio.com/) also supports a more advanced function named `rmarkdown::render()`, allows for output of multiple document formats. To accomplish this, we recommend saving a special file named `_output.yaml` in your project root. Here's an [example](https://github.com/hbc/bcbioRnaseq/blob/master/docs/downloads/_output.yaml) from our [bcbioRnaseq](https://github.com/hbc/bcbioRnaseq) package:

    rmarkdown::html_document:
        code_folding: hide
        df_print: kable
        highlight: pygments
        number_sections: false
        toc: true
    rmarkdown::pdf_document:
        number_sections: false
        toc: true
        toc_depth: 1

**Note**: PDF rendering is sometimes problematic, especially when running [R](https://www.r-project.org/) remotely, like on the Orchestra cluster. If you run into problems, it's likely an issue related to [pandoc](http://pandoc.org).

Working directory behavior
==========================

**Note**: [knitr](https://yihui.name/knitr/) redefines the working directory of an RMarkdown file in a manner that can be confusing. If you're working in [RStudio](https://www.rstudio.com/) with an [RMarkdown](http://rmarkdown.rstudio.com/) file that is not at the same location as the current R working directory (`getwd()`), you can run into problems with broken file paths. Let's say I have [RStudio](https://www.rstudio.com/) open without a project loaded. My working directory is set to `~/Users/mike`. Now, if I load an RMarkdown file from my desktop at `~/Users/mike/Desktop`, [knitr](https://yihui.name/knitr/) will set the working directory within chunks to be relative to my desktop. We advise against coding paths in a script to only work with [knitr](https://yihui.name/knitr/) and not base [R](https://www.r-project.org/).

A simple way to resolve this issue is by creating an R project for the analysis, and saving all RMarkdown files at the top level, to avoid running into unexpected problems related to this behavior.

Convert an R script to an RMarkdown knit report
===============================================

Let's convert our Mov10 DE analysis script into an RMarkdown report:

``` r
download.file(
    file.path("https://raw.githubusercontent.com",
              "hbctraining",
              "In-depth-NGS-Data-Analysis-Course",
              "may2017",
              "sessionVI",
              "scripts",
              "de_script_to_knit.R"),
    destfile = "de_script_to_knit.R")
```

Open this script as a new tab in RStudio. Let's create a new RMarkdown called `de_analysis.Rmd`. You can do this quickly in RStudio by going to `File -> New File -> R Markdown...`. Customize the RMarkdown YAML header and create chunks for each code block.

Additional links
================

-   [knitr in a knutshell](http://kbroman.org/knitr_knutshell/)
-   [knitr book](https://www.amazon.com/gp/product/1498716962)
-   [knitr examples](https://yihui.name/knitr/demos)
-   [knitr vignettes](https://github.com/yihui/knitr/tree/master/vignettes)
