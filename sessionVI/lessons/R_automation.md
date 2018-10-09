## Using R on a Unix system

For many analyses using R tools, the ability to utilize the resources of the cluster can greatly improve the speed and allocate greater memory to perform more efficient analyses, particularly for steps in single-cell RNA-seq (clustering, marker identification, and DE analysis) and ChIP-seq (ChIPQC and DiffBind) analyses. To run these analyses on the O2 cluster requires a set-up of a personal R library and X11 forwarding, if you want to interactively visualize output plots. 

Any non-base R packages needed for an analysis need to be downloaded to a personal R library prior to use. To create a personal R library, we can create a special directory. More information on personal R libraries on O2 are available on the [O2 Wiki](https://wiki.rc.hms.harvard.edu/display/O2/Personal+R+Packages).

Let's create our personal R library, if not already created, in our `HOME` directory. Be sure to include the R version number for any R library created.

```bash
$ mkdir -p ~/R/3.5.1/library
```

### Installing R packages

To add packages to our personal R library, we should have the path to our library designated in a special file in our home directory that R explores for settings/information during every R session, called `.Renviron`. Let's create/open this file:

```bash
vim ~/.Renviron
```

```bash
R_LIBS_USER="~/R/3.5.1/library"
```

It is possible to have multiple R libraries if you need to use different versions of R, for example, perhaps you performed an analysis a year ago and want to analyze it a bit more with the same versions of all the tools. You could just comment out your current library and point your `.Renviron` to the older library:

```bash
DO NOT RUN!

# R_LIBS_USER="~/R/3.5.1/library"

R_LIBS_USER="~/R/3.4.1/library"
```

This can help with making your research more reproducible.

## R on O2:

Once we have our library created, we can load the R module appropriate for the version of our R library and start R. 

```bash
$ module spider R
```

There are various versions of R available on O2.

```bash
$ module load R/3.5.1

$ R
```

The terminal window should now turn into the R console with the R prompt `>`. You can run all of the analyses performed on our laptops on the cluster, but there is no RStudio-like GUI.

Let's explore how to install packages by installing `tidyverse`, `Seurat` and `DESeq2`. 

To manually install a package on O2 from **CRAN**, we need to specify where our library is using the following: `install.packages("name-of-your-package", lib="~/R/3.5.1/library")`. 

Let's start with installing tidyverse:

```r
install.packages("tidyverse", lib="~/R/3.5.1/library")
```

Now, for Seurat:

```r
install.packages("Seurat", lib="~/R/3.5.1/library")
```

For **Bioconductor** packages nothing would change since we have already modified the environment variable to point to the library. 

```r
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```

You can quit R and get back to the `$` command prompt by typing `CTL + C`, and there is no need to save the workspace image.

## R scripts

In addition to running R interactively on the cluster, you can also run R scripts from the command prompt in Unix. These scripts are just like shell scripts, but with R code in them; we created a few in the past sessions. For running a script from the Unix command prompt, it will have to take into account the absolute or relative location of the files and folders that will be used. Also, your local environment will need to have all the packages installed and available. 

You can run an R script from the shell command prompt in several ways, 3 different ways are listed below for a script called `mean.R`:
**Do not run this**
	
```bash
$ R < mean.R
	
$ R CMD BATCH mean.R
	
$ Rscript mean.R
```
	
You can use any of the above ways to run an Rscript on O2. But, you will need a different shebang line:

```bash
#!/usr/bin/env Rscript
```

And, you can also submit it as a job to a Slurm queue as follows:

```bash
$ sbatch -q short -W 12:00 -R "rusage[mem=16000]" "Rscript mean.R" 
# note the high memory usage above
```

Talk to the folks at HMS RC to find out which packages are already installed. They have a [how-to guide available online]() for installing packages locally, if you feel comfortable trying it on your own.

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
