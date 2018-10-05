## Using R on a Unix system

For many analyses using R tools, the ability to utilize the resources of the cluster can greatly improve the speed and allocate greater memory to perform more efficient analyses, particularly for steps in single-cell RNA-seq (clustering, marker identification, and DE analysis) and ChIP-seq (ChIPQC and DiffBind) analyses. To run these analyses on the O2 cluster requires a set-up of a personal R library and X11 forwarding, if you want to interactively visualize output plots. 

Any non-base R packages needed for an analysis need to be downloaded to a personal R library prior to use. To create a personal R library, we can create a special directory. More information on personal R libraries on O2 are available on the [O2 Wiki](https://wiki.rc.hms.harvard.edu/display/O2/Personal+R+Packages).

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

It is possible to have multiple R libraries if you need to use different versions of R - perhaps you performed an analysis a year ago and want to analyze it a bit more with the same versions of all the tools. You could just comment out your current library and point your `.Renviron` to the older library:

```bash
DO NOT RUN!

# R_LIBS_USER="~/R/3.5.1/library"

R_LIBS_USER="~/R/3.4.1/library"
```

This can help with making your research more reproducible.

Once we have are library created, we can load the R module of appropriate version for our library and start R.

```bash
$ module load R/3.5.1

$ R
```

The terminal window should now turn into the R console with the R prompt `>`. 

Let's explore how to install packages by installing `tidyverse` and `Seurat`.

ADD INSTALLATIONS OF PACKAGES

If you were to manually install a package on Orchestra from CRAN we would have to specify where our library is using the following: `install.packages("name-of-your-package", lib="~/R/library")`. For Bioconductor packages nothing would change since we have already modified the environment variable to point to the library.You can also run R scripts from the command prompt in Unix. These scripts are just like shell scripts, but with R code in them; we created a few last session. For running a script from the Unix command prompt, it will have to take into account the absolute or relative location of the files and folders that will be used. Also, your local environment will need to have all the packages installed and available. 

You can run an R script from the shell command prompt in several ways, 3 different ways are listed below for a script called `mean.R`:
**Do not run this**
	
```bash
$ R < mean.R
	
$ R CMD BATCH mean.R
	
$ Rscript mean.R
```

### R on Orchestra:

R is available on Orchestra, and you can do all of the things we did on our laptops on the cluster instead. Let's try this out:

```bash
$ module avail stats/R
	
$ module load stats/R/3.2.5
	
$ R
```

As you can see, various versions of R are available on Orchestra, but there is no RStudio-like GUI. You can quit R and get back to the `$` command prompt by typing `q()`, no need to save the workspace image.
	
You can use any of the above ways to run an Rscript on Orchestra. But, you will need a different shebang line:

```bash
#!/usr/bin/env Rscript
```
And, you can also submit it as a job to the LSF queue as follows:

```bash
$ bsub -q short -W 12:00 -R "rusage[mem=16000]" "Rscript mean.R" 
# note the high memory usage above
```

Talk to the folks at HMS RC to find out which packages are already installed, and also about the best way to install R packages locally. They have a [how-to guide available online](https://wiki.med.harvard.edu/Orchestra/PersonalRPackages) for installing packages locally, if you feel comfortable trying it on your own.

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson was derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
