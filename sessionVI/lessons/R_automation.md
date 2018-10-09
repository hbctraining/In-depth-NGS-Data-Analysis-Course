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
$ module load gcc/6.2.0 R/3.5.1

$ R
```

The terminal window should now turn into the R console with the R prompt `>`. You can run all of the analyses performed on our laptops on the cluster, but there is no RStudio-like GUI.

As we know, to do many of the analyses performed throughout the course requires many different packages. To install packages on O2 is a bit different than installing on our laptops. Let's explore how to install packages by considering `dplyr`, `Seurat` and `DESeq2`. 

To manually install a package on O2 from **CRAN**, we would need to specify where our library is using the following: `install.packages("name-of-your-package", lib="~/R/3.5.1/library")`. 

For instance, for installing `dplyr`, we would run the following code:

```r
# DO NOT RUN
install.packages("dplyr", lib="~/R/3.5.1/library")
```

> **NOTE:** You will be prompted to choose a CRAN mirror or server to download from - try to pick a close location.

Similar code would be run to install Seurat:

```r
# DO NOT RUN
install.packages("Seurat", lib="~/R/3.5.1/library")
```

However, for **Bioconductor** packages we do not need to specify the library path since we have already modified the environment variable to point to the library. 

```r
# DO NOT RUN
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```

You can quit R and get back to the `$` command prompt by typing `CTL + C`, and there is no need to save the workspace image.

## R scripts

In addition to running R interactively on the cluster, you can also run R scripts from the command prompt in Unix. These scripts are just like shell scripts, but with R code in them; we created a few in the past sessions. For running a script from the Unix command prompt, it will have to take into account the absolute or relative location of the files and folders that will be used. Also, your local environment will need to have all the packages installed and available. 

Let's explore R scripts in a bit more detail by using the script we created for the single-cell RNA-seq marker identification. Let's make a new folder in our home directory called `Rscripts` and copy over the `marker_id.R` script and the `pbmcs_seurat_tsne.rds` .RData object.

```bash
# Create directory for lesson
$ mkdir ~/Rscripts

# Change directories into Rscripts
$ cd ~/Rscripts

# Copy over the marker identification script
$ cp /n/groups/hbctraining/ngs-data-analysis-longcourse/sessionVI/marker_id.R .

# Copy over the TSNE data
$ cp /n/groups/hbctraining/ngs-data-analysis-longcourse/sessionVI/pbmcs_seurat_tsne.rds .
```

Now let's open the R script to edit with vim:

```bash
$ vim marker_id.R
```


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
