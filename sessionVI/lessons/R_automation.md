## Using R on a Unix system

For many analyses using R tools, the ability to utilize the resources of the cluster can greatly improve the speed and allocate greater memory to perform more efficient analyses, particularly for steps in single-cell RNA-seq (clustering, marker identification, and DE analysis) and ChIP-seq (ChIPQC and DiffBind) analyses. To run these analyses on the O2 cluster requires a set-up of a personal R library and X11 forwarding, if you want to interactively visualize output plots. 

Any non-base R packages needed for an analysis need to be downloaded to a personal R library prior to use. To create a personal R library, we can create a special directory. More information on personal R libraries on O2 are available on the [O2 Wiki](https://wiki.rc.hms.harvard.edu/display/O2/Personal+R+Packages).

Let's create our personal R library, if not already created, in our `HOME` directory. Be sure to include the R version number for any R library created.

```bash
$ srun --pty -p interactive -t 0-12:00 --mem 36G /bin/bash

$ mkdir -p ~/R/3.5.1/library
```

> **NOTE:** If you were using X11 forwarding to view images, you could include the `--x11` flag in the interactive command:
`srun --pty -p interactive -t 0-12:00 --x11 --mem 36G /bin/bash`


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
$ module load gcc/6.2.0 R/3.5.1 hdf5/1.10.1
```

"The `hdf5` package provides an interface between HDF5 and R. HDF5's main features are the ability to store and access very large and/or complex datasets and a wide variety of metadata on mass storage (disk) through a completely portable file format. The rhdf5 package is thus suited for the exchange of large and/or complex datasets between R and other software package, and for letting R applications work on datasets that are larger than the available RAM" [[Bioconductor hdf5 package discription](https://bioconductor.org/packages/release/bioc/html/rhdf5.html)].

```bash
$ R
```

The terminal window should now turn into the R console with the R prompt `>`. You can run all of the analyses performed on our laptops on the cluster, but there is no RStudio-like GUI.

### R package installations

As we know, to do many of the analyses performed throughout the course requires many different R packages. To install packages on O2 is a bit different than installing on our laptops. While we will explore how to install the packages `dplyr`, `Seurat`, `AnnotationHub`, and `ensembldb`, it takes quite a bit of time to install, so we encourage you to do this on your own later. 

To manually install a package on O2 from **CRAN**, we would need to specify where our library is using the following: `install.packages("name-of-your-package", lib="~/R/3.5.1/library")`. 

For instance, for installing `dplyr` and `Seurat`, we would run the following code:

```r
# DO NOT RUN
install.packages("dplyr", lib="~/R/3.5.1/library")

install.packages("Seurat", lib="~/R/3.5.1/library")
```

> **NOTE:** You will be prompted to choose a CRAN mirror or server to download from - try to pick a relatively close location. 
> 
> Often during installations you might encounter errors such as: 
>
> `Error in loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]) : there is no package called ‘R.utils’`.
>
> Whenever dealing with these messages try to install the packages mentioned (i.e. `R.utils`).

However, for **Bioconductor** packages we do not need to specify the library path since we have already modified the environment variable to point to the library. 

```r
# DO NOT RUN
source("https://bioconductor.org/biocLite.R")
biocLite("AnnotationHub")

biocLite("ensembldb")
```

If typing in R and wanting to use `esc` to return the command prompt, on the command line we need to use `CTL + C` instead.

You can quit R and get back to the `$` command prompt by typing `q()`, and there is no need to save the workspace image.

> **NOTE:** Talk to the folks at HMS RC to find out which packages are already installed for each of the R modules available. 

For the rest of this session, we are going to use an R library in which we have already installed these libraries, so we are going to add this library to our `.Renviron` file. 

```bash
$ vim ~/.Renviron
```

```bash
R_LIBS_USER="/n/groups/hbctraining/R/R-3.5.1"
```

## R scripts

In addition to running R interactively on the cluster, you can also run R scripts from the command prompt in Unix. These scripts are just like shell scripts, but with R code in them; we created a few in the past sessions. For running a script from the Unix command prompt, it will have to take into account the absolute or relative location of the files and folders that will be used. Also, your local environment (library specified in the `.Renviron` file) will need to have all the packages installed and available. 

Let's explore R scripts in a bit more detail by using part of the script we created for the **single-cell RNA-seq marker identification**. Let's make a new folder in our home directory called `Rscripts` and run code from our `marker_id.R` script using the `pbmcs_seurat_tsne.rds` .RData object.

```bash
# Create directory for lesson
$ mkdir ~/Rscripts

# Change directories into Rscripts
$ cd ~/Rscripts

# Copy over the TSNE data
$ cp /n/groups/hbctraining/ngs-data-analysis-longcourse/sessionVI/pbmcs_seurat_tsne.rds .
```

Now let's create the R script with vim. We have added `print()` statements to the script to help with trouble-shooting, if necessary:

```bash
$ vim marker_id.R
```

```r
# Single-cell RNA-seq - Marker identification

# Load libraries
library(dplyr)
library(Seurat)
library(AnnotationHub)
library(ensembldb)

# Load Seurat clustered data
seurat <- readRDS("pbmcs_seurat_tsne.rds")

print("Identifying markers")
# Identify gene markers
all_markers <-FindAllMarkers(seurat,
                             min.pct =  0.25,
                             min.diff.pct = 0.25)

print("Acquiring annotations")
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)
	      
# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)
  
# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")
		     
# Merge annotations with all markers		     
all_markers <- dplyr::left_join(all_markers, annotations[ , c(1:3, 5)],
                         by = c("gene" = "gene_id"))

# Rearrange order of columns to make clearer
all_markers <- all_markers[, c(6:8, 1:5, 9:10)]

print("Writing to file")
# Write results to file
write.csv(all_markers, "all_markers.csv", quote = F, row.names = FALSE)
```
	
Before we can run the script we don't want to forget the shebang line, which is slightly different since we are running an R script:

```bash
#!/usr/bin/env Rscript
```

You can run an R script from the shell command prompt in several ways, such as:
	
```bash
# DO NOT RUN!
$ R < marker_id.R
	
$ R CMD BATCH marker_id.R
```

But we are going to run our `marker_id.R` script using the following method:

```r
$ Rscript marker_id.R
```

Instead of running interactively, you could also submit it as a job to a Slurm queue as follows:

```bash
# DO NOT RUN!
$ sbatch -p priority -t 0-12:00 --job-name sc_marker_id --mem 36G -o %j.out -e %j.err "Rscript marker_id.R" 
# note the high memory usage above
```

> **NOTE:** It is not uncommon to require additional memory if needed, especially for single-cell RNA-seq (some projects have required up to 128G memory). Since requesting so much memory can take a while to receive, submitting the job to the `priority` queue can be helpful.

### R positional parameters

Finally, it is helpful to know that these R scripts can take positional parameters as well. Therefore, we could use the same script to process different clustering analyses. We can add command line arguments to the script using the `commandArgs()` function:

```bash
$ vim marker_id.R
```

Go to top of file using `gg` in command mode and **delete the entire file** with `dG`.

Now paste in the following script with command line arguments and save.

```r
#!/usr/bin/env Rscript

# Usage: this Rscript is using Seurat to identify cell cluster markers. The input is an RDS file containing a Seurat object with clustering information contained within, and the output is a csv file containing the cluster markers. The script expects as a command line argument the path to the Seurat object and prefix to output file. To run:  Rscript marker_id.R "path/to/seurat.rds" "prefix_to_output_file"

# Single-cell RNA-seq - Marker identification

# Load libraries
library(dplyr)
library(Seurat)
library(AnnotationHub)
library(ensembldb)

#options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

# Load Seurat clustered data
seurat <- readRDS(args[1])


print("Identifying markers")
# Identify gene markers
all_markers <-FindAllMarkers(seurat,
                             min.pct =  0.25,
                             min.diff.pct = 0.25)

print("Acquiring annotations")
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)
	      
# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)
  
# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")
		     
# Merge annotations with all markers		     
all_markers <- dplyr::left_join(all_markers, annotations[ , c(1:3, 5)],
                         by = c("gene" = "gene_id"))

# Rearrange order of columns to make clearer
all_markers <- all_markers[, c(6:8, 1:5, 9:10)]

print("Writing to file")
# Write results to file
write.csv(all_markers,  paste0(args[2], "_all_markers.csv"), quote = F, row.names = FALSE)
```

Now to run we can provide the parameters:

```bash
Rscript marker_id.R "pbmcs_seurat_tsne.rds" "seurat_res0.8_all_clusters"
```

The print statements will output as the script proceeds. When completed take a quick peek at the markers file generated:

```bash
less seurat_res0.8_all_clusters_all_markers.csv
```

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
