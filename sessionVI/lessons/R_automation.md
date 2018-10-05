## Using R on a Unix system

Sleuth is an R package, and while some R packages are automatically available to us on Orchestra, some of the packages we need to run Sleuth are not. Therefore, to run Sleuth on Orchestra, we need to manually install these programs into our personal R library. If you haven't created a personal R library, you can do so by entering the following code ([Orchestra Wiki](https://wiki.med.harvard.edu/Orchestra/WebHome)):

```bash
$ mkdir -p ~/R/library
```

### Installing R packages

Since Sleuth was designed to use the output of Kallisto as input, our Salmon transcript abundance estimates need to be massaged into the format of the Kallisto output. To do this, we are going to use the package [Wasabi](https://github.com/COMBINE-lab/wasabi). 

We have installed all of these packages for you to copy to your personal libraries:

```bash

$ export R_LIBS_USER="/home/mp298/R/library"
```

Now, start R:

```bash
$ R
```

The terminal window should now turn into the R console with the R prompt `>`. 

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
