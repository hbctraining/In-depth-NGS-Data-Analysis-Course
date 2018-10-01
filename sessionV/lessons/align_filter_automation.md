---
title: "Automating QC and Alignment"
author: "Meeta Mistry, Radhika Khetani"
date: "March 20th, 2018"
---

Contributors: Meeta Mistry, Radhika Khetani

Approximate time: 80 minutes

## Learning Objectives

* Write a shell script to generate filtered BAM files for all samples
* Describe the difference between serial and parallel jobs

## Automating the ChIP-seq analysis path from sequence reads to BAM files

Once you have optimized tools and parameters using a single sample (using an interactive session), you can write a script to run the whole or a portion of the workflow on multiple samples in parallel (batch submission with a shell script).

Writing a reusable shell script ensures that every sample is run with the exact same parameters, and helps to keep track of all the tools and their versions. The shell script is like a lab notebook; in the future, you (or your colleagues) can go back and check the workflow for methods and versions, which goes a long way to making your work more efficient and reproducible.

Before we start with the script, let's check how many cores our interactive session has by using `sacct`. 

```bash
$ sacct
```

We need to have an interactive session with 6 cores, if you already have one you are set. If you have a session with fewer cores then `exit` out of your current interactive session and start a new one with `-n 6`.

```bash
$ srun --pty -p short -t 0-12:00 -n 6 --mem 8G --reservation=HBC /bin/bash
```

### More Flexibility with variables

We can write a shell script that will run on a specific file, but to make it more flexible and efficient we would prefer that it lets us give it an input fastq file when we run the script. To be able to provide an input to any shell script, we need to use **Positional Parameters**.

For example, we can refer to the components of the following command as numbered variables **within** the actual script:

```bash
# * DO NOT RUN *
sh  run_analysis.sh  input.fastq  input.gtf  12
```

`$0` => run_analysis.sh

`$1` => input.fastq

`$2` => input.gtf

`$3` => 12

The variables $1, $2, $3,...$9 and so on are **positional parameters** in the context of the shell script, and can be used within the script to refer to the files/number specified on the command line. Basically, the script is written with the expectation that $1 will be a fastq file and $2 will be a GTF file, and so on.

*There can be virtually unlimited numbers of inputs to a shell script, but it is wise to only have a few inputs to avoid errors and confusion when running a script that used positional parameters.*

> [This is an example of a simple script that used the concept of positional parameters and the associated variables](http://steve-parker.org/sh/eg/var3.sh.txt). You should try this script out after the class to get a better handle on positional parameters for shell scripting.

Let's use this new concept in the script we are writing. We want the first positional parameter ($1) to be the name of our fastq file. We could just use the variable `$1` throughout the script to refer to the fastq file, but this variable name is not intuitive, so we want to create a new variable called `fq` and copy the contents of `$1` into it.

First, we need to start a new script called `chipseq_analysis_on_input_file.sh` in the `~/chipseq/scripts/` directory:

```bash
$ cd ~/chipseq/scripts/

$ vim chipseq_analysis_on_input_file.sh
```

```bash
#!/bin/bash/

# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1
```

> When we set up variables we do not use the `$` before it, but when we *use the variable*, we always have to have the `$` before it. >
>
> For example: 
>
> initializing the `fq` variable => `fq=$1`
>
> using the `fq` variable => `fastqc $fq`

To ensure that all the output files from the workflow are properly named with sample IDs we should extract the "base name" (or sample ID) from the name of the input file.

```
# grab base of filename for naming outputs
base=`basename $fq _chr12.fastq`
echo "Sample name is $base"           
```

> **Remember `basename`?**
>
> 1. the `basename` command: this command takes a path or a name and trims away all the information before the last `\` and if you specify the string to clear away at the end, it will do that as well. In this case, if the variable `$fq` contains the path *"~/chipseq/raw_data/H1hesc_Nanog_Rep1_chr12.fastq"*, `basename $fq _chr12.fastq` will output "H1hesc_Nanog_Rep1".
> 2. to assign the value of the `basename` command to the `base` variable, we encapsulate the `basename...` command in backticks. This syntax is necessary for assigning the output of a command to a variable.

We'll create output directories, but with the `-p` option. This will make sure that `mkdir` will create the directory only if it does not exist, and it won't throw an error if it does exist.

```
# make all of the output directories
# The -p option means mkdir will create the whole path if it 
# does not exist and refrain from complaining if it does exist
mkdir -p ~/chipseq/results/fastqc
mkdir -p ~/chipseq/results/bowtie2/intermediate_bams
```

Now that we have already created our output directories, we can now create variables using those directories. 

```
# directory with bowtie genome index
genome=~/chipseq/reference_data/chr12

# set up output filenames and locations
fastqc_out=~/chipseq/results/fastqc/

## set up file names
align_out=~/chipseq/results/bowtie2/${base}_unsorted.sam
align_bam=~/chipseq/results/bowtie2/${base}_unsorted.bam
align_sorted=~/chipseq/results/bowtie2/${base}_sorted.bam
align_filtered=~/chipseq/results/bowtie2/${base}_aln.bam

## set up more variables for 2 additional directoties to help clean up the results folder
bowtie_results=~/chipseq/results/bowtie2
intermediate_bams=~/chipseq/results/bowtie2/intermediate_bams
```

Creating these variables makes it easier to see what is going on in a long command wherein we can use we can now use `align_out` instead of `~/chipseq/results/bowtie2/${base}_unsorted.sam`. In addition, if there is a need to change the output diretory or the genome being aligned to, the change needs to be made just in the one location instead of throughout the script. 

### Keeping track of tool versions

All of our variables are now staged. Next, let's make sure all the modules are loaded. This is also a good way to keep track of the versions of tools that you are using in the script:

```
# set up the software environment
module load fastqc/0.11.3
module load gcc/6.2.0  
module load bowtie2/2.2.9
module load samtools/1.3.1
export PATH=/n/app/bcbio/tools/bin:$PATH 	# for using 'sambamba'
```

### Preparing for future debugging

In the script, it is a good idea to use `echo` for debugging. `echo` basically displays the string of characters specified within the quotations. When you have strategically place `echo` commands specifying what stage of the analysis is next, in case of failure you can determine the last `echo` statement displayed to troubleshoot the script.

```
echo "Processing file $fq"
```

> You can also use `set -x`:
>
> `set -x` is a debugging tool that will make bash display the command before executing it. In case of an issue with the commands in the shell script, this type of debugging lets you quickly pinpoint the step that is throwing an error. Often, tools will display the error that caused the program to stop running, so keep this in mind for times when you are running into issues where this is not available.
> You can turn this functionality off by saying `set +x`

### Running the tools

Let's write up the commands to run the tools we have already tested, with a couple of modifications:
* use variable names instead of actual file names
* multithread when possible (`bowtie2`, `sambamba`)

```
# Run FastQC
fastqc $fq

# Run bowtie2
bowtie2 -p 6 -q --local -x $genome -U $fq -S $align_out

# Create BAM from SAM
samtools view -h -S -b -@ 6 -o $align_bam $align_out

# Sort BAM file by genomic coordinates
sambamba sort -t 6 -o $align_sorted $align_bam

# Filter out multi-mappers and duplicates
sambamba view -h -t 6 -f bam -F "[XS] == null and not unmapped and not duplicate" $align_sorted > $align_filtered

# Move intermediate files out of the bowtie2 directory
mv $bowtie_results/${base}*sorted* $intermediate_bams
```

### Last addition to the script

It is best practice to have the script **usage** specified at the top any script. This should have information such that when your future self, or a co-worker, uses the script they know what it will do and what input(s) are needed. For our script, we should have the following lines of comments right at the top after `#!/bin/bash/`:

```
# This script takes a fastq file of ChIP-seq data, runs FastQC and outputs a BAM file for it that is ready for peak calling. Bowtie2 is the aligner used, and the outputted BAM file is sorted by genomic coordinates and has duplicate reads removed using sambamba.
# USAGE: sh chipseq_analysis_on_input_file.sh <path to the fastq file>
```

It is okay to specify this after everything else is set up, since you will have most clarity about the script only once it is fully done.

Your script should now look like this:

```
#!/bin/bash/

# This script takes a fastq file of ChIP-seq data, runs FastQC and outputs a BAM file for it that is ready for peak calling. Bowtie2 is the aligner used, and the outputted BAM file is sorted by genomic coordinates and has multi-mappers and duplicate reads removed using sambamba.
# USAGE: sh chipseq_analysis_on_input_file.sh <path to the fastq file>

# initialize a variable with an intuitive name to store the name of the input fastq file
fq=$1

# grab base of filename for naming outputs
base=`basename $fq .fastq`
echo "Sample name is $base"    

# directory with bowtie genome index
genome=~/chipseq/reference_data/chr12

# make all of the output directories
# The -p option means mkdir will create the whole path if it 
# does not exist and refrain from complaining if it does exist
mkdir -p ~/chipseq/results/fastqc
mkdir -p ~/chipseq/results/bowtie2/intermediate_bams

# set up output filenames and locations
fastqc_out=~/chipseq/results/fastqc/

## set up file names
align_out=~/chipseq/results/bowtie2/${base}_unsorted.sam
align_bam=~/chipseq/results/bowtie2/${base}_unsorted.bam
align_sorted=~/chipseq/results/bowtie2/${base}_sorted.bam
align_filtered=~/chipseq/results/bowtie2/${base}_aln.bam

## set up more variables for 2 additional directoties to help clean up the results folder
bowtie_results=~/chipseq/results/bowtie2
intermediate_bams=~/chipseq/results/bowtie2/intermediate_bams

# set up the software environment
module load fastqc/0.11.3
module load gcc/6.2.0  
module load bowtie2/2.2.9
module load samtools/1.3.1
export PATH=/n/app/bcbio/tools/bin:$PATH 	# for using 'sambamba'

echo "Processing file $fq"

# Run FastQC and move output to the appropriate folder
fastqc $fq
mv *_fastqc.* ~/chipseq/results/fastqc/

# Run bowtie2
bowtie2 -p 6 -q --local -x $genome -U $fq -S $align_out

# Create BAM from SAM
samtools view -h -S -b -@ 6 -o $align_bam $align_out

# Sort BAM file by genomic coordinates
sambamba sort -t 6 -o $align_sorted $align_bam

# Filter out duplicates
sambamba view -h -t 6 -f bam -F "[XS] == null and not unmapped " $align_sorted > $align_filtered

# Move intermediate files out of the bowtie2 directory
mv $bowtie_results/${base}*sorted* $intermediate_bams
```

### Saving and running script

We should all have an interactive session with 6 cores, so we can run the script as follows:

```bash
$ sh chipseq_analysis_on_input_file.sh ~/chipseq/raw_data/H1hesc_Nanog_Rep1_chr12.fastq
```

## Submitting jobs **in serial** to the SLURM scheduler

The above script will run in an interactive session **one file at a time**. But the whole point of writing this script was to run it on all files at once. How do you think we can do this?

To run the above script **"in serial"** for all of the files on a worker node via the job scheduler, we can create a separate submission script that will need 2 components:

1. **SLURM directives** at the **beginning** of the script. This is so that the scheduler knows what resources we need in order to run our job on the compute node(s).
2. a **`for`** loop that iterates through and runs the above script for all the fastq files.

Below is what this second script (`chipseq_analysis_on_allfiles.slurm`) would look like **\[DO NOT RUN THIS\]**:

```
# **\[DO NOT RUN THIS\]**

#!/bin/bash

#SBATCH -p short 		# partition name
#SBATCH -t 0-2:00 		# hours:minutes runlimit after which job will be killed
#SBATCH -n 6 		# number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --job-name Nanog_Rep1 		# Job name
#SBATCH -o %j.out			# File to which standard out will be written
#SBATCH -e %j.err 		# File to which standard error will be written

# this `for` loop, will take chip-seq fastq files as input and output filtered BAM files ready for peak calling.

for fq in ~/chipseq/raw_data/*.fastq
do
  echo "running analysis on $fq"
  sh chipseq_analysis_on_input_file.sh $fq
done
```

**But we don't want to run the analysis on these 6 samples one after the other!** We want to run them "in parallel" as 6 separate jobs. 

> **Note:** If you create and run the above script, or something similar to it, i.e. with SLURM directives at the top, you should give the script name `.run` or `.slurm` as the extension. This will make it obvious that it is meant to submit jobs to the SLURM scheduler. 

To the run the above script you would have used the following command: `sbatch chipseq_analysis_on_allfiles.slurm`.  

## Submitting jobs **in parallel** to the SLURM scheduler

Parallelization will save you a lot of time with real (large) datasets. To parallelize our analysis, we will still need to write a second script that will call the original script we just wrote. We will still use a `for` loop, but we will be creating a regular shell script and we will be specifying the SLURM directives a little differently. 

Use `vim` to start a new shell script called `chipseq_analysis_on_allfiles-for_slurm.sh`: 

```bash
$ vim chipseq_analysis_on_allfiles_for-slurm.sh
```

This script loops through the same files as in the previous (demo) script, but the command being submitted within the `for` loop is `sbatch` with SLURM directives specified on the same line:

```bash
#! /bin/bash

for fq in ~/chipseq/raw_data/*.fastq
do

sbatch -p short -t 0-2:00 -n 6 --job-name chipseq-analysis -o %j.out -e %j.err \
--wrap="sh ~/chipseq/scripts/chipseq_analysis_on_input_file.sh $fq"

sleep 1	    # wait 1 second between each job submission
  
done
```
> Please note that after the `sbatch` directives the command `sh ~/chipseq/scripts/chipseq_analysis_on_input_file.sh $fq` is in quotes.

What you should see on the output of your screen would be the jobIDs that are returned from the scheduler for each of the jobs that your script submitted.

You can use `sacct login_ID` to check progress.

Don't forget about the `scancel` command, should something go wrong and you need to cancel your jobs.

> **NOTE:** All job schedulers are similar, but not the same. Once you understand how one works, you can transition to another one without too much trouble. They all have their pros and cons which are considered by the system administrators when picking one for a given HPC environment. 

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
