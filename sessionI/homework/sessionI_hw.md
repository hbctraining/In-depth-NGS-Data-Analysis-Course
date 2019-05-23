# Session I: UNIX, HPC and FastQC Homework


 

## Logging in and starting an interactive job

1. Log into O2 and start a single core interactive job.


## Practicing Unix commands and using `vim`

2. Make a new directory called `ngs-course-homework` inside your `home` directory (~). Change directories to your `ngs-course-homework` directory.

* What is the absolute path to your `ngs-course-homework` directory?
 
 
3. Using vim, create a file entitled `mov10_fastq_reads.txt`.  On the first line type the sentence “Read counts from irrelevant siRNA, Mov10 knock down and Mov10 over-expression samples.” Save and exit vim.


4. Change directories into the `/n/groups/hbctraining/ngs-data-analysis-longcourse` directory. Locate and enter the `rnaseq` folder, then enter the `other` folder.

* What is the absolute path of your current directory?
* Copy both of the Mov10_kd files to the `ngs-course-homework` directory inside of your `home` directory. 
 

5. Change directories to your home directory. Locate and enter the `raw_data` folder (*Hint: this was from the rna-seq portion of Session I*).

* What is the absolute path of your current directory?
* Copy all of the fastq files from your current directory to your `ngs-course-homework` directory. You should now have 8 fastq files in your `ngs-course-homework` directory. 

6. Change directories to your `ngs-course-homework` directory. Pipe together the `grep` and `wc` commands to determine the number of reads generated for `Mov10_oe_1.subset.fastq`.  Look at your FastQC report for this sample and examine the statistics to determine whether your pipe command generated the correct number of reads for the sample. If not, fix your pipe command so that it gives you the same number of reads as the FastQC.

* How many reads were obtained for this sample?
* What was the command you used for determining the number of reads in your file?


## Practicing more advanced commands, including for_loops, redirection, scripts.

7. Using the pipe command from the previous question, write a `for_loop` to determine the number of reads obtained for every fastq sample in the `ngs-course-homework` folder.

* Within the for_loop, echo the filename, then, use your pipe command from Question #6 to determine the read counts for each file.

Your output should be in the format:

```bash
Mov10_oe_1.fastq

100

Mov10_oe_2.fastq

50…

```

8. Open vim and write a new script called `fastq_read_counts.sh`. This script will contain the `for_loop` you created in the previous question.

 * How many reads were generated for each of the files?
 

9. You should have seen the results of your script displayed in your shell. Now open your `fastq_read_counts.sh` script in vim. Edit your script, such that for every file you append the results to the file you created in Question #3 of this homework: `mov10_fastq_reads.txt`.

**Hint**: *Use “>>” to redirect the output. What does redirecting with “>>” do (google is your friend)?*

* Run your script for all fastq files in your `ngs-course-homework` directory. Your output in the `mov10_fastq_reads.txt` should look something like this:

```bash
Read counts from irrelevant siRNA, Mov10 knock down and Mov10 over-expression samples.

Mov10_oe_1.fastq

100

Mov10_oe_2.fastq

50…
```

**Upload the `mov10_fastq_reads.txt` script.**

## Practicing FastQC

10. Exit the interactive session (compute node) and go back to the login node.

* Make a copy of the mov10_fastqc.run script created during Session I, into your `ngs-course-homework` directory and give it the name hw_fastqc.run. Modify this newly copied script to perform FastQC on the following samples **in parallel**:

Irrel_kd_2.subset.fq, Mov10_kd_2.subset.fq, Mov10_oe_2.subset.fq

* Use Filezilla to download the zipped files to your own computer and examine the FastQC reports.
* Which samples have quality scores that drop below a score of 25?
* Do any samples have adapter contamination based on the information in the FastQC report?

**Upload the `hw_fastqc.run` script.**

## Metadata

11. Enter the `rnaseq` folder, and create a README file (hint: use vim to create the file), as described in the Exercise section of the [Data Management lesson](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/01_data_organization.html). Give a short description of the project and brief descriptions of the types of files you will be storing within each of the sub-directories.

12. Complete the metadata table in the [Experimental considerations lesson](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/experimental_planning_considerations.html), Question #1. What group did you assign to each sample?

13. Upload your answers to the questions on this assignment.
