# Session II: Abundance estimation and R Homework 

### Exploring the STAR parameters when creating an index using the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf), investigate the parameters of the STAR command to answer the questions following each STAR command script.

**Q1.** Below is the script we ran in class to create the **genome index** for alignment with STAR. We have added line numbers to the script to answer the questions that follow (you will need to refer to the online manual for STAR to answer these questions):
    
```
(1)  #!/bin/bash

(2)  #SBATCH -p priority # partition name
(3)  #SBATCH -t 0-2:00 # days-hours:minutes runlimit after which job will be killed.
(4)  #SBATCH -c 6 # number of cores requested
(5)  #SBATCH --job-name STAR_index         # Job name
(6)  #SBATCH -o %j.out       # File to which standard out will be written
(7)  #SBATCH -e %j.err       # File to which standard err will be written

(8)  STAR --runThreadN 6 \
(9)  --runMode genomeGenerate \
(10) --genomeDir my_genome_index \
(11) --genomeFastaFiles chr1.fa \
(12) --sjdbGTFfile chr1-hg19_genes.gtf \
(13) --sjdbOverhang 99
```

a. Which line number(s) in the script would you change if you were aligning your reads to the entire genome (instead of chr1)?

b. Provide the modified lines you would use in the script to use reference files from `/groups/shared_databases` on O2. (Hint: see [markdown](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/03_QC_STAR_and_Qualimap_run.html) for path information)

c. Could you run this indexing/alignment in STAR if you had no GTF annotation file? If so, which line(s) would change, and how would the line(s) change? 

d. Which line(s) in the script would you change if the length of your reads were 150nt? How would the line(s) change? 


### Exploring the Salmon quasi-alignment parameters

**Q2.** Below is the script containing the Salmon command we ran in class to **obtain abundance estimates** for our sequencing reads to the genome. We have added line numbers to the script to answer the questions that follow (you may need to refer to the online manual for [Salmon](https://salmon.readthedocs.io/en/latest/) to answer these questions):

```
(1)  #!/bin/bash
(2)  #SBATCH -p priority       # Partition to submit
(3)  #SBATCH -c 6                  # Number of cores
(4)  #SBATCH  -t 0-1:30               # Runtime in D-HH:MM (or use minutes)
(5)  #SBATCH  —job-name  salmon_mov10         # Job name
(6)  #SBATCH -o %j.out       # File to which standard out will be written
(7)  #SBATCH -e %j.err       # File to which standard err will be written

(8) salmon quant -i /n/groups/hbctraining/ngs-data-analysis-longcourse/rnaseq/salmon.ensembl38.idx \
(9)  -l A \
(10) -r ~/rnaseq/raw_data/Mov10_oe_1.subset.fq \
(11) -o Mov10_oe_1.subset.salmon \
(12) -p 6 \
(13) —writeMappings=salmon.out \
(14) —seqBias \
(15) —useVBOpt
```
    
a. Which line(s) in the script would you change if you wanted to change the number of cores you use to use only 3 cores? How would the line(s) change? 
    
b. Would any line(s) in the script change if you used an unstranded library preparation method for sequencing? If so, how would the line(s) change? 
    
c. Would any line(s) in the script change if you had paired-end data? If so, which line(s) would change, and how would the line(s) change?  
    
d. After performing QC you realized that you have GC bias present in your data. Is there anything you would change to address this when mapping and quantifying with Salmon?
    

### Positional parameters for scripting

**Q3.** Use vim to create a script by copying and pasting the contents of [this page](https://steve-parker.org/sh/eg/var3.sh.txt). Call it `pos_param_test.sh`.

Run the new script as follows and report back on what the contents of the variable `$0`, `$1`, `$2` and `$@` are in each case.
    
a. `sh pos_param_test.sh ngs_course 2018 "4" 6`
    
b. `sh pos_param_test.sh this is easy to understand`
    
c. `sh pos_param_test.sh 3 x 5 = 15 , 5 x 5 = 25`
    
**Q4.** Open up vim to create a shell script called run_salmon_single sample.sh. Add a shebang line to the top of your script. Copy and paste the Salmon command we used in class into this script when running it on `Mov10_oe_1.subset.fq`. You should have something like what is shown below:

```
#!/bin/bash
salmon quant -i /n/groups/hbctraining/ngs-data-analysis-longcourse/rnaseq/salmon.ensembl38.idx \
-l A \
-r ~/rnaseq/raw_data/Mov10_oe_1.subset.fq \
-o Mov10_oe_1.subset.salmon \
--writeMappings=salmon.out \
--seqBias \
--useVBOpt 
```
     
Now make the following modifications:
    
   - Remove the argument that writes the mappings to a SAM file
   - Before the Salmon command add a line to create a variable called fq and assign it the value of $1 (the filename that we expect the user to provide)
   - Change your command in the appropriate place(s) to now accept the value from this positional parameter as input to Salmon (which is now stored in the fq variable) **Hint: be sure to reference the variable with a $**
   - Use the fq variable to create a new variable called base that stores a prefix for the output. Change the Salmon command to utilize the base variable for naming out output directory
   - Save and exit vim. Upload this script.
    
**Q5.** From the homework in Session I, you should have the Mov10 knockdown 2 and Mov10 knockdown 3 (`Mov10_kd_2` and `Mov10_kd_3`) sample FASTQ files in your homework directory. Use the shell script you created in #4 to run Salmon on each of these samples.

### Introduction to R practice

#### Creating vectors/factors and dataframes

**Q6.** We are performing RNA-Seq on cancer samples being treated with three different types of treatment (A, B, and P). You have 12 samples total, with 4 replicates per treatment. Write the R code you would use to construct your metadata table as described below.

   - Create the vectors/factors for each column (Hint: you can type out each vector/factor, or if you want the process go faster try exploring the `rep()` function).
   - Put them together into a dataframe called `meta`.
   - Use the `rownames()` function to assign row names to the dataframe (Hint: you can type out the row names as a vector, or if you want the process go faster try exploring the `paste()` function).

Your finished metadata table should have information for the variables `sex`, `stage`, `treatment`, and `myc` levels: 

| |sex	| stage	| treatment	| myc |
|:--:|:--: | :--:	| :------:	| :--: |
|sample1|	M	|I	|A	|2343|
|sample2|	F	|II	|A	|457|
|sample3	|M	|II	|A	|4593|
|sample4	|F	|I	|A	|9035|
|sample5|	M	|II	|B	|3450|
|sample6|	F|	II|	B|	3524|
|sample7|	M|	I|	B|	958|
|sample8|	F|	II|	B|	1053|
|sample9|	M|	II|	P|	8674|
|sample10	|F|	I	|P	|3424|
|sample11|	M	|II	|P	|463|
|sample12|	F|	II|	P|	5105|


#### Subsetting vectors/factors and dataframes

**Q7.** Using the `meta` data frame from question #1, write out the R code you would use to perform the following operations (questions **DO NOT** build upon each other):

   - return only the `treatment` and `sex` columns using `[]`:
   - return the `treatment` values for samples 5, 7, 9, and 10 using `[]`:
   - use `filter()` to return all data for those samples receiving treatment `P`:
   - use `filter()`/`select()`to return only the `stage` and `treatment` columns for those samples with `myc` > 5000:
   - remove the `treatment` column from the dataset using `[]`:
   - remove samples 7, 8 and 9 from the dataset using `[]`:
   - keep only samples 1-6 using `[]`:
   - add a column called `pre_treatment` to the beginning of the dataframe with the values T, F, F, F, T, T, F, T, F, F, T, T (Hint: use `cbind()`): 
   - change the names of the columns to: "A", "B", "C", "D":

#### Releveling Factors

**Q8.** Use the `samplegroup` factor we created in class, and relevel it such that KO is the first level followed by CTL and OE.

#### More practice with R and "Nesting" (optional)

Let's derive some nested functions similar to those we will use in our RNA-Seq analysis. The following dataframes, `value_table` and `meta`, should be used to address the questions below (you do not actually need to create these dataframes):

**value_table**

| |MX1|	MX2|	MX3|
|:--: |:--:|	:--:|	:--:|
|KD.2	|-222517.197	|-21756.82	|-16036.035|
|KD.3	|17453.907	|-30058.14	|-25837.482|
|OE.1	|-31247.923|	73061.38	|7019.940|
|OE.2	|-4184.355	|61994.47	|1777.858|
|OE.3|	147391.709	|11970.45	|-18663.686|
|IR.1|	-32247.617	|-27896.01	|29383.153|
|IR.2	|25456.820|	-30714.29	|19148.752|
|IR.3	|99894.656|	-36601.04|	3207.501|

**meta**

| |sampletype|	MOVexpr|
|:--: |:--:|	:--:|
|KD.2|	MOV10_knockdown	|low|
|KD.3	|MOV10_knockdown|	low|
|OE.1	|MOV10_overexpression	|high|
|OE.2|	MOV10_overexpression|	high|
|OE.3	|MOV10_overexpression	|high|
|IR.1	|siRNA|	normal|
|IR.2	|siRNA|	normal|
|IR.3|	siRNA	|normal|


**Q9.** We would like to count the number of samples which have normal Mov10 expression (`MOVexpr`) in the `meta` dataset. Let's do this in steps:
 
   - Write the R code you would run to return the row numbers of the samples with `MOVexpr` equal to "normal".

   - Write the R code you would run to determine the number of elements in the `MOVexpr` column.

   - Now, try to combine your first two actions into a single line of code using nested functions to determine the number of elements in the MOVexpr column with expression levels of MOV10 being normal. 

**Q10.** Now we would like to add the `MX1` and `MX3` columns to the `meta` data frame. Let's do this in steps:

   - Write the R code you would run to extract columns `MX1` and `MX3` from the `value_table` and to save it to a variable `mx` (hint: you will need to use the `c()` function to specify the columns you want to extract). 

   - Using the `cbind()` function, write the R code you would use to add the columns in your `mx` variable to the end of your `meta` dataset.
 
   - Now, try to combine your first two actions into a single line of code using nested functions (hint: you do not need to generate the `mx` variable) to add the `MX1` and `MX3` columns to the `meta` file.
   
