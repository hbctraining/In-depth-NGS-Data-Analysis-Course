---
title: "The Shell: Searching and Redirection"
author: Sheldon  McKay, Bob Freeman, Mary Piper, Michael Steinbaugh
date: 2017-05-31
duration: 45
---

- Learn how to search for characters or patterns in a text file using the `grep` command
- Learn about output redirection
- Learn how to use the pipe (`|`) character to chain together commands



# Searching files

We went over how to search within a file using `less`. We can also
search within files without even opening them, using `grep`. Grep is a command-line
utility for searching plain-text data sets for lines matching a pattern or regular expression (regex).
Let's give it a try!

We are going to practice searching with `grep` using our fastq files, which contain the sequencing reads (nucleotide sequences) output from a sequencing facility. Each sequencing read in a FASTQ file is associated with four lines of output, with the first line (header line) always starting with an `@` symbol. A whole fastq record for a single read should appear similar to the following:

	@HWI-ST330:304:H045HADXX:1:1101:1111:61397
	CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
	+
	@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################

Suppose we want to see how many reads in our file `Mov10_oe_1.subset.fq` are really bad, with 10 consecutive Ns (`NNNNNNNNNN`).

```bash
cd ~/ngs_course/unix_lesson/raw_fastq

grep NNNNNNNNNN Mov10_oe_1.subset.fq
```

We get back a lot of lines.  What if we want to see the whole fastq record for each of these reads? We can use the `-B` and `-A` arguments for grep to return the matched line plus one before (`-B 1`) and two lines after (`-A 2`). Since each record is four lines and the second line is the sequence, this should return the whole record.

```bash
grep -B 1 -A 2 NNNNNNNNNN Mov10_oe_1.subset.fq
```

```
@HWI-ST330:304:H045HADXX:1:1101:1111:61397
CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
+
@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
--
@HWI-ST330:304:H045HADXX:1:1101:1106:89824
CACAAATCGGCTCAGGAGGCTTGTAGAAAAGCTCAGCTTGACANNNNNNNNNNNNNNNNNGNGNACGAAACNNNNGNNNNNNNNNNNNNNNNNNNGTTGG
+
?@@DDDDDB1@?:E?;3A:1?9?E9?<?DGCDGBBDBF@;8DF#########################################################
```

## Exercise

1. Search for the sequence CTCAATGA in `Mov10_oe_1.subset.fq`.
In addition to finding the sequence, have your search also return
the name of the sequence.
2. Search for that sequence in all Mov10 replicate fastq files.



# Redirection

We're excited we have all these sequences that we care about that we
just got from the FASTQ files. That is a really important motif
that is going to help us answer our important question. But all those
sequences just went whizzing by with grep. How can we capture them?

We can do that with something called "redirection". The idea is that
we're redirecting the output from the terminal (all the stuff that went
whizzing by) to something else. In this case, we want to print it
to a file, so that we can look at it later.

The redirection command for putting something in a file is `>`

Let's try it out and put all the sequences that contain 'NNNNNNNNNN'
from all the files in to another file called `bad_reads.txt`.

```bash
grep -B 1 -A 2 NNNNNNNNNN Mov10_oe_1.subset.fq > bad_reads.txt
```

The prompt should sit there a little bit, and then it should look like nothing
happened. 

```bash
ls -l
```

You should have a new file called `bad_reads.txt`. Take a look at it and see if it has what you think it should.

If we use `>>`, it will append to rather than overwrite a file.  This can be useful for saving more than one search, for example.
    
```bash
grep -B 1 -A 2 NNNNNNNNNN Mov10_oe_2.subset.fq >> bad_reads.txt

ls -l
```

Since our `bad_reads.txt` file isn't a raw fastq file, we should move it to a different location within our directory. We decide to create a new folder called `other`, and move the `bad_reads.txt` to this `other` folder using the command `mv`. 

```bash
mkdir ../other/

mv bad_reads.txt ../other/
```

There's one more useful redirection command that we're going to show, and that's
called the pipe command, and it is `|`. It's probably not a key on
used very much on your keyboard. What `|` does is take the output that
scrolling by on the terminal and then can run it through another command.
When it was all whizzing by before, we wished we could just slow it down and
look at it, like we can with `less`. Well it turns out that we can! We pipe
the `grep` command to `less`.

```bash
grep -B 1 -A 2 NNNNNNNNNN Mov10_oe_1.subset.fq | less
```

Now we can use the arrows to scroll up and down and use `q` to get out.

We can also do something tricky and use the command `wc`. `wc` stands for
*word count*. It counts the number of lines or characters. So, we can use
it to count the number of lines we're getting back from our `grep` command.
And that will magically tell us how many sequences we're finding.

```bash
grep NNNNNNNNNN Mov10_oe_1.subset.fq | wc
```

This command tells us the number of lines, words and characters in the file. If we
just want the number of lines, we can use the `-l` flag for `lines`.

```bash
grep NNNNNNNNNN Mov10_oe_1.subset.fq | wc -l
```

Redirecting is not super intuitive, but it's powerful for stringing together these different commands, so you can do whatever you need to do.

The philosophy behind these commands is that none of them really do anything all that impressive. BUT when you start chaining them together, you can do some really powerful things 
efficiently. If you want to be proficient at using the shell, you must learn to become proficient with the pipe and redirection operators: `|`, `>`, `>>`.



# Practice with searching and redirection

Finally, let's use the new tools in our kit and a few new ones to examine our gene annotation file, **chr1-hg19_genes.gtf**, which we will be using later to find the genomic coordinates of all known exons on chromosome 1.

```bash
cd ../reference_data/
```

Let's explore our `chr1-hg19_genes.gtf` file a bit. What information does it contain?

```bash
less chr1-hg19_genes.gtf
```

```
chr1    unknown exon    14362   14829   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    14970   15038   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    15796   15947   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    16607   16765   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
chr1    unknown exon    16858   17055   .       -       .       gene_id "WASH7P"; gene_name "WASH7P"; transcript_id "NR_024540"; tss_id "TSS7245";
```

<div class="note">
The GTF file is a tab-delimited gene annotation file often used in NGS analyses. For more information on this file format, check out the [Ensembl site](http://useast.ensembl.org/info/website/upload/gff.html).
</div>

The columns in the **GTF file contain the genomic coordinates of gene features (exon, start_codon, stop_codon, CDS) and the gene_names, transcript_ids and protein_ids (p_id) associated with these features**. Note that sometimes an exon can be associated with multiple different transcripts or gene isoforms.

```bash
grep PLEKHN1 chr1-hg19_genes.gtf | head -n 5
```

This search returns two different transcripts of the same gene, NM_001160184 and NM_032129, that contain the same exon.

Now that we know what type of information is inside of our gtf file, let's explore our commands to answer a simple question about our data: **how many total exons are present on chromosome 1 using `chr1-hg19_genes.gtf`?**

To determine the number of total exons on chromosome 1, we are going to perform a series of steps:
	
1. Subset the dataset to only include the feature type and genomic location information
2. Extract only the genomic coordinates of exon features
3. Remove duplicate exons
4. Count the total number of exons

## Subsetting dataset

We will define an exon by it's genomic coordinates. Therefore, we only need the feature type and the genomic location (chr, start, stop, and strand) information to find the total number of exons. The columns corresponding to this information are 1, 3, 4, 5, and 7. 

'cut' is a program that will extract columns from files.  It is a very good command to know.  Let's first try out the 'cut' command on a small dataset (just the first 5 lines of chr1-hg19_genes.gtf) to make sure we have the command correct.

```bash
cut -f1,3,4,5,7 chr1-hg19_genes.gtf | head -n 5
```
   
`-f1,3,4,5,7` means to cut these fields (columns) from the dataset.

```
chr1	exon	14362	14829	-
chr1	exon	14970	15038	-
chr1	exon	15796	15947	-
chr1	exon	16607	16765	-
chr1	exon	16858	17055	-
```

The `cut` command assumes our data columns are separated by tabs (i.e. tab-delimited). The `chr1-hg19_genes.gtf` is a tab-delimited file, so the default `cut` command works for us. However, data can be separated by other types of delimiters. Another common delimiter is the comma, which separates data in comma-separated value (csv) files. If your data is not tab delimited, there is a `cut` command argument (-d) to specify the delimiter.

Our output looks good, so let's cut these columns from the whole dataset (not just the first 5 lines) and save it as a file, **`chr1-hg19genes_cut`**:

```bash
cut -f1,3,4,5,7 chr1-hg19_genes.gtf > chr1-hg19genes_cut
```

Check the cut file to make sure that it looks good using `less`. 

## Extracting genomic coordinates of exon features

We only want the exons (not CDS or start_codon features), so let's use `grep` to only keep the exon lines and save to file, **`chr1_exons`**:

```bash
grep exon chr1-hg19genes_cut > chr1_exons
```

## Removing duplicate exons

Now, we need to remove those exons that show up multiple times for different transcripts.    

We can use a new tool, `sort`, to remove exons that show up more than once.  We can use the `sort` command with the `-u` option to return only unique lines.

```bash
sort -u chr1_exons | head -n 5
```

## Counting the total number of exons

Now, to count how many unique exons are on chromosome 1, we need to pipe the output to `wc -l`:

```bash
sort -u chr1_exons | wc -l
```

## Exercise

How could have you have determined the number of total exons by combining all of the previous commands (starting with the original chr1-hg19_genes.gtf), into a single command (no intermediate files) using pipes?



# Where can I learn more about the shell?

- Software Carpentry tutorial: [The Unix shell](http://software-carpentry.org/v4/shell/index.html)
- The shell handout: [Command Reference](http://files.fosswire.com/2007/08/fwunixref.pdf)
- [explainshell.com](http://explainshell.com)
- http://tldp.org/HOWTO/Bash-Prog-Intro-HOWTO.html
- man bash
- Google: if you don't know how to do something, try Googling it. Other people
have probably had the same question.
- Learn by doing. There's no real other way to learn this than by trying it
out.  Write your next paper in nano (really emacs or vim), open pdfs from
the command line, automate something you don't really need to automate.



# Commands, options, and keystrokes covered in this lesson

```
grep
>       # output redirection
>>      # output redirection, append
|       # output redirection, pipe
wc
cut
sort
uniq
```



* * *



{% include copyright.md %}

Adapted from the lesson by Tracy Teal. Contributors: Paul Wilson, Milad Fatenejad, Sasha Wood, and Radhika Khetani for [Software Carpentry](http://software-carpentry.org/).
