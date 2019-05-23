## Session VI: Variant Calling, Version Control, Intermediate Shell

### Description

Session VI will explore variant calling methods and tools for resequencing and identifying sequence variations. We will explore pre-call sequence processing, methods for variant calling and workflows, and functional interpretation of variant calls.

**Day1:** We will start the first day with an introduction to variant calling and a guest lecture from Dr. Brad Chapman. Following the lecture, we will discuss alignment tools and pre-call processing. After performing variant calling with FreeBayes we will filter our variant call data for quality. Next, we will annotate the high-confidence variants with known information and use SnpEff to obtain a functional annotation report, and we will finish the day finish with variant prioritization using Gemini.

**Day2:** The second day of this session will start with visualizing the variants identified on Day1 using IGV. Following the variant calling session, NGS workflow commonalities and suggestions for troubleshooting will be discussed (bring any questions regarding any of the NGS workflows!). After the coffee break, students will explore Git and Github for version control. This will be followed by several short topics including, incorporating version control into R, using Rmarkdowns to generate reports in RStudio, and regular expressions in shell along with a couple of bash commands that can help to improve efficiency.   

### Lessons
[Click here for the schedule with links to the lessons.](schedule)

### Learning Objectives
* Define appropriate BWA alignment parameters and clean-up steps
* Call variants with Freebayes
* Demonstrate how to access and add information to the Variant Call Format (VCF)
* Use vcftools to perform some simple filtering on the variants in the VCF file
* Incorporate annotation information to filter out important variants
* Explore variant information through the GEMINI framework
* Utilize git and Github for version control of scripts and R projects
* Create professional R analysis reports with Rmarkdown and knitr
* Demonstrate how to access genomic data, including genome builds, SRA, and other FTP sites
