---
layout: default
title: "Introduction to Versioning"
author: "Daniel van Strien, Radhika Khetani, Bob Freeman, Amir Karger"
output: html_document
---

#  Versioning your Data and Scripts

>> NOTE: Materials used in these lessons are derived/adapted from [Daniel van Strien's "An Introduction to Version Control Using GitHub Desktop," Programming Historian, (17 June 2016)](http://programminghistorian.org/lessons/getting-started-with-github-desktop). Licensing information available at the bottom of this page.

## What is Version Control?

We'll start by exploring how version control can be used to keep track of what one person did and when. Even if you aren't collaborating with other people, it is better than the scenario like this:

```
mydocument.txt
mydocumentversion2.txt
mydocumentwithrevision.txt
mydocumentfinal.txt
```
The system used for naming files may be more or less systematic. Adding dates makes it slightly easier to follow when changes were made:

```
mydocument2016-01-06.txt
mydocument2016-01-08.txt
```

Some word processors let us deal with this a little better, without creating a new file for every "save", such as Microsoft Word's "Track Changes" or Google Docs' [version history](https://support.google.com/docs/answer/190843?hl=en).

Version control systems start with a base version of the document and then save just the changes you made at each step of the way by taking a so-called "snapshot". A snapshot records information about when the it was taken, but also about what changes occurred between different snapshots. You decide when these snapshots are collected, and this allows you to ‘rewind’ your file to an older version. 

<img src="img/play-changes.png" width="600" align="center">

Once you think of changes as separate from the document itself, you can then think about "playing back" different sets of changes onto the base document and getting different versions of the document. For example, two users can make independent sets of changes based on the same document.

<img src="img/versions.png" width="400" align="center">

If there aren't conflicts, you can even play two sets of changes onto the same base document.

<img src="img/merged_example.png" width="400" align="center">

A Version Control system is a tool that keeps track of these changes for us and helps us version and merge our files.

## Why use Version Control?

The 2 main reasons to use version control are to:

* Manage data/code effectively 
* Collaborate efficiently

Though version control was originally designed for dealing with code (`.R`, `.pl`. `.py`) there are many benefits to using it to with ***text files*** too (`.txt`, `.csv`, `.tsv`). For substantial work such as articles, books, or dissertations, version control makes a lot of sense.

Though not all of these benefits will be covered in this lesson, version controlling your document allows you to:

* Track developments and changes in your documents (removing mental clutter!)
* Coordinate sets of changes as one unit
* Revert to previous versions of your document
* ‘Merge’ versions of a document and manage conflicts between versions
* Experiment with different versions of a document while maintaining the original version by creating branches

> Note: Different Version Control systems handle different non-text files differently. In most cases Word documents, graphics files, data objects from R or STATA, etc., can be included but most tools have limited capabilities for these.

Version control is particularly useful for facilitating collaboration. One of the original motivations behind version control systems was to allow different people to work on large projects together, and in the case of Git, to manage the Linux kernel source code. 

Benefits of collaborating with Version Control include:

* Flexibility and control
* Multiple people can simultaneously work on a document
* Easy conflict resolution 
* Easy to revert to an older version

## Why Not use Dropbox or Google Drive?

Dropbox, Google Drive and other services offer some form of version control in their systems. There are times when this may be sufficient for your needs. However there are a number of advantages to using a version control system like Git:

* **Language support**: Git supports both text and programming languages. As research moves to include more digital techniques and tools it becomes increasingly important to have a way of managing and sharing both the ‘traditional’ outputs (journal articles, books, etc.) but also these newer outputs (code, datasets etc.)
* **More control**: a proper version control systems gives you a much greater deal of control over how you manage changes in a document, including the ability to comment on every change making it easier to revert. This is especially true when sets of changes across documents need to be coordinated as one unit.
* **Useful history**: using version control systems like Git will allow you to produce a history of your document in which different stages of the documents can be navigated easily both by yourself and by others.

## What are Git and GitHub?

Though often used synonymously, Git and GitHub are two different things:

* Git is a particular implementation of version control originally designed by Linus Torvalds as a way of managing the Linux source code. Git can be used to refer both to a particular approach taken to version control and the software underlying it.

* GitHub is a company which hosts Git repositories (more on this below) and provides software for using Git. This includes ‘GitHub Desktop’ which will be covered in this tutorial. GitHub is currently the most popular host of open source projects by [number of projects and number of users](https://en.wikipedia.org/wiki/Comparison_of_source_code_hosting_facilities#Popularity). But other hosts exist, including [SourceForge](https://sourceforge.net/), [BitBucket](https://bitbucket.org/), and [Gitlab](https://about.gitlab.com/), to name a few.

Although GitHub’s focus is primarily on source code, other projects are increasingly making use of version control systems like GitHub to manage the work-flows of journal publishing, open textbooks, other humanities projects, and teaching materials.

Becoming familiar with GitHub will be useful not only for version controlling your own documents but will also make it easier to contribute and draw upon other projects which use GitHub. 

In this lesson the focus will be on gaining an understanding of the basic aims and principles of Version Control by uploading and version controlling a plain text document using Github. This lesson will not cover everything but will provide a starting point to using Version Control (Git/Github or [other](https://en.wikipedia.org/wiki/Comparison_of_version_control_software)).

***

* Materials used in these lessons are derived from Daniel van Strien's ["An Introduction to Version Control Using GitHub Desktop,"](http://programminghistorian.org/lessons/getting-started-with-github-desktop), Programming Historian, (17 June 2016). [The Programming Historian ISSN 2397-2068](http://programminghistorian.org/), is released under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*

* Materials are also derived from [Software Carpentry instructional material](https://swcarpentry.github.io/git-novice/). These materials are also licensed under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
***
