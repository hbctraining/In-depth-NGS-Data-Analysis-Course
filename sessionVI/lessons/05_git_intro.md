---
title: Version Control with Git
subtitle: Introduction to version control systems
duration: 20
---

Learning Objectives:

- Understand the benefits of an automated version control system.



# Introduction to version control systems

Version control systems are designed to manage changes to files or sets of files. Changes made to files that are committed to version control are never truly lost. Since all old versions of files are saved, it’s always possible to go back in time to see exactly what was written on a particular day, or what version of a program was used to generate a particular set of results, and, if needed, we can revert to a previous version, much like the “undo” feature in an editor.

Version control can also be used when a group of people can edit the same files, and it can record who made what changes when. When several people collaborate in the same project, it’s possible to accidentally overlook or overwrite someone’s changes: the version control system automatically notifies users whenever there’s a conflict between one person’s work and another’s.

Teams are not the only ones to benefit from version control: lone researchers can benefit immensely. Keeping a record of what was changed, when, and why is extremely useful for all researchers if they ever need to come back to the project later on (e.g., a year later, when memory has faded).

Version control is the lab notebook of the digital world: it’s what professionals use to keep track of what they’ve done and to collaborate with other people. Every large software development project relies on it, and most programmers use it for their small jobs as well. And it isn’t just for software: books, papers, small data sets, and anything that changes over time or needs to be shared can and should be stored in a version control system.


## Exploring version control

We'll start by exploring how version control can be used
to keep track of what one person did and when.
Even if you aren't collaborating with other people,
automated version control is much better than this situation:

[![Piled Higher and Deeper by Jorge Cham, http://www.phdcomics.com](../img/phd101212s.gif)](http://www.phdcomics.com)

"Piled Higher and Deeper" by Jorge Cham, http://www.phdcomics.com

We've all been in this situation before: it seems ridiculous to have multiple nearly-identical versions of the same document. Some word processors let us deal with this a little better, such as Microsoft Word's "Track Changes" or Google Docs' version history.

Version control systems start with a base version of the document and then save just the changes you made at each step of the way. You can think of it as a tape: if you rewind the tape and start at the base document, then you can play back each change and end up with your latest version.

![Changes are saved sequentially](https://cdn.rawgit.com/hbc/NGS_Data_Analysis_Course/master/sessionVI/img/play-changes.svg)

Once you think of changes as separate from the document itself, you can then think about "playing back" different sets of changes onto the base document and getting different versions of the document. For example, two users can make independent sets of changes based on the same document.

![Different versions can be saved](https://cdn.rawgit.com/hbc/NGS_Data_Analysis_Course/master/sessionVI/img/versions.svg)

If there aren't conflicts, you can even try to play two sets of changes onto the same base document.

![Multiple versions can be merged](https://cdn.rawgit.com/hbc/NGS_Data_Analysis_Course/master/sessionVI/img/merge.svg)

A version control system is a tool that keeps track of these changes for us and
helps us version and merge our files. It allows you to
decide which changes make up the next version, called a
[commit](../reference.html#commit), and keeps useful metadata about them. The
complete history of commits for a particular project and their metadata make up
a [repository](../reference.html#repository). Repositories can be kept in sync
across different computers facilitating collaboration among different people.

> ### The long history of version control systems
>
> Automated version control systems are nothing new. 
> Tools like RCS, CVS, or Subversion have been around since the early 1980s and are used by many large companies.
> However, many of these are now becoming considered as legacy systems due to various limitations in their capabilities.
> In particular, the more modern systems, such as Git and [Mercurial](http://swcarpentry.github.io/hg-novice/) 
> are *distributed*, meaning that they do not need a centralized server to host the repository.
> These modern systems also include powerful merging tools that make it possible for multiple authors to work within 
> the same files concurrently.



# Git version control system

Git is a widely-used, free and open-source version control system. The basic command structure for saving versions of files in Git is relatively easy to learn, with features allowing you to examine the new versions of your files prior to saving. In addition, Git has more powerful features to aid in group project work. [Git user documentation](https://git-scm.com/book) is accessible and thorough, and simple guides for quick command look-ups are widely available, such as [Git - the simple guide](http://rogerdudler.github.io/git-guide/) and [Bitbucket's Basic Git command documentation](https://confluence.atlassian.com/bitbucketserver/basic-git-commands-776639767.html). So while the intricacies of Git may be a bit overwhelming at first, there is a lot of documentation to help you along the way. 

While Git is fast, with most actions being performed on your local computer, we often use Git with a web-based Git repository hosting service such [Github](https://github.com) or [Bitbucket](https://bitbucket.org) to store our file repositories and to share them with others. By the end of these lessons, you should be able to perform routine file version back-ups and collaborate simply with collaborators using Git and Github.
