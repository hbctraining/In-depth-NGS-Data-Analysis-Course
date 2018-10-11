---
title: "Bash_extras"
author: "Radhika Khetani"
date: 2017-07-06
duration: 30
---

## Overview

* [Creating shortcuts or `alias`](#alias)
* [Copying files with `rsync`](#rsync) 
* [Symbolic Links or "sym links"](#symlink)

***

## Setting up some `alias`es <a name="alias"></a>

On your local machine do the following:

```bash
$ cd

$ ls -l

$ ll
```

`ll` should not work for you, but it works on my computer, why? It's because I have set up an alias for my bash environment, using the `alias` command, such that it knows that I want to actually do `ls -l` when I say `ll`. Let's set it up for your environment.

```bash
$ alias ll='ls -l'

$ ll
```

This alias is only going to be available to you while that Terminal window is open. If you wanted to use that alias all the time, what would you do? 

You would add it to `~/.bashrc` or `~/.bash_profile`!

Let's open either `~/.bash_profile` or `~/.bashrc` files on your laptop (*not on orchestra*), and add a few commands to it.

```bash
alias ll='ls -l'

alias o2='ssh <your_ecommons_ID>@o2.hms.harvard.edu'
```

Now, open a new Terminal window, and try these out! You will still need to add your password, if you want to set up some "ssh keys" so that you don't have to enter your password you can find more information [within the O2 documentation](https://wiki.rc.hms.harvard.edu/display/O2/How+to+Generate+SSH+Keys).

You can now create an alias to run an interactive session on O2!

```bash
alias interactive=`srun --pty -p interactive bash`

# and/OR

alias interactive6='srun --pty -p interactive -n 6 --mem 8G bash'
```

You can not try one of them out from the login node.

Similar to what we did above, you can put this (or a similar) command in the `.bashrc` or `.bashprofile` files so it is available when you log on next time.


## Copying files with `rsync` <a name="rsync"></a>

`rsync` is used to copy or synchronize data between directories. It has many advantages over `cp`, `scp` etc. It works in a specific direction, i.e. from the first diretory **to** the second directory, similar to `cp`.


### Between directories on the same machine

```bash
#DO NOT RUN
$ rsync -av ~/large_dataset/. /n/groups/dir/groupdata/
```

### Between different machines

When copying over large datasets to or from a remote machine, `rsync` works similarly to `scp`.

```bash
#DO NOT RUN
$ rsync -av -e ssh testfile <your_ecommons_ID>@transfer.o2.hms.harvard.edu:~/large_files/
```

> Please do not use Orchestra’s login servers for heavy I/O jobs like rsync or sftp. When transfering large files to and from O2, use their transfer server `transfer.o2.hms.harvard.edu`.

**Salient Features of `rsync`**

* If the command (or transfer) is interrupted, you can start it again and *it will restart from where it was interrupted*.
* Once a folder has been synced between 2 locations, the next time you run `rsync` it will *only update and not copy everything over again*. 
* It runs a check to ensure that every file it is "syncing" over is the exact same in both locations. This check is run using a version of ["checksum"](https://en.wikipedia.org/wiki/Checksum). 

> You can run the checksum function yourself when transferring large datasets without `rsync` using one of the following commands (or similar): `md5`, `md5sum`.


## Symbolic Links or "sym links" <a name="symlink"></a>

Symbolic links are like shortcuts you may create on mac. Make a file or directory and then a symlink to it called myshortcut:

```bash
ls -l /n/app/bcbio/tools/bin/
```


```bash
$ cd

$ ln -s /n/groups/hbctraining/ngs-data-analysis-longcourse/sessionIV_hmwk/ rnaseq_homework_NGScourse

$ ls -l
```

We recommend that you create something like this for your raw data so it does not accidentally get corrupted or overwritten. 

> Note: a “hard” link (just `ln` without the `-s` option) is very different. Always use “ln -s” unless you really know what you’re doing!

## Additional topics

If you are interested in learning more about regular expressions (regex) and the tools `awk` and `sed1`, you can find more information in the ["extra_bash_tools"](extra_bash_tools.md) lesson.


***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
