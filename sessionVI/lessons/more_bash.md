---
title: "Bash_extras"
author: "Radhika Khetani", "Meeta Mistry"
date: 2019-10-29
duration: 30
---

## Overview

* [Creating shortcuts or `alias`](#alias)
* [Copying files using `scp` and `rsync`](#rsync) 
* [Symbolic Links or "sym links"](#symlink)

***

## Setting up some `alias`es <a name="alias"></a>

In your terminal, do the following:

```bash
$ cd

$ ls -l

$ ll
```

`ll` should have output the same thing as `ls -l`. Why does it work this way? This is because the HMS-RC folks have internally setup what is called an **alias**. 

A **shell alias is a shortcut to reference a command**. It can be used to avoid typing long commands. For common patterns it can reduce keystrokes and improve efficiency. A simple example is setting default options on commands to avoid having to type them each time a command is run.

For example suppose that because you are just starting out on the cluster, and you prefer to confirm deleting a file before using the `rm` command. Remember that the `rm` command supports this with the `-i` option. To avoid forgetting to use the `-i` option each time, an alias can be created so that each time `rm` is run it will use the `-i` option and prompt the user to confirm.


```bash
$ alias rm='rm -i'
```

However, this alias is only going to be available to you while that Terminal window is open. If you wanted to **use that alias all the time, what would you do?** 

You would add it to `~/.bashrc`! Let's open `~/.bashrc` and add a few commands to it. At the bottom of the file you should see a header titled "User specific aliases". Under that header go ahead and add the alias.

```bash
$ vim ~/.bashrc

# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions
```

Now, we can source the `.bashrc` file for the alias to take effect and we can try it out. You should see the question `
remove draft.txt?` and here you can answer `n` for No.

```bash
$ source ~/.bashrc

$ rm  ~/unix_lesson/other/draft.txt 
```

As we mentioned, aliases are super helpful for long commands that we are repeatedly having to tyoe out. A good example of this is the `srun` command for starting and interactive session. **First exit the interactive session and get on a login node, if you are not there already.**

```bash
$ alias o2i='srun --pty -p interactive -t 0-12:00 --mem 2G --reservation=HBC /bin/bash'
```

Now you can test it out!

```bash
$ o2i
```

Similar to what we did above, you can put this (or a similar) command in the `.bashrc` file so it is available when you log on next time.

> ### `.bashrc` versus `.bash_profile`
> `.bash_profile` is executed for login shells, while `.bashrc` is executed for interactive non-login shells. When you login (type username and password) to O2 the `.bash_profile` is executed. So if you want the alias available when you login, you will want to put it in your `.bash_profile`

## Copying files to and from the cluster <a name="rsync"></a>

So far we have used FileZilla to copy files over from O2, but there are other way to do so using the command line interface. When you obtain your data from the sequencing facility, it will likely be stored on some remote computer and they will give you login credentials which will allow you to access it. There are various commands that can be used to help you copy those files from the remote computer over to 1) your local computer, 2) O2, or 3) whatever cluster environment you plan to work on. We present a few options here.

### `scp`

Similar to the `cp` command to copy there is a command that allows you to **securely copy files between computers**. The command is called `scp` and allows files to be copied to, from, or between different hosts. It uses ssh for data transfer and provides the same authentication and same level of security as ssh. 

In the example below, the first argument is the **location on the remote server** and the second argument is the **destination on your local machine**. 

> *You can also do this in the opposite direction by swapping the arguments.*

```bash
$ scp username@transfer.rc.hms.harvard.edu:/path/to/file_on_O2 Path/to/directory/local_machine
```

Let's try copying over the README file from your `unix_lesson` folder:

```bash
$ scp rc_trainingXX@transfer.rc.hms.harvard.edu:~/unix_lesson/other/draft.txt  ~/Desktop/
```

### `rsync` 

`rsync` is used to copy or synchronize data between directories. It has many advantages over `cp`, `scp` etc. It works in a specific direction, i.e. from the first directory **to** the second directory, similar to `cp`.

**Salient Features of `rsync`**

* If the command (or transfer) is interrupted, you can start it again and *it will restart from where it was interrupted*.
* Once a folder has been synced between 2 locations, the next time you run `rsync` it will *only update and not copy everything over again*. 
* It runs a check to ensure that every file it is "syncing" over is the exact same in both locations. This check is run using a version of ["checksum"](https://en.wikipedia.org/wiki/Checksum) which ensures the data integrity during the data transfer process. 

> You can run the checksum function yourself when transferring large datasets without `rsync` using one of the following commands (or similar): `md5`, `md5sum`.


### Between directories on the same machine

```bash
#DO NOT RUN
$ rsync -av ~/large_dataset/. /n/groups/dir/groupdata/
```

### Between different machines

When copying over large datasets to or from a remote machine, `rsync` works similarly to `scp`.

```bash
#DO NOT RUN
$ rsync -av -e ssh testfile username@transfer.rc.hms.harvard.edu:~/large_files/
```

> Please do not use O2’s login nodes for transferring large datasets (like fastq files) between your computer and O2 with `rsync` or `scp`. Instead, use the transfer nodes `ssh eCommons@transfer.rc.hms.harvard.edu.


## Symbolic Links or "sym links" <a name="symlink"></a>

Symbolic links are like shortcuts you may create on your laptop. Let's check out an example of a folder with lots of symlinks.


```bash
ls -l /n/app/bcbio/tools/bin/
```

Now, let's create a shortcut in our home directory for that folder with the scripts for session IV homework.

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
