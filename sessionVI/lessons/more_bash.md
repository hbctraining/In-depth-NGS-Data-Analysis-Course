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

### Setting up some `alias`es <a name="alias"></a>

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

alias o2='ssh <your_ecommons_ID>@orchestra.med.harvard.edu'
```

Now, open a new Terminal window, and try these out! You will still need to add your password, if you want to set up some "ssh keys" so that you don't have to enter your password you can find more information [within the O2 documentation](https://wiki.rc.hms.harvard.edu/display/O2/How+to+Generate+SSH+Keys).

### Copying files with `rsync` <a name="rsync"></a>

### Symbolic Links or "sym links" <a name="symlink"></a>

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
