---
title: "Easy access to Orchestra using sshfs"
author: "Michael Steinbaugh, Lorena Pantano, Radhika Khetani"
date: "Wednesday, July 28th, 2017"
---

Approximate time: 30 minutes

## Learning Objectives

* Access a remote server as a folder on a local computer

## For Windows

> *Note that these are untested*.

#### OPTION 1 (Try first)

[http://www.sftpnetdrive.com/download-thanks](http://www.sftpnetdrive.com/download-thanks)

#### OPTION 2 (Try only if Option 1 fails!)

Download the .exe file for Dokan from [https://github.com/dokan-dev/dokany/releases/tag/v0.7.4](https://github.com/dokan-dev/dokany/releases/tag/v0.7.4) and run it.

Next, download the .zip file for win-sshfs [https://github.com/Foreveryone-cz/win-sshfs/releases/tag/1.5.12.8](https://github.com/Foreveryone-cz/win-sshfs/releases/tag/1.5.12.8), unzip it and run the .exe file within.


## For OSX

To have orchestra accessible on your laptop/desktop as a folder, you need to use something called [`sshfs`](https://en.wikipedia.org/wiki/SSHFS) (ssh filesystem). This is a command that is not native to OSXand you need to go through several steps in order to get it. Once you have `sshfs`, then you need to set up ssh keys to connect Orchestra to your laptop without having to type in a password. 

### 1. Installing sshfs on OSX

Below are 2 ways to get sshfs on OSX or macs, and I am listing both since one might work better on some versions of OSX than others.

#### OPTION 1 (Try first)

Download OSXfuse from [https://github.com/osxfuse/osxfuse/releases](https://github.com/osxfuse/osxfuse/releases/download/osxfuse-3.6.0/osxfuse-3.6.0.dmg), and install it.

Download sshfs from [https://github.com/osxfuse/sshfs/releases](https://github.com/osxfuse/sshfs/releases/download/osxfuse-sshfs-2.5.0/sshfs-2.5.0.pkg), and install it.

#### OPTION 2 (Try only if Option 1 fails!)

Step 1. Install [Xcode](https://developer.apple.com/xcode/)
```bash
$ xcode-select --install
```

Step 2. Install Homebrew using ruby (from Xcode)
```bash
$ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

	# Uninstall Homebrew
	# /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall)"
```

Step 2.1. Check to make sure that Homebrew is working properly
```bash
$ brew doctor
```

Step 3. Install Cask from Homebrew's caskroom
```bash
$ brew tap caskroom/cask
```

Step 4. Install OSXfuse using Cask
```bash
$ brew cask install osxfuse
```

Step 5. Install sshfs from fuse
```bash
$ brew install homebrew/fuse/sshfs
```

### 2. Set up "ssh keys"

Now, we have installed `sshfs`, the next step is to connect Orchestra (or a remote server) to our laptops. To make this process seamless, we will first set up ssh keys which can be used to connect to the server without having to type in a password everytime.

```bash
# set up ssh keys
$ ssh-keygen -t rsa -b 4096 -f ~/.ssh/id_rsa -C "ecommonsID"
$ ssh-add -K ~/.ssh/id_rsa
```

Arguments for `ssh-keygen`:
* `-t` = Specifies the type of key to create. The possible values are "rsa1" for protocol version 1 and "rsa" or "dsa" for protocol version 2. *We want rsa.*
* `-b` = Specifies the number of bits in the key to create. For RSA keys, the minimum size is 768 bits and the default is 2048 bits. *We want 4096*
* `-f` = name of output "keyfile"
* `-C` = Provides a new comment

Arguments for `ssh-add`:
* `-K` = Store passphrases in your keychain

```bash
# copy the contents of `id_rsa.pub` to ~/.ssh/authorized_keys on Orchestra
$ cat ~/.ssh/id_rsa.pub | pbcopy
```

> `pbcopy` puts the output of `cat` into the clipboard (in other words, it is equivalent to copying with <kbd>ctrl + c</kbd>) so you can just paste it as usual with <kbd>ctrl + v</kbd>.


Log into Orchestra and use vim to open `~/.ssh/authorized_keys` and copy the contents from your computer to this file and save it. 


### 3. Mount Orchestra using sshfs

Now, let's set up for running `sshfs` on our laptops (local machines), by creating a folder with an intuitive name for your home directory on the cluster to be mounted in.

```bash
$ mkdir ~/Orchestra
```

Finally, let's run the `sshfs` command to have Orchestra mount as a folder in the above space.
```bash
$ sshfs ecommonsID@transfer.orchestra.med.harvard.edu:. ~/Orchestra -o volname="Orchestra" -o follow_symlinks
```

Now we can browse through our home directory on Orchestra as though it was a folder on our laptop. 

> If you want to access your lab's directory in `/groups/` or your directory in `/n/scratch2`, you can create sym links to those in your home directory and you will be able to access those as well.

Once you are done with it, you can cancel the connection using `umount` and the name of the folder.

```bash
$ umount ~/Orchestra 
```

#### Create an "alias" for mounting and logging into orchestra

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

alias orchestra='ssh ecommonsID@orchestra.med.harvard.edu'

alias orch_mount='sshfs ecommonsID@transfer.orchestra.med.harvard.edu:. ~/Orchestra -o volname="Orchestra" -o follow_symlinks'
```

Now, open a new Terminal window, or source the file you just modified, and try these out!

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
