---
title: "Permissions and Environment Variables"
author: "Christina Koch, Radhika Khetani"
date: "Sunday, May 28, 2017"
---

Approximate time: 60 minutes

## Learning Objectives
 
* How to grant or restrict access to files on a multi-user UNIX system
* What is an "Environment Variable" in a shell.
* What is $PATH, and why I should care.

## Permissions

UNIX controls who can read, modify, and run files using *permissions*.

Let's start with how users are identified in a shared, multi-user system.
We all have a unique username, e.g. rsk27 and a userid 124292.

Find out yours:

```bash
$ id <username>
```

Users of a multi-user UNIX system can belong to any number of groups, each of which has a unique group name, and a numeric group ID.

The list of who is in what group is usually stored in the system file `/etc/group`.

Let's see what groups we all belong to:

```bash
$ groups
```

Depending on our affiliation, we all belong to at least a couple of groups. I belong to 5 groups,
* rsk27
* bcbio
* hbctraining
* Domain_Users
* genomebrowser-uploads 
* shen

As you can imagine, on a shared system it is important to protect each user's data. To start, every file and directory on a Unix computer belongs to one owner and one group. Along with each file's content, the operating system stores the numeric IDs of the user and group that own it, which is the "metadata" for a given file.

The user-and-group model means that for each file every user on the system falls into one of three categories:

* the owner of the file,
* someone in the file's group,
* and everyone else.

For each of these three categories, the computer keeps track of whether people in that category can read the file, write to the file, or execute the file (i.e., run it if it is a program).

For example, if a file had the following set of permissions:

<table class="table table-striped">
<tr><td></td><th>user</th><th>group</th><th>all</th></tr>
<tr><th>read</th><td>yes</td><td>yes</td><td>no</td></tr>
<tr><th>write</th><td>yes</td><td>no</td><td>no</td></tr>
<tr><th>execute</th><td>no</td><td>no</td><td>no</td></tr>
</table>

it would mean that:

*   the file's owner can read and write it, but not run it;
*   other people in the file's group can read it, but not modify it or run it; and
*   everybody else can do nothing with it at all.

Let's look at this model in action.

If we say,

```bash
$ ls -l /bin/ls
```

`-rwxr-xr-x. 1 root root 109208 Oct 15  2014 /bin/ls`. 
 
So, `ls` is an executable file that belong to user root and group root, and only they can modify (write) it.

> ### Necessary But Not Sufficient
>
> The fact that something is marked as executable doesn't actually mean it contains or is a program of some kind. We could easily mark the `~/ngs_course/unix_lesson/raw_fastq/Irrel_kd_1.subset.fq` file as executable using the commands that are introduced below. Depending on the operating system we're using, trying to "run" it will fail (because it doesn't contain instructions the computer recognizes).

Now let's run the following command:

```bash	
$ ls -l ~/ngs_course/unix_lesson
```
```
drwxrwsr-x 2 rsk27 rsk27  78 Oct  6 10:29 genomics_data
drwxrwsr-x 2 rsk27 rsk27 228 Oct  6 10:28 raw_fastq
-rw-rw-r-- 1 rsk27 rsk27 377 Oct  6 10:28 README.txt
drwxrwsr-x 2 rsk27 rsk27 238 Oct  6 10:28 reference_data
```

The `-l` flag tells `ls` to give us a long-form listing. It's a lot of information, so let's go through the columns from right to left.

On the right side we have the file names and to the left of them are the times and dates these files were last modified. Backup systems and other tools use this information in a variety of ways, but you can use it to tell when you (or anyone else with permission) last changed a file.

To the left of the modification time is the file's size in bytes and the names of the user and group that owns it (in this case, `rsk27` and `rsk27` respectively). We'll skip over the second column for now (the one showing `1` for each file),  because it's the first column that we care about most. This shows the file's permissions, i.e., who can read, write, or execute it.

Let's have a closer look at one of those permission strings for README.txt:
	
```
-rw-rw-r--
```

The first character tells us what type of thing this is: '-' means it's a regular file, while 'd' means it's a directory, and other characters mean more esoteric things.

The next three characters tell us what permissions the file's owner has. Here, the owner can read and write the file: `rw-`.

The middle triplet shows us the group's permissions. If the permission is turned off, we see a dash, so `rw-` means "read and write, but not execute". (In this case the group and the owner are the same so it makes sense that this is the same for both categories.)

The final triplet shows us what everyone who isn't the file's owner, or in the file's group, can do. In this case, it's `r--` again, so everyone on the system can look at the file's contents.

To change permissions, we use the `chmod` command (whose name stands for "change mode"). Let's make our README.txt file **inaccessible** to all users other than you and your group, currently they are able to read it:

```bash
$ ls -l ~/ngs_course/unix_lesson/README.txt

-rw-rw-r-- 1 rsk27 rsk27 377 Oct  6 10:28 /home/rsk27/ngs_course/unix_lesson/README.txt
```

```bash
$ chmod o-rw ~/ngs_course/unix_lesson/README.txt         # the "-" after o denotes removing that permission
```

```bash
$ ls -l ~/ngs_course/unix_lesson/README.txt

-rw-rw---- 1 rsk27 rsk27 377 Oct  6 10:28 /home/rsk27/ngs_course/unix_lesson/README.txt
```

The `o` signals that we're changing the privileges of "others", and the `-` indicates that we are removing read and write permissions.

Let's change it back to allow it to be readable by others, i.e. add read permission:
	
```bash
$ chmod o+r ~/ngs_course/unix_lesson/README.txt         # the "+" after o denotes adding/giving that permission
```

```bash
$ ls -l ~/ngs_course/unix_lesson/README.txt

-rw-rw-r-- 1 rsk27 rsk27 377 Oct  6 10:28 /home/rsk27/ngs_course/unix_lesson/README.txt
```

If we wanted to make this an executable file for ourselves (the file's owners) we would say `chmod u+rwx`, where the `u` signals that we are changing permission for the file's owner. To change permissions for a whole group, you'd use the letter `g`, e.g. `chmod g-w`. 

Before we go any further,
let's run `ls -l` on the `~/ngs_course/unix_lesson` directory to get a long-form listing:

```bash
$ ls -l

drwxrwsr-x 2 rsk27 rsk27  78 Oct  6 10:29 genomics_data
drwxrwsr-x 2 rsk27 rsk27 228 Oct  6 10:28 raw_fastq
-rw-rw-r-- 1 rsk27 rsk27 377 Oct  6 10:28 README.txt
drwxrwsr-x 2 rsk27 rsk27 238 Oct  6 10:28 reference_data
```

Look at the permissions for directories (`drwxrwsr-x`): the 'x' indicates that "execute" is turned on. What does that mean? A directory isn't a program or an executable file, we can't "run" it.

Well, 'x' means something different for directories. It gives someone the right to *traverse* the directory, but not to look at its contents. The distinction is subtle, so let's have a look at an example.

Dr. Vlad Smith's home directory has three subdirectories called `venus`, `mars`, and `pluto`:

![execute](../img/permission-directory.png "Execute Permission for Directories")

Each of these has a subdirectory in turn called `notes`, and those sub-subdirectories contain various files.
If a user's permissions on `venus` are 'r-x', then if she tries to see the contents of `venus` and `venus/notes` using `ls`, the computer lets her see both.
If her permissions on `mars` are just 'r--', then she is allowed to read the contents of both `mars` and `mars/notes`.
But if her permissions on `pluto` are only '--x', she cannot see what's in the `pluto` directory: `ls pluto` will tell her she doesn't have permission to view its contents.
If she tries to look in `pluto/notes`, though, the computer will let her do that.
She's allowed to go through `pluto`, but not to look at what's there. She will be able to do this, only if she knows that there is a file called `notes` in the directory, since she cannot list what is in there.

This trick gives people a way to make some of their directories visible to the world as a whole without opening up everything else.

****
### Exercise

If `ls -l myfile.php` returns the following details:

```bash
-rwxr-xr-- 1 caro zoo  2312  2014-10-25 18:30 myfile.php
```
 
Which of the following statements is true?
 
1. caro (the owner) can read, write, and execute myfile.php
2. caro (the owner) cannot write to myfile.php
3. members of caro (a group) can read, write, and execute myfile.php
4. members of zoo (a group) cannot execute myfile.php
****

## Environment Variables

Every time a shell session spawns, a process takes place to gather and compile information to determine its behavior and access to resources. One way that the shell keeps track of all of these settings and details is through an area it maintains called the **environment**.

The environment is built by the shell every time that it starts a session. The environment is defined by **environment variables** as they define the system properties of the environment. Two commonly encountered variables are `HOME` and `PATH`.

* `HOME` defines the home directory for a user.
* `PATH` defines a list of directories to search through when looking for a command to execute.

In the context of the shell the Environment variables are usually all upper case.

First, let's see our list of environment variables:

```bash
$ env
```

Let's see what is stored in these variables:

```bash
$ echo $HOME

/home/rsk27
```

Variables, in most systems, are called/denoted with a "$" before the variable name

```bash
$ echo $PATH

/opt/lsf/7.0/linux2.6-glibc2.3-x86_64/bin:/groups/bcbio/bcbio/anaconda/bin:/opt/bcbio/local/bin:/opt/lsf/7.0/linux2.6-glibc2.3-x86_64/etc:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin
```
I have a lot of full/absolute paths in my $PATH variable, which are separated from each other by a "**:**"; here is the list in a more readable format:

* /opt/lsf/7.0/linux2.6-glibc2.3-x86_64/bin
* /groups/bcbio/bcbio/anaconda/bin
* /opt/bcbio/local/bin
* /opt/lsf/7.0/linux2.6-glibc2.3-x86_64/etc
* /usr/local/bin
* /bin
* /usr/bin
* /usr/local/sbin
* /usr/sbin
* /sbin

These are the directories that the shell will look in (in the same order as they are listed) for an executable file that you type on the command prompt. 

When someone says a command or an executable file is "*in you path*", they mean that the parent directory for that command/file is contained in the list in the PATH variable. 

For any command you execute on the command prompt, you can find out where they are located using the which command.

Try it on a few of the basic commands we have learned so far:
```bash
$ which ls
$ which <your favorite command>
$ which <your favorite command>
```

> #### Modifying Environment Variables
>
> If you are interested in adding a new entry to the path variable, the command to use is `export`. This command is usually executed as follows: 
`export PATH=$PATH:~/opt/bin`, which tells the shell to add the ~/opt/bin directory to the end of the preexisiting list within $PATH. Alternatively, if you use `export PATH=~/opt/bin:$PATH`, the same directory will be added to the beginning of the list. The order determines where the shell will look first.

#### Closer look at the inner workings of the shell, in the context of $PATH
 
The $PATH variable is reset to a set of defaults (/bin:/usr/bin and so on), each time you start a new shell terminal. To make sure that a command/program you need is always at your fingertips, you have to put it in one of 2 special shell scripts that are always run when you start a new terminal. These are hidden files in your home directory called `.bashrc` and `.bash_profile`. You can create them if they don't exist, and shell will use them.

Check what hidden files exist in our home directory:
```bash
$ ls -al ~/
```

Open the .bashrc file and at the end of the file add the export command that adds a specific location to the list in $PATH. This way when you start a new shell, that location will always be in your path. 

The location we want to add to the beginning of the list is `/opt/bcbio/centos/bin`, we need this for later in the course.

```bash
$ vim ~/.bashrc

# at the end of the file type in the following - "export PATH=/opt/bcbio/centos/bin:$PATH"
# Don't forget the ":" between!
```

**In closing, permissions and environment variables, especially $PATH, are very useful and important concepts to understand in the context of UNIX and HPC.**

### [[Extra, extra!!]] Modifying command prompt content

To change the command prompt you would modify the `PS1` environment variable. Let's take a look at what the current contents are.

```bash
echo $PS1
```

To change the command prompt to contain the full path of the directory you are in, open the `~/.bashrc` file using "vim" and add the following line to the end of the file.

```bash
PS1="\u@\h:\w\$ "
```
The `\u`, `\h`, `\w` and `\W` are "bash prompt special characters." You can find more information about what those characters mean [at this link](https://linuxconfig.org/bash-prompt-basics).

After you add it, run the `source` command on the file to make the modification go into effect.

```bash
$ source ~/.bashrc
```

These steps will modify the command prompt when you start an interactive session on Orchestra to display the full path. If you want to make the same change to the command prompt when you log in, you will have to add the `PS1` reassignment line to `~/.bash_profile`.

---

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson was derived from work that is Copyright © Software Carpentry (http://software-carpentry.org/). 
All Software Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
