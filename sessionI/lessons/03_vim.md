---
title: Introduction to Vim
author: "Mary Piper"
date: "Wednesday October 26, 2016"
---

Approximate time: 30 min

## Learning Objectives

* Learn basic operations using the Vim text editor

## Writing files

We've been able to do a lot of work with files that already exist, but what if we want to write our own files. Obviously, we're not going to type in a FASTA file, but you'll see as we go, there are a lot of reasons we'll want to write/create a file or edit an existing file.

To create or edit files we will need to use a **text editor**. When we say, "text editor," we really do mean "text": these editors can
only work with plain character data, not tables, images, or any other
media. The types of text editors available can generally be grouped into **graphical user interface (GUI) text editors** and **command-line editors**.

### GUI text editors

A GUI is an interface that has buttons and menus that you can click on to issue commands to the computer and you can move about the interface just by pointing and clicking. You might be familar with GUI text editors, such as [TextWrangler](http://www.barebones.com/products/textwrangler/), [Sublime](http://www.sublimetext.com/), and [Notepad++](http://notepad-plus-plus.org/), which allow you to write and edit plain text documents. These editors often have features to easily search text, extract text, and highlight syntax from multiple programming languages. They are great tools, but since they are 'point-and-click', we cannot efficiently use them from the command line remotely on a compute cluster.

### Command-line editors

When working remotely, we need a text editor that functions from the command line interface. Within these editors, since you cannot 'point-and-click', you must navigate the interface using the arrow keys and shortcuts. 

While there are simpler editors available for use (i.e. [nano](http://www.nano-editor.org/)), most computational scientists tend to favor editors that have greater functionality. Some popular editors include [Emacs](http://www.gnu.org/software/emacs/), [Vim](http://www.vim.org/), or a graphical editor such as [Gedit](http://projects.gnome.org/gedit/). These are editors which are generally available for use on high-performance compute clusters.

### Introduction to Vim 

To write and edit files, we're going to use a text editor called 'Vim'. Vim is a very powerful text editor, and it offers extensive text editing options. However, in this introduction we are going to focus on exploring some of the more basic functions. There is a lot of functionality that we are not going to cover during this session, but encourage you to go further as you become more comfortable using it. To help you remember some of the keyboard shortcuts that are introduced below and to allow you to explore additional functionality on your own, we have compiled a [cheatsheet](../../resources/VI_CommandReference.pdf).


### Vim Interface

You can create a document by calling a text editor and providing the name of the document you wish to create. Change directories to the `unix_lesson/other` folder and create a document using `vim` entitled `draft.txt`:

```bash
$ cd ~/ngs_course/unix_lesson/other
	
$ vim draft.txt
```

Notice the `"draft.txt" [New File]` typed at the bottom left-hand section of the screen. This tells you that you just created a new file in vim. 


### Vim Modes
Vim has **_two basic modes_** that will allow you to create documents and edit your text:   

- **_command mode (default mode):_** will allow you to save and quit the program (and execute other more advanced commands).  

- **_insert (or edit) mode:_** will allow you to write and edit text


Upon creation of a file, vim is automatically in command mode. Let's _change to insert mode_ by typing <kbd>i</kbd>. Notice the `--INSERT--` at the bottom left hand of the screen. Now type in a few lines of text:

![vim-insert-mode](../img/vim_insert.png)

After you have finished typing, press <kbd>esc</kbd> to enter command mode. Notice the `--INSERT--` disappeared from the bottom of the screen.

### Vim Saving and Quitting
To **write to file (save)**, type <kbd>:w</kbd>. You can see the commands you type in the bottom left-hand corner of the screen. 

![vim-save](../img/vim_save.png)

After you have saved the file, the total number of lines and characters in the file will print out at the bottom left-hand section of the screen.

![vim-postsave](../img/vim_postsave.png)

Alternatively, we can **write to file (save) and quit**. Let's do that by typing <kbd>:wq</kbd>. Now, you should have exited vim and returned back to your terminal window.

To edit your `draft.txt` document, open up the file again by calling vim and entering the file name: `vim draft.txt`. Change to insert mode and type a few more lines (you can move around the lines using the arrows on the keyboard). This time we decide to **quit without saving** by typing <kbd>:q!</kbd>
 
![vim-quit](../img/vim_quit.png)

### Vim Editing
Create the document "spider.txt" in vim. Enter the text as follows: 

![image](../img/vim_spider.png)

To make it easier to refer to distinct lines, we can add line numbers by typing <kbd>:set number</kbd>. **Save the document.** Later, if you choose to remove the line numbers you can type <kbd>:set nonumber</kbd>.

![image](../img/vim_spider_number.png)

While we cannot point and click to navigate the document, we can use the arrow keys to move around. Navigating with arrow keys can be very slow, so Vim has shortcuts (which are completely unituitive, but very useful as you get used to them over time). Check to see what mode you are currently in. While in command mode, try moving around the screen and familarizing yourself with some of these shortcuts:    

| key              | action                 |
| ---------------- | ---------------------- |
| <kbd>gg</kbd>     | to move to top of file |
| <kbd>G</kbd>     | to move to bottom of file     |
| <kbd>$</kbd>     | to move to end of line |
| <kbd>0</kbd>     | to move to beginning of line     |
| <kbd>w</kbd>     | to move to next word     |
| <kbd>b</kbd>     | to move to previous word     |


In addition to shortcuts for navigation, vim also offers editing shortcuts such as:

| key              | action                 |
| ---------------- | ---------------------- |
| <kbd>dw</kbd>     | to delete word |
| <kbd>dd</kbd>     | to delete line     |
| <kbd>u</kbd>     | to undo |
| <kbd>Ctrl + r</kbd>     | to redo     |
	
Practice some of the editing shortcuts, then quit the document without saving any changes.

*** 

**Exercise**

We have covered some basic commands in vim, but practice is key for getting comfortable with the program. Let's
practice what we just learned in a brief challenge.

1. Open spider.txt, and delete the word "water" from line #2.
2. Quit without saving.
3. Open `spider.txt` again, and delete: "Down came the rain." 
4. Save the file.
5. Undo your previous deletion.
6. Redo your previous deletion.
7. Delete the first and last words from each of the lines.
8. Save the file and see whether your results match your neighbors.

***

### Overview of vim commands

**Vim modes:**


| key              | action                 |
| ---------------- | ---------------------- |
| <kbd>i</kbd>     | insert mode - to write and edit text |
| <kbd>esc</kbd>     | command mode - to issue commands / shortcuts  |


**Saving and quiting:**


| key              | action                 |
| ---------------- | ---------------------- |
| <kbd>:w</kbd>     | to write to file (save) |
| <kbd>:wq</kbd>     | to write to file and quit     |
| <kbd>:q!</kbd>     | to quit without saving |


**Shortcuts for navigation:**


| key              | action                 |
| ---------------- | ---------------------- |
| <kbd>gg</kbd>     | to move to top of file |
| <kbd>G</kbd>     | to move to bottom of file     |
| <kbd>$</kbd>     | to move to end of line |
| <kbd>0</kbd>     | to move to beginning of line     |
| <kbd>w</kbd>     | to move to next word     |
| <kbd>b</kbd>     | to move to previous word     |

**Shortcuts for editing:**

| key              | action                 |
| ---------------- | ---------------------- |
| <kbd>dw</kbd>     | to delete word |
| <kbd>dd</kbd>     | to delete line     |
| <kbd>u</kbd>     | to undo |
| <kbd>Ctrl + r</kbd>     | to redo     |
| <kbd>:set number</kbd>     | to number lines |
| <kbd>:set nonumber</kbd>     | to remove line numbers    |
	
---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
