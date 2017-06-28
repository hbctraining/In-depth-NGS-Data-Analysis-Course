To have orchestra accessible on your laptop/desktop as a folder, you need to use something called [sshfs](https://en.wikipedia.org/wiki/SSHFS) (ssh filesystem). This is a command that is not native to OSX or Windows and you need to go through several steps in order to get it. 

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
$ brew install caskroom/cask/brew-cask
```

Step 4. Install OSXfuse using Cask
```bash
$ brew cask install osxfuse
```

Step 5. Install sshfs from fuse
```bash
$ brew install homebrew/fuse/sshfs
```

Now, we have installed `sshfs`, the next step is to connect Orchestra (or a remote server) to our laptops. To make this process seamless, we have to set up ssh keys which can be used to connect to the server without having to type in a password everytime.

```bash
# set up ssh keys
$ ssh-keygen -t rsa -b 4096 -f ~/.ssh/id_rsa -C "EMAIL"
$ ssh-add -K ~/.ssh/id_rsa

# copy the contents of `id_rsa.pub` to ~/.ssh/authorized_keys on Orchestra
$ cat ~/.ssh/id_rsa.pub | pbcopy

# pbcopy puts the contents into the clipboard (in other words it is equivalent to copying with "ctrl + c") so you can just paste it as usual with "ctrl + v"
```
Use vim to open `~/.ssh/authorized_keys` and copy the contents from your computer to this file and save it. 

Now, let's set up for running `sshfs`, by creating a folder with an intuitive name for your home directory on the cluster to be mounted in.

```bash
$ mkdir ~/Orchestra
```

Finally, let's run the `sshfs` command to have Orchestra mount as a folder in the above space.
```bash
$ sshfs USER@transfer.orchestra.med.harvard.edu:. ~/Orchestra -o volname="Orchestra" -o follow_symlinks
```

Now we can browse through our home directory on Orchestra as though it was a folder on our laptop. 

> If you want to access your lab's directory in `/groups/` or your directory in `/n/scratch2`, you can create sym links to those in your home directory and you will be able to access those as well.

Once you are done with it, you can cancel the connection using `umount` and the name of the folder.

```bash
$ umount ~/Orchestra 
```

You can also create a shell script and anytime you want to mount Orchestra, you can run that script.

```bash
# mount orchestra with the appropriate login
sshfs USER@transfer.orchestra.med.harvard.edu:. ~/Orchestra -o volname="Orchestra" -o follow_symlinks

# for unmounting use `umount ~/Orchestra`
```
