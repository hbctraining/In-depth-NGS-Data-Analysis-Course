```bash
# Install Xcode CLT
$ xcode-select --install
```
```bash
# Install Homebrew
$ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

	# Uninstall Homebrew
	# /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall)"
```
```bash
# check that Homebrew is working properly
$ brew doctor
```
```bash
# Install Homebrew Cask
$ brew install caskroom/cask/brew-cask
```
```bash
# Use Cask to install OSXfuse
$ brew cask install osxfuse
```
```bash
# Use fuse to install sshfs
$ brew install homebrew/fuse/sshfs
```
```bash
# Create the folder with an intuitive name for the cluster to be mounted in
$ mkdir -p ~/Orchestra
```
```bash
# set up ssh keys
$ ssh-keygen -t rsa -b 4096 -f ~/.ssh/id_rsa -C "EMAIL"
$ ssh-add -K ~/.ssh/id_rsa
```
```bash
# copy to ~/.ssh/authorized_keys on Orchestra
$ cat ~/.ssh/id_rsa.pub | pbcopy
```
```bash
# create a shell script on your laptop

# mount orchestra with the appropriate login
sshfs USER@transfer.orchestra.med.harvard.edu:. ~/Orchestra -o volname="Orchestra" -o follow_symlinks

# for unmounting use `umount Orchestra` in `~/`
```
