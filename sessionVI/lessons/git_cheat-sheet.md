---
layout: page
title: Version Control with Git
subtitle: Reference
---
## Git Cheatsheets for Quick Reference 
* A great [printable git cheatsheet](https://training.github.com/kit/downloads/github-git-cheat-sheet.pdf) is available in PDF from the GitHub training website.
* An [interactive one-page visualisation](http://ndpsoftware.com/git-cheatsheet.html) about the relationships between workspace, staging area, local repository, upstream repository, and the commands associated with each (with explanations).
* Both resources are also available in other languages e.g. Spanish, French, and more

## [Setting Up Git](02-setup.html)

*   Use `git config` to configure
    a user name, email address, editor, and other preferences once per machine.

## [Creating a Repository](03-create.html)

*   `git init` initializes a repository.

## [Tracking Changes](04-changes.html)

*   `git status` shows the status of a repository.
*   Files can be stored in a project's working directory (which users see),
    the staging area (where the next commit is being built up)
    and the local repository (where commits are permanently recorded).
*   `git add` puts files in the staging area.
*   `git commit` saves the staged content as a new commit in the local repository.
*   Always write a log message when committing changes.

## [Exploring History](05-history.html)

*   `git diff` displays differences between commits.
*   `git checkout` recovers old versions of files.

## [Ignoring Things](06-ignore.html)

*   The `.gitignore` file tells Git what files to ignore.

## [Remotes in GitHub](07-github.html)

*   A local Git repository can be connected to one or more remote repositories.
*   Use the HTTPS protocol to connect to remote repositories until you have learned how to set up SSH.
*   `git push` copies changes from a local repository to a remote repository.
*   `git pull` copies changes from a remote repository to a local repository.

## [Collaborating](08-collab.html)

*   `git clone` copies a remote repository to create a local repository
    with a remote called `origin` automatically set up.

## [Conflicts](09-conflict.html)

*   Conflicts occur when two or more people change the same file(s) at the same time.
*   The version control system does not allow people to blindly overwrite each other's changes.
    Instead, it highlights conflicts so that they can be resolved.

## [Open Science](10-open.html)

*   Open scientific work is more useful and more highly cited than closed.

## [Licensing](11-licensing.html)

*   People who incorporate GPL'd software into their own software must
    make their software also open under the GPL license; most other
    open licenses do not require this.
*   The Creative Commons family of licenses allow people to mix and
    match requirements and restrictions on attribution, creation of
    derivative works, further sharing, and commercialization.
*   People who are not lawyers should not try to write licenses from
    scratch.

## [Hosting](13-hosting.html)

*   Projects can be hosted on university servers,
    on personal domains,
    or on public forges.
*   Rules regarding intellectual property and storage of sensitive information apply
    no matter where code and data are hosted.

## Glossary

change set
:   A group of changes to one or more files that are or will be added
    to a single [commit](#commit) in a [version control](#version-control)
    [repository](#repository).

commit
:   To record the current state of a set of files (a [change set](#change-set))
    in a [version control](#version-control) [repository](#repository). As a noun, 
    the result of commiting, i.e. a recorded change set in a repository.
    If a commit contains changes to multiple files,
    all of the changes are recorded together.

conflict
:   A change made by one user of a [version control system](#version-control)
    that is incompatible with changes made by other users.
    Helping users [resolve](#resolve) conflicts
    is one of version control's major tasks.

HTTP
:   The Hypertext Transfer [Protocol](#protocol) used for sharing web pages and other data
    on the World Wide Web.

infective license
:   A license, such as the [GPL](http://opensource.org/licenses/GPL-3.0),
    that legally requires people who incorporate material under the
    infective license
    into their own work to also release under the same infective license
    (eg. under the GPL license).

merge
:   (a repository): To reconcile two sets of changes to a
    [repository](#repository).

protocol
:   A set of rules that define how one computer communicates with another.
    Common protocols on the Internet include [HTTP](#http) and [SSH](#ssh).

remote
:   (of a repository) A version control [repository](#repository) connected to another,
    in such way that both can be kept in sync exchanging [commits](#commit).

repository
:   A storage area where a [version control](#version-control) system
    stores the full history of [commits](#commit) of a project and information
    about who changed what, when.

resolve
:   To eliminate the [conflicts](#conflict) between two or more incompatible changes to a file or set of files
    being managed by a [version control](#version-control) system.

revision
:   A synonym for [commit](#commit).

SHA-1
:   [SHA-1 hashes](http://en.wikipedia.org/wiki/SHA-1) is what Git uses to compute identifiers, including for commits.
    To compute these, Git uses not only the actual change of a commit, but also its metadata (such as date, author,
    message), including the identifiers of all commits of preceding changes. This makes Git commit IDs virtually unique.
    I.e., the likelihood that two commits made independently, even of the same change, receive the same ID is exceedingly
    small.

SSH
:   The Secure Shell [protocol](#protocol) used for secure communication between computers.

timestamp
:   A record of when a particular event occurred.

version control
:   A tool for managing changes to a set of files.
    Each set of changes creates a new [commit](#commit) of the files;
    the version control system allows users to recover old commits reliably,
    and helps manage conflicting changes made by different users.
