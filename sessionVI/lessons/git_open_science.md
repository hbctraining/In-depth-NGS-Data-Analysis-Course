---
title: Version Control with Git
subtitle: Open Science
duration: 20
---

## Learning Objectives {.objectives}

* Explain how a version control system can be leveraged as an electronic lab notebook for computational work.
* Explain different options for hosting scientific work.
* Explain why adding licensing and citation information to a repository is important.
* Learn how to choose a proper license.

## Using Git for open science

> The opposite of "open" isn't "closed".
> The opposite of "open" is "broken".
>
> --- John Wilbanks

Free sharing of information might be the ideal in science,
but the reality is often more complicated.
Normal practice today looks something like this:

*   A scientist collects some data and stores it on a machine
    that is occasionally backed up by her department.
*   She then writes or modifies a few small programs
    (which also reside on her machine)
    to analyze that data.
*   Once she has some results,
    she writes them up and submits her paper.
    She might include her data&mdash;a growing number of journals require this&mdash;but
    she probably doesn't include her code.
*   Time passes.
*   The journal sends her reviews written anonymously by a handful of other people in her field.
    She revises her paper to satisfy them,
    during which time she might also modify the scripts she wrote earlier,
    and resubmits.
*   More time passes.
*   The paper is eventually published.
    It might include a link to an online copy of her data,
    but the paper itself will be behind a paywall:
    only people who have personal or institutional access
    will be able to read it.

For a growing number of scientists,
though,
the process looks like this:

*   The data that the scientist collects is stored in an open access repository
    like [figshare](http://figshare.com/) or
    [Zenodo](http://zenodo.org), possibly as soon as it's collected,
    and given its own 
    [Digital Object Identifier](https://en.wikipedia.org/wiki/Digital_object_identifier) (DOI).
    Or the data was already published and is stored in
    [Dryad](http://datadryad.org/).
*   The scientist creates a new repository on GitHub to hold her work.
*   As she does her analysis,
    she pushes changes to her scripts
    (and possibly some output files)
    to that repository.
    She also uses the repository for her paper;
    that repository is then the hub for collaboration with her colleagues.
*   When she's happy with the state of her paper,
    she posts a version to [arXiv](http://arxiv.org/)
    or some other preprint server
    to invite feedback from peers.
*   Based on that feedback,
    she may post several revisions
    before finally submitting her paper to a journal.
*   The published paper includes links to her preprint
    and to her code and data repositories,
    which  makes it much easier for other scientists
    to use her work as starting point for their own research.

This open model accelerates discovery:
the more open work is,
[the more widely it is cited and re-used](http://dx.doi.org/10.1371/journal.pone.0000308).
However,
people who want to work this way need to make some decisions
about what exactly "open" means and how to do it.

This is one of the (many) reasons we teach version control.
When used diligently,
it answers the "how" question
by acting as a shareable electronic lab notebook for computational work:

*   The conceptual stages of your work are documented, including who did
    what and when. Every step is stamped with an identifier (the commit ID)
    that is for most intents and purposes unique.
*   You can tie documentation of rationale, ideas, and other
    intellectual work directly to the changes that spring from them.
*   You can refer to what you used in your research to obtain your
    computational results in a way that is unique and recoverable.
*   With a distributed version control system such as Git, the version
    control repository is easy to archive for perpetuity, and contains
    the entire history.

> ### Making Code Citable {.callout}
>
> [This short guide](https://guides.github.com/activities/citable-code/) from GitHub
> explains how to create a Digital Object Identifier (DOI) for your code,
> your papers,
> or anything else hosted in a version control repository.

> ### How Reproducible Is My Work? {.challenge}
>
> Ask one of your labmates to reproduce a result you recently obtained
> using only what they can find in your papers or on the web.
> Try to do the same for one of their results,
> then try to do it for a result from a lab you work with.

## Hosting
The second big question for groups that want to open up their work
is where to host their code and data.
One option is for the lab, the department, or the university to provide a server,
manage accounts and backups,
and so on.
The main benefit of this is that it clarifies who owns what,
which is particularly important if any of the material is sensitive
(i.e.,
relates to experiments involving human subjects
or may be used in a patent application).
The main drawbacks are the cost of providing the service and its longevity:
a scientist who has spent ten years collecting data
would like to be sure that data will still be available ten years from now,
but that's well beyond the lifespan of most of the grants that fund academic infrastructure.

Another option is to purchase a domain
and pay an Internet service provider (ISP) to host it.
This gives the individual or group more control,
and sidesteps problems that can arise when moving from one institution to another,
but requires more time and effort to set up than either
the option above or the option below.

The third option is to use a public hosting service like
[GitHub](http://github.com), [BitBucket](http://bitbucket.org), or
[SourceForge](http://sourceforge.net).  Each of these services provides a web
interface that enables people to create, view, and edit their code repositories.
These services also provide communication and project management tools including
issue tracking, wiki pages, email notifications, and code reviews.  These
services benefit from economies of scale and network effects: it's easier to run
one large service well than to run many smaller services to the same standard.
It's also easier for people to collaborate.  Using a popular service can help
connect your project with communities already using the same service.

As an example,
Software Carpentry [is on GitHub](https://github.com/swcarpentry/)
where you can find the [source for this page](https://github.com/swcarpentry/git-novice/blob/gh-pages/13-hosting.md).  Anyone with a GitHub account can suggest changes to this text.

Using large, well-established services can also help you quickly take advantage
of powerful tools.  One such tool, continuous integration (CI), can
automatically run software builds and tests whenever code is committed or pull
requests are submitted.  Direct integration of CI with an online hosting service
means this information is present in any pull request, and helps maintain code
integrity and quality standards.  While CI is still available in self-hosted
situations, there is much less setup and maintenance involved with using an
online service.  Furthermore, such tools are often provided free of charge to
open source projects, and are also available for private repositories for a fee.

> ### Institutional Barriers {.callout}
>
> Sharing is the ideal for science,
> but many institutions place restrictions on sharing,
> for example to protect potentially patentable intellectual property.
> If you encounter such restrictions,
> it can be productive to inquire about the underlying motivations
> either to request an exception for a specific project or domain,
> or to push more broadly for institutional reform to support more open science.

> ### Can My Work Be Public? {.challenge}
>
> Find out whether you are allowed to host your work openly on a public forge.
> Can you do this unilaterally,
> or do you need permission from someone in your institution?
> If so, who?

## Licensing

When a repository with source code, a manuscript or other creative
works becomes public, it should include a file `LICENSE` or
`LICENSE.txt` in the base directory of the repository that clearly
states under which license the content is being made available. This
is because creative works are automatically eligible for intellectual
property (and thus copyright) protection. Reusing creative works
without a license is dangerous, because the copyright holders could
sue you for copyright infringement.

A license solves this problem by granting rights to others (the
licensees) that they would otherwise not have. What rights are being
granted under which conditions differs, often only slightly, from one
license to another. In practice, a few licenses are by far the most
popular, and [choosealicense.com](http://choosealicense.com/) will
help you find a common license that suits your needs.  Important
considerations include:

* Whether you want to address patent rights.
* Whether you require people distributing derivative works to also
  distribute their source code.
* Whether the content you are licensing is source code.
* Whether you want to license the code at all.

Choosing a licence that is in common use makes life easier for
contributors and users, because they are more likely to already be
familiar with the license and don't have to wade through a bunch of
jargon to decide if they're ok with it.  The [Open Source
Inititative](http://opensource.org/licenses) and [Free Software
Foundation](http://www.gnu.org/licenses/license-list.html) both
maintain lists of licenses which are good choices.

The following article provides an excellent overview of licensing and
licensing options from the perspective of scientists who also write
code:

> Morin, A., Urban, J., and Sliz, P. “[A Quick Guide to Software
> Licensing for the Scientist-Programmer](http://dx.doi.org/10.1371/journal.pcbi.1002598)” PLoS Computational Biology
> 8(7) (2012): e1002598.

At the end of the day what matters is that there is a clear statement
as to what the license is. Also, the license is best chosen from the
get-go, even if for a repository that is not public. Pushing off the
decision only makes it more complicated later, because each time a new
collaborator starts contributing, they, too, hold copyright and will
thus need to be asked for approval once a license is chosen.

> ### Can I Use Open License? {.challenge}
>
> Find out whether you are allowed to apply an open license to your software.
> Can you do this unilaterally,
> or do you need permission from someone in your institution?
> If so, who?

> ### What licenses have I already accepted? {.challenge}
>
> Many of the software tools we use on a daily basis (including in this workshop) are
> released as open-source software. Pick a project on GitHub from the list below, or
> one of your own choosing. Find its license (usually in a file called `LICENSE` or 
> `COPYING`) and talk about how it restricts your use of the software. Is it one of
> the licenses discussed in this session? How is it different?
> - [Git](https://github.com/git/git), the source-code management tool
> - [CPython](https://github.com/python/cpython), the standard implementation of the Python language
> - [Jupyter](https://github.com/jupyter), the project behind the web-based Python notebooks we'll be using
> - [EtherPad](https://github.com/ether/etherpad-lite), a real-time collaborative editor

## Making your work easy to cite

You may want to include a file called `CITATION` or `CITATION.txt`
that describes how to reference your project;
the [one for Software
Carpentry](https://github.com/swcarpentry/website/blob/gh-pages/CITATION)
states:

~~~
To reference Software Carpentry in publications, please cite both of the following:

Greg Wilson: "Software Carpentry: Getting Scientists to Write Better
Code by Making Them More Productive".  Computing in Science &
Engineering, Nov-Dec 2006.

Greg Wilson: "Software Carpentry: Lessons Learned". arXiv:1307.5448,
July 2013.

@article{wilson-software-carpentry-2006,
    author =  {Greg Wilson},
    title =   {Software Carpentry: Getting Scientists to Write Better Code by Making Them More Productive},
    journal = {Computing in Science \& Engineering},
    month =   {November--December},
    year =    {2006},
}

@online{wilson-software-carpentry-2013,
  author      = {Greg Wilson},
  title       = {Software Carpentry: Lessons Learned},
  version     = {1},
  date        = {2013-07-20},
  eprinttype  = {arxiv},
  eprint      = {1307.5448}
}
~~~
