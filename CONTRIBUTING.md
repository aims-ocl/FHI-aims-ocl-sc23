# Contributing to FHI-aims #

## Getting Started with Coding in FHI-aims ##

We highly recommend that all users and developers join the FHI-aims Slack
channel at <http://fhi-aims.slack.com/> .  Instructions for joining the Slack
channel may be found in the [README.md](README.md) file located in the root
directory of this repo.

New developers should also read the [AIMS\_CODE\_GUIDELINES.md](AIMS_CODE_GUIDELINES.md)
file located in the root directory of this repo.  This file is a summary of best
practices and our general philosophy towards writing code in the FHI-aims
project.

## Using git with FHI-aims ##

FHI-aims tracks changes to the code using git, a version control tool which is
universally used for software development.  A free e-book describing git is
available at <https://git-scm.com/book/> .

   * Don't memorize the whole book.  Only the first few chapters (in version 2
     of the book: "Getting Starting", "Git Basics", and "Git Branching") are
     needed to understand how git works.

To make changes to FHI-aims, you must be using the git version downloaded from
our GitLab at <https://aims-git.rz-berlin.mpg.de/aims/FHIaims> .

When coding, remember: "Commit early, commit often."  It is highly
recommended that each individual git commit be small and have a well-defined
purpose:  it is not uncommon to commit three or more times a day when doing
heavy coding work.  By keeping commits small, it is easier to track changes
in your code over time as well as track down exactly when bugs were
introduced into the code.

Commit messages are important, because they serve as documentation of what
you did.  A commonly accepted format for commit messages is:

   * A single "title" line of 50 characters or less summarizing the git commit,
     following by
   * An empty line, and then
   * A "body" containing multiple lines, each line containing 72 characters or
     less.

This format ensures that your commit is readable, and many git utilities can
parse this format to make your message look "pretty".  Note that only the title
is mandatory; for small commits, the empty line and body is commonly omitted.
The title should be written in the imperative mood, i.e. use "Fix typo in
documentation" instead of "This commit fixes typo in documentation".

An example commit message written in this format is:

```
[CI] Remove absolute paths from all builds

This commit extends commit 7bb55dbe to all builds and has been verified
to fix issue #11.
```

## Submitting Code to the Main Repo ##

So you've made some changes to FHI-aims and you'd like to include them in future
releases.  Great!  Here's what you should do to make sure other users have
access to your code:

1. Run the `git status` command and check that all the files that you intend to
   distribute are indeed included.
    * This may seem petty, but the most common mistake that developers make is
      forgetting to add a file to a commit before pushing.  When this happens,
      only the original developer will have a copy of the updated file, breaking
      the code for all other developers.

2. Run the regression tests (located in the `regression_tests/` folder in this
   package) using an FHI-aims executable that you compiled based on your code
   changes.
    * The regression tests run a series of FHI-aims calculations using a
      user-provided executable and compare the results from that executable to
      previous "known-good" calculations to ensure that no functionality was
      broken (hopefully).
    * If you need help running the regression tests, don't hesitate to ask for
      help from other developers on the Slack channel.

3. Push your changes to the mainline repo via `git push`.

4. See if your push passed the continuous integration (CI) platforms, which test
   FHI-aims on a a variety of compilers, libraries, and build settings.
   FHI-aims currently employs two different CIs:
    * A Buildbot at <http://www.theo.ch.tum.de/bbot/>, which is maintained by
      the Reuter group at Technische Universität München, and
    * A GitLab CI, which is directly integrated into the FHIaims GitLab repo at
      <https://aims-git.rz-berlin.mpg.de/aims/FHIaims> .

## Adding New Functionality to FHI-aims ##

If you've added new functionality to FHI-aims, as opposed to tweaking existing
functionality, there are two additional steps you should do:

1. Update the manual (located in the `doc/` directory in this package) with
   instructions how to use your new functionality.  You are encouraged to
   shamelessly include references to papers of yours that people should cite if
   they're using your functionality.

2. Add a new test to the regression tests which tests the functionality you
   added.  Adding your own test to the regression suite has a number of
   advantages, including:
    * Other developers won't accidentally break your code,
    * Your code will be automatically tested by the continuous integration in
      FHI-aims, and
    * People will be reminded that your functionality exists every time they run
      the regression tests.

As always, feel free to contact other developers for more information how to
modify the manual and update the regression tests.
