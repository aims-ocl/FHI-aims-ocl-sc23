### Summary

(Summarize the bug encountered concisely.  This should be something that developers can quickly glance at to understand the basic details of your bug report, no more than a paragraph.)

### Build environment

(Many of the bugs we encounter in HPC applications like FHI-aims are the result of specific versions of compilers and/or libraries.  Thus, it is important that we know as much as possible about how your executable was built and on what computing resource you are running.  A list of important compilation choices includes:

* Commit of FHI-aims used:  specify what git commit of FHI-aims was used, and whether you've modified it.  You may find this information at the top of FHI-aims' output in the "Git rev." section of the following code block, such as shown here:
```
Git rev. (modified): ba8e920 sum over states approximation for lorentz [...]
```
* Computing resource used: e.g. MPG's Hydra, OLCF's Titan, ALCF's Theta, my Linux Mint 18.3 laptop, my group's/university's compute cluster
* Compiler used:  e.g. Intel 2018.0.1, Cray 2.5.12, IBM XL 13.1.5, GNU 7.1.0, PGI 17.9.0
* MPI libraries used: e.g. Intel MPI 2018.0.1, Cray MPICH 7.6.3, Open MPI 1.7.5, MPICH 3.1.4
* Math libraries used: e.g. Intel MKL 2018.0.1, Cray LibSci 16.11.1, ATLAS 3.10.3

Please try to include as much of this information as possible, but feel free to leave it out things you don't know.  We would rather have a partial bug report than no bug report at all.  Alternatively, if you have already narrowed the bug down to a specific compilation choice, feel free to list only that compilation choice and the steps you took to determine it was responsible for the bug.
)

### Steps to reproduce

(How one can reproduce the issue - commonly, this will include information about the compilation settings for FHI-aims and keywords in your control.in file.  If you are having difficulty figuring out what the steps to reproduce may be, a list of things you may want to think about includes:

* Can the bug be fixed by recompiling FHI-aims in a different way?  (See previous section for examples)
* Did this bug only appear when you upgraded to a newer version of FHI-aims?
* Are you running a non-periodic (molecular) or periodic (solid) calculation?
* How many atoms does your system contain?  Is it a very small calculation (2 atoms) or a large calculation (1000s of atoms)?
* How many MPI tasks are you using?
* If you are running on an HPC platform, how many nodes are you using?
* If you are running a periodic calculation, how many k-points are you using?
* Which basis set are you using?  (Tier 1, Tier 2, etc.)
* Are you doing anything beyond a ``standard'' SCF cycle? (e.g. spin-orbit coupling, band structures, GW, implicit solvation, NMR calculations, etc.)

Please only report the bare minimum set of steps necessary to reproduce the bug.
)

### What is the current *bug* behavior?

(What actually happens)

### What is the expected *correct* behavior?

(What you should see instead)

### Relevant logs and/or screenshots

(Please include the output of your aims executable as an attachment to this bug report, as this will allow us to better reproduce the bug.  If you would like to insert text here, please use code blocks (```) to format console output, logs, and code as it's very hard to read otherwise.)

### Possible fixes

(Feel free to speculate about what you think may be causing the problem here.  If you can, link to the file and line of code that might be responsible for the problem)

/label ~Bug
