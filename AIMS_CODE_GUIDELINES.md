# First Rule of Coding in FHI-aims, Above All Else: #

** NEVER copy other people's code into FHI-aims without their permission, in
violation of their license. **
```
  Other people work hard to create good scientific code. So do we.  We do not
  want anyone to take our code against our intentions.  So do they.

  It is not ok to take their code without telling them, putting it into
  FHI-aims, hoping that no one will notice. This is illegal and it is also not
  nice to your fellow scientists who wrote the code.  Please consider this.

  Never forget that there is always the option to talk to the original authors.
  Write to them and ask if they would be ok to see their work included into
  FHI-aims. Get written permission and include that permission in your version
  in FHI-aims, including attribution.

  When in doubt, please contact me (VB).  Normally, a solution can be found.
  Just, please do not ignore the problem.
```

In particular, everyone should understand that code under the GPL ("Free
Software") is NOT code without a license. In fact, the license is extremely
strict. Please read the license if you use such code.

# Where to Develop #

If you choose to develop any additional functionality, corrections, bugfixes, etc.
(please do! FHI-aims, at its heart, really is a community code - there is no strong
separation between "developers and users" - what gets developed does get developed
because someone needs it) - if you do choose to make any changes, please strongly
consider doing so in the main development version of FHI-aims, located at

https://aims-git.rz-berlin.mpg.de/ .

The main development version is routinely tested and usually very stable - many of
us use it in production. However, more importantly, making your changes there will make
sure that this improvement will remain available for yourself and others in the future.

# Commit Early, Commit Often #

Also, we encourage a strategy of "commit early, commit often" - i.e., consider not
creating a huge changeset, accumulated over months and years, and then committing all
at once. This will lead to code divergences, inability of others to review code changes
and point out improvements early, etc. Instead, try to break down your changes into
small, simpler chunks and push them into the main repository reasonably often. (Of
course, please do make sure that your changes do not break any existing functionality.)

Thank you!

# Code Conventions #

FHI-aims has never had strictly regulated set of "coding conventions". FHI-aims
is a product of its individual developers, users, and their enthusiasm. By
writing down a set of guidelines below, we do not wish to create any unnecessary
burdens. In the end, working code, and working science, counts. If you intend to
contribute code, that would be great. What we hope for is that any code can be
as cleanly structured and as well documented as you can - with reasonable effort
for you.

Thus, the following page contains a number of "guidelines" that code in FHI-aims
should follow. Most of them are motivated by a simple constraint - our desires
to keep even novice users able to compile a working, optimal version of FHI-aims
on their preferred platform. It is important to remember that, while we have
many great external tools available today, not all of them work on every
platform - and that especially the most expensive computers, those that do the
heaviest scientific work, do not always come with standardized, Gnu versions of
all external tools.

This document is also open for discussion and improvement. If you take issue
with any of these guidelines, or if you need to break one of them for clear
reasons, please talk about the issue in the forum. We will usually find a way.

# Some simple guidelines #

* Code that is linked to FHI-aims should be written in Fortran, if that is
possible.
```
  Reason: We are certainly aware of the virtues of C or C++ especially when it
  comes to dealing with system calls or similar. However, remember that FHI-aims
  should remain compilable even by users who might only have a very rudimentary
  idea of what a compiler is, let alone where it lives on their system. But
  every single new dependency in a Makefile multiplies the probability that all
  goes well by (say) ~1/2. So keeping things as simple as possible (but not
  simpler) can save someone else's day.

  The truth is that Fortran has been called outdated by some for 20+ years, and
  yet it is still as good for scientific computing as it ever was. It simplifies
  a few things; it is also very flexible. Keep an open mind. The real challenge
  in scientific computing is not the programming language. The real challenge is
  to understand the science before attempting to read or write the code. (This
  is not trivial for any of us.)
```

* Please use conservative code elements if at all possible.
```
  Reason: We understand about object oriented programming, new language
  features, and updated code standards.  However, language features do not
  automatically become usable by every compiler just because they are written
  down in a recent standard. Many current compilers, especially those on the
  largest and most powerful platforms, may not support every new language
  feature in the exact same way, and what happens is that the latest and
  greatest feature that works on one machine prevents someone else from using
  another machine altogether, 2 months later and halfway around the globe.

  What we do in science, though, is mainly "Formula Translation". We add and
  subtract numbers according to a recipe and those are integer, real, and
  complex numbers of a certain and usually entirely predictable type. Using
  intrinsic types that are supported by all compilers (like real*8) is therefore
  totally fine and will perform very well on most platforms. In fact, many
  compiler optimizers will recognize such variables and, in case of doubt,
  produce faster code.  However, if that number is encapsulated in some
  complicated new type construct it is not at all guaranteed that any compiler
  will understand how to produce optimized code for that.

  In FHI-aims, we use Fortran 2003 and that is good for pretty much everything.
  Fortran 2008 is also around, but is in fact much newer and can lead to problems.
  Although 2008 sounds like a long time ago, it isn't indicative of the year
  the standard went into effect, however, or the time it took to shake out most
  compiler bugs on all compilers. In fact, as of May 2019, it is still the case
  that *most* Fortran compilers do not fully support the Fortran 2008 standard.

  If you critically need a new language feature, please discuss, but if at all
  possible, try to find a way to work with the traditional language elements
  used in numerical programming. This may end up saving someone else a lot of
  time.
```

* If you need a dependency on an external library, please make it optional, if
at all possible. Follow the existing model of "stubs" that we provide for
external dependencies, and possibly a separate Makefile.
```
  Reason: Again, we are trying to keep the compilation process of FHI-aims as
  simple as possible, in terms of the prerequisites needed.  As it is, we have
  two definite library dependencies that we can not circumvent:

  (1) Lapack
  (2) BLAS

  because these libraries are so critical for numerical efficiency that we can
  not afford to just compile them ourselves.  Furthermore,

  (3) MPI is practically indispensable, as is the infrastructure provided by
  (4) Scalapack
  (5) BLACS

  Yet, even (3)-(5) can be circumvented by using appropriate targets in the
  Makefile. In that case, the library calls will be replaced by 'stub' routines
  supplied by ourselves. Thus, a serial version of FHI-aims will never call the
  'stub' routine in question, and all that the stub routine does is complain and
  stop if it gets called by mistake.

  Why do we not use ifdef's?

  First, because with the stub model, we can still control entirely at runtime
  (by if's) what we do.  This is often more intuitive.

  Second, because "-D" means very little to most users. Plenty of "-D"'s just
  make it yet much harder to understand the variables needed to compile the
  code.

  Many of our developers have added excellent and valuable external library
  calls to the code, padded with stub libraries that can serve as a great
  example of how to do this. What we need to do next is to limit the growth of
  too many disconnected Makefiles by a more transparent infrastructure. But it
  is actually a virtue that, at this point, one can compile a down-to-earth
  version of FHI-aims without having to worry about those extra dependencies at
  all, if needed.
```

* No "ifdef"s in the code. See above.
```
  Yes, we know what ifdef's are and what a preprocessor is. They were already
  old technology when even the oldest FHI-aims developers were young. The
  decision to avoid ifdef's was deliberate - see above - and this actually
  worked. Please try without them.
```

* Use meaningful, reasonably verbose variable names for everything you write.
Yes, even for counters.
```
  Using verbose name means more typing, but it helps others understand your
  code. Knowing that your loop runs over i_grid, i_basis is eventually so much
  more meaningful than having it run over i, j. But there is another reason.
  Have you ever tried to use a search function to search for the variable called
  'n'?

  We have some code where authors did not choose to use long names in FHI-aims,
  but that code remains hard to read, in the opinion of this writer (VB), even
  after working with it for a long time. If possible, use variable names that
  are intuitive for a physicist. If it sounds like a textbook quantity, chances
  increase dramatically that someone with textbook knowledge will understand it.
```

* Please use "implicit none" at the beginning of each file.
```
  Reason: It is always safer to declare variables properly. Saves you from type
  mismatches, typos, uninitialized variables, etc. And even more importantly, it
  may save future developers down the road from running into trouble with your
  code because a new compiler suddenly exposes an issue with your code.
```

* Always specify which variables you are importing from a module via
`use module_name, only:` statements.
```
  Reason: Declaring which variables you're importing from a module makes the
  code structure easier to understand for other developers by documenting the
  implicit interface introduced by "use" statements, similar to how "implicit
  none" statements makes code easier to understand by forcing developers to
  declare variables within a subroutine.
```

* Please write properly indented code.
```
  Reason: It just makes the code so much more legible. However, let people use
  the indentation scheme that they prefer. Arguing over two indents vs three
  intents vs four will not help us do great science.

  That being said, whenver you are modifying pre-existing files, it is
  strongly encouraged that you use whatever indentation scheme was already in
  use in the file.  Otherwise, over time, the file will become difficult to
  read thanks to frequent, jarring style changes.
```

* Please do not exceed line lengths of 130 characters.
```
  The gfortran compiler enforces this length by default for reasons that are not
  entirely clear. Of course, there is a documented gfortran flag that all of us
  know about to fix the issue, but most new users do no know this. If they
  encounter a line that is too long, the impression is "FHI-aims does not
  compile".  After years of this, the decision was to please keep lines
  reasonably short. This also helps developers who do not work on extra-wide
  editors.

  In fact, 80 characters is not that bad for most purposes, and most coding
  guidelines recommend using 80 characters whenever possible.
```

* If you do cleanups to an entire file, be careful to understand the
consequences.
```
  Reason: If you (for example) add a single space to every line of a source
  file, the 'git' version control system is conservative and concludes that
  every single line has changed. Now imagine someone else who has uncommitted
  modifications to that file. For better or for worse, git will no longer
  silently (and correctly) merge their changes. Instead, everyone will have to
  merge their changes by hand, which creates a significant risk of introducing
  unintended errors. So, while code cleanup is a good thing, please be
  conservative.
```

* Please refrain from using pointers.
```
  Reason: Most of what pointers do can be done in Fortran without them. However,
  pointers have the unpleasant side effect of being able to bring about memory
  leaks, something that is normally not a problem of Fortran. They can also
  bring out subtle compiler bugs. There are still a very few pointers that got
  into FHI-aims over the objections of the author of these lines (VB), for the
  simple reason that they did obvious and lasting good. But the compiler issue
  is real, and has struck us even with these very few pointers.
```

* For parallel implementations, please use MPI.
```
  Reason: MPI is still the most generic parallel infrastructure out there, and
  much of why FHI-aims is fast comes down to someone having figured out a way
  how to organize memory and communication explicitly and optimally - which is
  exactly what MPI forces you to do. Many basic communication interfaces of our
  own already exist in the synchronize_mpi.f90 and synchronize_mpi_basic.f90
  modules. The code contains plenty of excellent examples how these can be used.
```

* Please use the variable 'use\_unit', not '6' or `'*'` for output. If possible,
please use the infrastructure provided in `localorb_io.f90` .
```
  Reason: Multiple instances of FHI-aims can actually be used side by side as
  library routines, side by side under a parallel wrapper. Each instance gets
  its own MPI communicator and does not know of the others. This is ongoing work
  with great contributions from Matt Farrow and also Andrew Logsdail in London,
  and makes the code that much more powerful. However, the standard output
  streams of each instance of FHI-aims must then be kept separate. That is
  exactly what 'use_unit' does.
```

* Please write ALL comments and variable names in English - even when you never
intend to commit that code in question.
```
  Reason: Much code that is never intended to be shared turns out to be great,
  and extremely useful. Pushing it to the mainline is suddenly essential. Yet,
  you have that paper to write and that degree to finish at the same point in
  time. No one, not even with the best of intentions, will have the time to
  change all those variable names and comments written in a Frankonian dialect.
  Using English right away is simpler.
```

* Use the `aims_stop()` subroutine (located in the `mpi_tasks` module) to stop
an FHI-aims calculation gracefully.
```
  Reason:  In parallel calculations, the Fortran built-in STOP statement
  called on any MPI task will kill the calculation across all MPI tasks.  When
  many MPI tasks are being used, and they're all at different stages of the
  calculation, the abrupt nature of the STOP statement will often lead to
  output being abruptly cut off, making it difficult to understand and debug
  what went wrong.

  To alleviate this, we have provide the aims_stop() subroutine which allows
  the developer to output an error message, ensures that all output is properly
  accounted for, and gracefully exits the calculation.
```

* Use meaningful error messages, especially when the code needs to stop.
```
  Reason: Having FHI-aims stop with a message such as "VB - Error 10" will make
  it very hard for anyone else to understand what happened. If you create a code
  path where a stop is needed (for example, due to conflicting input requests) -
  please take the time to explain in physical terms what is needed, and what the
  user could do to prevent this error.
```

* Do not use static (a.k.a. `save`) variables.
```
  Reason: Static/save variables in subroutines should be avoided whenever
  possible, as they are a legacy feature inherited from older Fortran standards
  which can introduce bugs that are difficult to find.  As of Fortran 95, it is
  far more preferable to place variables whose state needs to be saved into a
  module.

  One important thing to note is Fortran will implicitly mark any variable that
  is initialized when it is declared as a static variable whether or not the
  save keyword is used, e.g. statements of the form:

  integer :: i = 5

  Instead, it is better to split the declaration and initialization onto
  different lines:

  integer :: i
  i = 5
```

* Use a subroutine to initialize default values for modules.
```
  Reason:  Many people choose to run FHI-aims as a library subroutine called
  as part of a workflow, rather than as a stand-alone executable.  When
  running FHI-aims as a library subroutine, the state of the calculation
  between each call of FHI-aims needs to be reset, otherwise information from
  previous calculations will leak into future calculations.

  Fortran provides the ability to set default values for module variables when
  they're declared, however this only applies to the first time FHI-aims is run
  as a library subroutine.  Subsequent runs will use the value from the previous
  run and not the specified default value.

  To overcome this behavior, we strongly recommend that all developers use
  subroutines to initialize default values for module variables.  For those of
  you familiar with object-oriented programming, think of this as writing a
  default constructor for the module.

  To help with setting default values for modules, there is a
  set_aims_defaults() subroutine which is called at the beginning of an
  FHI-aims calculation and serves as a wrapper around other subroutines to set
  default values for module variables.  You may choose to set defaults in a
  different location, depending on what works for you, but the important thing
  is that you do it somewhere.
```

* Deallocate all allocatable arrays defined in modules at the end of every
FHI-aims calculation.
```
  Reason:  When running FHI-aims as a library subroutine, allocatable arrays
  that are defined in modules will not be deallocate automatically at the end
  of the FHI-aims calculation.  This creates a memory leak and has the potential
  to change the behavior of subsequent FHI-aims calculations.

  We strongly recommend that all developers write a subroutine to deallocate
  all allocatable arrays within a module.  For those of familiar with
  object-oriented programming, think of this as writing a destructor for the
  module.  (It's also not a bad idea to set the values for module variables back
  to their default values.)

  To help with deallocating allocatable arrays for modules, there is a
  final_deallocations() subroutine which is called at the end of an FHI-aims
  calculation and serves as a wrapper around other subroutines to deallocate
  allocatable arrays for modules.  You may choose to deallocate arrays in a
  different location, depending on what works for you, but the important thing
  is that you do it somewhere.
```

* Use the subroutines in the `aims_memory_tracking` module to allocate and
deallocate memory.
```
  Reason:  A common misperception about scientific computing is that we are
  primarily concerned with the runtime of our calculations.  While this was true
  in the past and continues to be an important factor in the present,
  increasingly our major concern when writing code is the memory consumption, as
  this often has a greater effect in limiting the achievable sizes of our
  calculations.  (As well as the parallelism of our algorithms, but that's a
  topic for another time.)

  We have written an aims_memory_tracking module which we use to estimate the
  memory overhead for the FHI-aims calculations by tracking memory allocations
  and deallocations using aims_allocate() and aims_deallocate() subroutines.
  We strongly encourage all developers to use this framework, as it has aleady
  proven extremely useful in tracking down memory bottlenecks and potential
  memory leaks.
```

* Be practical. Do not let yourself be bogged down by some rule.
```
  The code guidelines above are all there because they have made evident sense
  in at least one particular context, or else they would not be there. However,
  if one of them absolutely does not make sense in your case - talk to us.
  Chances are you might be right, if only for your specific needs. We would much
  rather have great code that changes science than code that follows rules for
  the sake of looking good.
```
