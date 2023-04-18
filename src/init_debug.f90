!****s* FHI-aims/init_debug
!  NAME
!   init_debug
!  SYNOPSIS
      subroutine init_debug()
!  PURPOSE
!     This function helps you to define a debug module. After defining your
!     module here with
!
!        register_debugmodule("your_tag")
!
!     Afterwards, you can conveniently activate debugging in the inputfile by
!     using the flag:
!
!        debug_module your_tag
!
!     This tells the code to enable debugging for your code. Just call
!
!        debugprint(message, your_tag)
!
!     in your code to only print the message if debugging for your module is
!     enabled. For more complicated debug functionality, there also is
!
!        module_is_debugged(your_tag)
!
!     This function returns a logical value and thus can be used as condition
!     in an enclosing if-block.
!
!     That's already everything you need to do to enable conditional debugging
!     in your modules!
!  USES
      use debugmanager, only: init_debuglist, register_debugmodule
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2014).
!  SOURCE
      implicit none

      call init_debuglist()
      call register_debugmodule("DFPT")
      call register_debugmodule("friction")

      end subroutine init_debug
