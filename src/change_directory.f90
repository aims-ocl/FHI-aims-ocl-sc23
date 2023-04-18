!****s* FHI-aims/change_directory
!  NAME
!    change_directory -- Subroutine to ti interface against the C based change
!                        of current directory
!  SYNOPSIS

   subroutine change_directory(path)

!  PURPOSE
!
!    This is FHI-aims code
!
!    In principle, FHI-aims can also be built as a library. In that case, aims
!    is called as a subroutine aims() (see main.f90). The present program embeds
!    this subroutine for the standalone executable.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUT
!    o path   The path which schould be set to the current directory
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!     Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement.
!  SOURCE
    
!  USES
      use localorb_io
      use mpi_tasks, only: aims_stop

!  VARIABLES
      implicit none
      
      !Arguments
      character(*),  intent(in)  :: path

      ! Local variables
      integer        :: error
      character*200  :: message
  
      interface
         subroutine change_directory_c(path,error) &
            bind(C,name="change_directory_c")
            use, intrinsic :: iso_c_binding 
            implicit none
            
            character(kind=C_CHAR), dimension(*), intent(in)  :: path
            integer(kind=C_INT),                  intent(out) :: error
         
         end subroutine change_directory_c
      end interface

!  START WORK

      ! Remove leading and following blanks and add terminating
      ! /0 character for C compatibility
      call change_directory_c(TRIM(ADJUSTL(path))//CHAR(0), error)

      if (error /= 0) then
         write(message,'(2A)') "Could not change to directory : ", &
            TRIM(ADJUSTL(path))
         call aims_stop(message,"change_directory.f90")
       end if

    end subroutine change_directory
