!****s* FHI-aims/read_specified_grid
!  NAME
!   read_specified_grid
!  SYNOPSIS

   subroutine read_specified_grid ( i_species )

!  PURPOSE
!  Reads user-specified angular grid settings for a species.
!
!  USES

     use dimensions
     use constants
     use species_data
     use localorb_io
     use mpi_tasks, only: aims_stop
     implicit none

!  ARGUMENTS

     integer :: i_species

!  INPUTS
!   o i_species -- species where the grid belongs
!  OUTPUT
!    none
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
!    Release version, FHI-aims (2008).
!  SOURCE





     ! Local variables

     integer :: i_separator
     integer :: i_code
     character*132 :: inputline

     logical :: flag_outer

     integer :: n_outer_points

     character*20 :: desc_str

     character*100 :: info_str

     ! Begin work

     i_separator = 0
     flag_outer = .false.

     lineloop: do

       read (7,'(A)',iostat=i_code) inputline

       if (i_code<0) then
          backspace(7)
          exit lineloop
       endif
       if (i_code>0) then
          backspace(7)
          return
       endif

       read (inputline,*,iostat=i_code) desc_str
       if (i_code/=0) cycle lineloop ! empty line
       if (desc_str(1:1).eq."#") cycle lineloop ! comment

       if (desc_str.eq."division") then
         i_separator = i_separator + 1

         read(inputline,*,end=88,err=99) desc_str, &
           r_ang_shell(i_separator,i_species), n_ang_points(i_separator,i_species)

       else if (desc_str.eq."outer_grid") then
         read(inputline,*,end=88,err=99) desc_str, n_outer_points
         flag_outer = .true.

       else
         ! unknown keyword - must be the end of the grid specification
         ! read keyword again in read_species_data

         backspace(7)
         exit lineloop

       end if

     enddo lineloop

     ! Verify input grid consistency

     if (flag_outer) then
       n_ang_points(i_separator+1,i_species) = n_outer_points
     else
       if (i_separator.ge.1) then

         write(info_str,'(2X,A)') &
         "| Angular grid for outermost region not specified."
         call localorb_info(info_str)
         write(info_str,'(2X,A,I5,A)') &
         "| Defaulting to ", n_ang_points (i_separator, i_species), " points." 
         call localorb_info(info_str)

         n_ang_points(i_separator+1,i_species) = n_ang_points(i_separator,i_species)

       else

         write(info_str,'(1X,A)') & 
         "* Error: Angular grid specification not found. Please correct."
         call localorb_info(info_str)
         call aims_stop ("No angular grid specified","read_specified_grid") 

       end if
     end if

     n_ang_shells (i_species) = i_separator+1

     ! check whether the order of radii given makes sense
     do i_separator = 2, n_ang_shells(i_species)-1, 1
       if ( r_ang_shell(i_separator-1,i_species) .gt. r_ang_shell(i_separator,i_species) ) then
         
         write(info_str,'(1X,A,I3,A,I3,A)') & 
         "* Error: Incorrect order of grid divisions ", i_separator-1, ", ", i_separator, "."
         call localorb_info(info_str)
         call aims_stop("Incorrect order of grid divisions",&
                        "read_specified_grid") 

       end if
     enddo

     ! convert to Bohrs
     do i_separator = 1, n_ang_shells(i_species)-1, 1
       r_ang_shell(i_separator,i_species) = r_ang_shell(i_separator,i_species) / bohr
     enddo

     write(info_str,'(2X,A,I5,A)') &
     "| Specified grid contains ", n_ang_shells (i_species), " separate shells." 
     call localorb_info(info_str)
     write(info_str,'(2X,A)') &
     "| Check grid settings after all constraints further below." 
     call localorb_info(info_str)

     return

  88 continue
     write(use_unit,*) "Syntax error reading 'control.in' (missing arguments)"
     write(use_unit,*) "line: '"//trim(inputline)//"'"
     stop

  99 continue
     write(use_unit,*) "Syntax error reading 'control.in'"
     write(use_unit,*) "line: '"//trim(inputline)//"'"
     stop

   end subroutine read_specified_grid
!******	
