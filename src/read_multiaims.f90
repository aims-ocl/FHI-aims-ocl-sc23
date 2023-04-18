!****s* FHI-aims/read_multiaims
!  NAME
!    read_control
!  SYNOPSIS

subroutine read_multiaims (tasks_per_subjob, start_id)

!  PURPOSE
!    This is the central subroutine that reads and checks all information found
!    in input file "multiaims.in". 
!    The centerpiece is simply a line-by-line reading
!    of multiaims.in, with a large "if" statement deciding what to do if a given
!    keyword is found at the beginning of a line.
!
!  USES

   use localorb_io
   use mpi_utilities

!  INPUTS
!    none
!  OUTPUT
!    o tasks_per_subjob  Number of tasks per subjob
!    o start_id          Id of job to start with
!  AUTHOR
!    FHI-aims team.
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!  e.V. Please note that any use of the "FHI-aims-Software" is subject to
!  the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

   implicit none

   ! Arguments
   integer, intent(out) :: tasks_per_subjob
   integer, intent(out) :: start_id

   ! local variables
   character*140 :: info_str
   character(*), parameter :: func = 'read_multiaims'

   !  begin work

   call localorb_info('',use_unit)
   call localorb_info( &
   "------------------------------------------------------------")
   call localorb_info("          Reading file multiaims.in.")
   call localorb_info( &
   "------------------------------------------------------------")

   call initialize
   call readfile
   call checkinput

   call localorb_info( &
   "------------------------------------------------------------")
   call localorb_info('')

   return

 contains

   subroutine initialize

      tasks_per_subjob = -1
      start_id         = -1

   end subroutine initialize



   subroutine readfile
      implicit none

!  local variables

!  i_code : Flag to determine i/o status
!  desc_str: line descriptor

      integer :: i_code, linecount

      character*132 inputline
      character*40 desc_str

      open (7, FILE="multiaims.in")
      linecount = 0
      do

         read(7,'(A)',iostat=i_code) inputline
         if(i_code<0) exit               ! end of file reached
         if(i_code>0) then
            call aims_stop_coll("Unknown error reading file 'multiaims.in'", &
                  func)
         endif

         linecount = linecount + 1

         read(inputline,*,iostat=i_code) desc_str
         if(i_code/=0) cycle              ! skip empty line

         if (desc_str(1:1).eq.'#') cycle  ! skip comment

!     decode contents of current line

         select case (desc_str)

         case('tasks_per_subjob')
            read(inputline,*,end=88,err=99) desc_str, tasks_per_subjob
            write(info_str, '(2X,I5,A)') &
               & tasks_per_subjob, " tasks per subjob requested."
            call localorb_info(info_str)

         case('start_id')
            read(inputline,*,end=88,err=99) desc_str, start_id
            write(info_str, '(2X,I5,A)') &
               & start_id, " is first job in line."
            call localorb_info(info_str)
        
         case default
            write (info_str,*) &
               "Syntax error reading 'control.in' (unknown keyword)"
            call localorb_info(info_str)
            write (info_str,*) "line: '"//trim(inputline)//"'"
            call localorb_info(info_str)
            call aims_stop_coll('UNKNOWN KREYWORD in multiaims.in', func)

         end select

      end do
 
!  close multiaims.in

      if (linecount == 0) then
         call aims_stop_coll('Attention - empty input file multiaims.in.', func)
      end if

      close(7)

      write(info_str,'(2X,A)') &
      "Finished reading input file 'multiaims.in'."
      call localorb_info(info_str)
      write(info_str,'(2X,A)') &
      "Consistency checks are next."
      call localorb_info(info_str)

      return

   88 continue
      write (info_str,*) &
         "Syntax error reading 'multiaims.in' (missing arguments)"
      call localorb_info(info_str)
      write (info_str,*) "line: '"//trim(inputline)//"'"
      call localorb_info(info_str)
      call aims_stop_coll('MISSING ARGUMENT in multiaims.in',func)

   99 continue
      write (info_str,*) "Syntax error reading 'control.in'"
      call localorb_info(info_str)
      write (info_str,*) "line: '"//trim(inputline)//"'"
      call localorb_info(info_str)
      call aims_stop_coll('READING ERROR in multiaims.in', func)

   end subroutine readfile


   subroutine checkinput
   
     implicit none

    !  check whether relevant quantities have been set

    if (tasks_per_subjob.eq.-1) &
       call aims_stop(' Tasks per subjob not specified. What shall I do?')
    if (tasks_per_subjob.gt.n_tasks) &
       call aims_stop(' Tasks per subjob is greater than total tasks.')

    if (start_id.eq.-1) &
       call aims_stop(' Start ID for job not specified. Where shall I start?')

    write(info_str,'(2X,A)') &
       "All input parameters seem consistent."
    call localorb_info(info_str)
   end subroutine checkinput


end subroutine read_multiaims
!******
