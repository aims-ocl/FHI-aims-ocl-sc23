!****h* FHI-aims/applicable_citations
!  NAME
!    timing
!  SYNOPSIS
      module applicable_citations
!  USES
      implicit none
!  PURPOSE
!    This module gathers citations that apply to a given FHI-aims run
!    along the way and writes them out at the end of the code.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu, 
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio molecular simulations with numeric atom-centered orbitals"
!    Computer Physics Communications 180, 2175-2196 (2009).
!    http://dx.doi.org/10.1016/j.cpc.2009.06.022
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2014).
!  SOURCE
!     global variable declarations - exported to other program parts
!     these are the individual time stamps as collected throughout the code

      ! all module-wide symbols following this statement will not
      ! be made public

      logical, public :: warn_citation_unknown

      public :: init_citationlist, cite_reference, register_citation, &
         print_references, cleanup_citationlist

      private

      type citation
         character*150 :: tag
         character(len=130) :: text(15)
         integer :: lines
         integer :: priority
         logical :: active
      end type citation

      type(citation), allocatable :: citelist(:)
      integer :: citecount
      character(len=2000) :: info_str

      contains
!******	

!****s* applicable_citations/init_citationlist
!  NAME
!    init_citationlist
!  SYNOPSIS
      subroutine init_citationlist()
!  PURPOSE
!
!  INPUT
      use mpi_tasks, only: check_allocation
      implicit none
!     none
!  OUTPUT
!     none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  SEE ALSO
!    FHI-aims aimsclub Wiki.
!  SOURCE
      integer :: info
      character(*), parameter :: funcname = "init_citationlist"

      if (allocated(citelist)) then
         deallocate(citelist)
      end if
      allocate(citelist(100), stat=info)
      call check_allocation(info, funcname)
      citecount = 0
      warn_citation_unknown = .false.

      end subroutine

!****s* applicable_citations/register_citation
!  NAME
!    register_citation
!  SYNOPSIS
      subroutine register_citation(tag, reference, priority)
!  PURPOSE
!  USES
      use mpi_tasks, only: aims_stop
      implicit none
!  INPUT
      character(len=*), intent(in) :: tag
      character(len=130), intent(inout) :: reference(15)
      integer :: lines
      character(len=*), intent(in) :: priority
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  SEE ALSO
!    FHI-aims aimsclub Wiki.
!  SOURCE
      character(*), parameter :: funcname = "register_citation"
      integer :: i, j
      type(citation) :: newcite

      select case (priority)
         case ("core")
            newcite%priority = 0
         case ("feature")
            newcite%priority = 1
         case ("technical")
            newcite%priority = 2
         case default
            write(info_str,'(1X,A,A)') &
            "* Error - citation ", tag, " requested an unknown priority!"
            call aims_stop(info_str, funcname)
      end select
      do i = 1, 15
         if (trim(reference(i)) == "____") exit
         newcite%lines = i
      enddo
      newcite%tag = tag
      newcite%text(:) = reference(:)
      newcite%active = .false.
      do i = 1, citecount+1
         if ((newcite%priority < citelist(i)%priority) .or. &
         (newcite%priority == citelist(i)%priority .and. &
         newcite%tag < citelist(i)%tag) ) then
            do j = citecount+1, i+1, -1
               citelist(j) = citelist(j-1)
            enddo
            citelist(i) = newcite
            exit
         elseif (i == citecount+1) then
            citelist(i) = newcite
         endif
      enddo
      citecount = citecount + 1
      reference(:) = "____"

      end subroutine register_citation

!****s* applicable_citations/cite_reference
!  NAME
!    cite_reference
!  SYNOPSIS
      subroutine cite_reference(tag)
!  PURPOSE
!
!  INPUT
      use localorb_io, only: localorb_info
      implicit none
      character(len=*), intent(in) :: tag
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  SEE ALSO
!    FHI-aims aimsclub Wiki.
!  SOURCE
      integer i

      do i = 1, citecount
         if (citelist(i)%tag == tag) then
            citelist(i)%active = .true.
            exit
         elseif (i == citecount) then
            write(info_str,'(1X,A)') &
            "* Warning - a code part requested a citation which was not known"
            call localorb_info(info_str)
            write(info_str,'(1X,A)') &
               "* please add it to init_citations.f90."
            call localorb_info(info_str)
            write(info_str,'(1X,A,A)') "* Missing reference tag: ", tag
            call localorb_info(info_str)
            write(info_str,'(A)') ""
            call localorb_info(info_str)
            warn_citation_unknown = .true.
         endif
      enddo

      end subroutine cite_reference

!****s* applicable_citations/print_references
!  NAME
!    print_references
!  SYNOPSIS
     subroutine print_references()
!  PURPOSE
!
!  USES
      use localorb_io, only: localorb_info
      implicit none
!  INPUT
!  none
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  SEE ALSO
!    FHI-aims aimsclub Wiki.
!  SOURCE
      integer :: i, j

      call localorb_info("  Methods described in the following list of references were used in this FHI-aims run.")
      call localorb_info("  If you publish the results, please make sure to cite these reference if they apply.")
      call localorb_info("  FHI-aims is an academic code, and for our developers (often, Ph.D. students")
      call localorb_info("  and postdocs), scientific credit in the community is essential.")
      call localorb_info("  Thank you for helping us!")
      !call localorb_info(info_str)
      do i = 1, citecount
         if (citelist(i)%active) then
            call localorb_info("")
            do j = 1, citelist(i)%lines
               write(info_str,'(A)') trim(citelist(i)%text(j))
               call localorb_info(info_str)
            end do
            call localorb_info("")
         endif
      enddo
      call localorb_info("  Of course, there are many other important community references, e.g., those cited in the")
      call localorb_info("  above references. Our list is limited to references that describe implementations in the")
      call localorb_info("  FHI-aims code. The reason is purely practical (length of this list) - please credit others as well.")
      !call localorb_info(info_str)

      end subroutine print_references

!****s* applicable_citations/cleanup_citationlist
!  NAME
!    cleanup_citationlist
!  SYNOPSIS
      subroutine cleanup_citationlist()
!  PURPOSE
!
!  INPUT
      implicit none
!     none
!  OUTPUT
!     none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  SEE ALSO
!    FHI-aims aimsclub Wiki.
!  SOURCE
      integer :: info
      character(*), parameter :: funcname = "cleanup_citationlist"

      if (allocated(citelist)) deallocate(citelist)
      end subroutine
!******

   end module applicable_citations
