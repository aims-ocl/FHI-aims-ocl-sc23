!****s* FHI-aims/generate_aims_uuid
!  NAME
!   generate_aims_uuid
!  SYNOPSIS
!
module generate_aims_uuid
!  PURPOSE
!  Generates a random number based unique Identifier for every FHI-aims run
!
!  AUTHOR
!    Lydia Nemec, FHI-aims team, Chair for Theoretical Chemistry, Technical University Munich
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180, 2175-2196 (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2016).
!
!  Uses

   implicit none

!  SOURCE

! By default all variables below are private
   private

!  PURPOSE
!  The Routine generate a unique Identifier for ever FHI-aims run
!  Ideally the the aims_uuid will be written in any output file
!  Allowing to sort different Aims output files to a specific
!  FHI-aims run.
!  The identifier is given in the usual format according to RFC 4122
!  type 4 Standard with 5 hexadecimal numbers as a 37 character long string
!  where 4 stands for the uuid type and 'a' is a reserved character
!  ********-****-4***-a***-************
!  Details can be found at https://tools.ietf.org/html/rfc4122
!  The remaining numbers '*' are set by random numbers.
!  The seed for the pseudo random number generator is generated
!  by urandom where available or alternatively by the CPU-time.
!
   character(LEN=37) :: aims_uuid
   logical :: FirstCall = .TRUE.

   public :: write_aims_uuid
   public :: write_bare_aims_uuid

   contains

   integer function lcg(s)

      implicit none
      integer (kind=8) :: s

      if (s == 0) then
         s = 104729
      else
         s = mod(s, int(huge(0_2)*2,kind=8))
      end if

      s = mod(s*int(huge(0_2),kind=8), int(huge(0_2)*2,kind=8))
      lcg = int(mod(s, int(huge(0),kind=8)), kind(0))

   end function lcg

   integer function newunit(unit) result(i_io)

     implicit none

     integer, intent(out), optional :: unit

     ! local
     integer, parameter :: io_min=10, io_max=100
     logical :: unit_open

     ! begin

     do i_io=io_min,io_max
        inquire(unit=i_io,opened=unit_open)
        if (.not. unit_open) then
          if (present(unit)) then
            unit=i_io
          end if
          return
        end if
     end do
  
     i_io = -1

   end function newunit

   subroutine write_aims_uuid(info_str)

     implicit none

     character(LEN=*), intent(inout) :: info_str
     integer :: i

     if (FirstCall) then
       FirstCall = .FALSE.
       call gen_aims_uuid(aims_uuid)
     end if

     write(info_str,'(A,X,A)') "aims_uuid :", aims_uuid

   end subroutine write_aims_uuid

   ! This is functionally identical to write_aims_uuid, but the string it
   ! returns contains only the 36-character uuid plus padding space
   subroutine write_bare_aims_uuid(info_str)

     implicit none

     character(LEN=*), intent(inout) :: info_str
     integer :: i

     if (FirstCall) then
       FirstCall = .FALSE.
       call gen_aims_uuid(aims_uuid)
     end if

     write(info_str,'(A)') trim(aims_uuid)

   end subroutine write_bare_aims_uuid


   subroutine gen_aims_uuid(aims_uuid)
      use mpi_tasks
      implicit none
      character(LEN=37) :: aims_uuid

      integer :: i3, i4
      real :: r(8)
      integer :: i_entry, error
      character (LEN=3) :: s3
      character (LEN=4) :: s4

      i3 = 4095
      i4 = 65535

      call init_random_seed()
      call RANDOM_NUMBER(r)

      if (myid == 0) then
         do i_entry=1,8
            write(s3,"(Z3.3)") TRANSFER(int(r(i_entry)*i3),16)
            write(s4,"(Z4.4)") TRANSFER(int(r(i_entry)*i4),16)
            if (i_entry == 1) then
               write(aims_uuid,'(A)') s4
            else if (i_entry == 2) then
               write(aims_uuid,'(A,A)') trim(aims_uuid), s4         
            else if (i_entry==3) then
               write(aims_uuid,'(A,A,A)') trim(aims_uuid), '-', s4
            else if (i_entry==4) then        
               write(aims_uuid,'(A,A,A)') trim(aims_uuid), '-4', s3
            else if (i_entry==5) then
               write(aims_uuid,'(A,A,A,A)') trim(aims_uuid), '-A', s3, '-'
            else    
               write(aims_uuid,'(A,A)') trim(aims_uuid), s4
            end if
         end do
      end if

      if (use_mpi) then
        call MPI_BCAST(aims_uuid, len(aims_uuid), MPI_CHARACTER, 0, mpi_comm_global, error)
      end if

   end subroutine gen_aims_uuid

   subroutine init_random_seed()

      implicit none

      integer (KIND=4), allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8)
      integer (kind=8) :: t

      call random_seed(size = n)
      allocate(seed(n))

      ! First try if the OS provides a random number generator
      ! The following line has been tested on a Microsoft 10 naitive 
      ! environment. On a Windows maschine it cannot find the file
      ! urandom and will use the system clock as a seed instead.
      open(unit=newunit(un), file="/dev/urandom", access="stream", &
           form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
         read(un) seed
         close(un)
      else

         ! Using the current time to generate a seed

         call system_clock(t)
         ! In case that the system_clock is 0 and therefore
         ! not a helpful choice to generate a seed

         if (t == 0) then
            call date_and_time(values=dt)
            t = (dt(1) - 1970) * 365 * 24 * 60 * 60 * 1000 &
                 + dt(2) * 31 * 24 * 60 * 60 * 1000 &
                 + dt(3) * 24 * 60 * 60 * 1000 &
                 + dt(5) * 60 * 60 * 1000 &
                 + dt(6) * 60 * 1000 + dt(7) * 1000 &
                 + dt(8)
         end if

         ! Writting the array with seeds
         do i = 1, n
            seed(i) = lcg(t)
         end do
      end if

       ! Finally setting the random seed
      call random_seed(put=seed)

      end subroutine init_random_seed

end module generate_aims_uuid
