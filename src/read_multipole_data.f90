!       Subroutine read_multipole_data reads the multipole elements
!       from geometry.in. it is called in read_geo.f

      subroutine read_multipole_data &
       (i_multipole, multipole_order, multipole_data, num_of_elem)

      use constants
      use mpi_tasks
      use localorb_io,only:use_unit

      implicit none

!     input
      integer i_multipole
      integer multipole_order
      integer num_of_elem   ! number of multipole tensor elements

!     output
      real*8 multipole_data(num_of_elem)

!     local scope
      character*20 desc_str ! contains subtags
      integer i_index       ! temporary index
      integer i_code        ! read status flag
      logical flag_eod      ! end of data for multipoles
      logical flag_eof      ! end of file
      logical got_data      ! required multipole data obtained?

!     init
      flag_eod =.false.
      flag_eof =.false.
      got_data=.false.

!     multipole order supported?
      if (multipole_order.gt.1) then
        if (myid.eq.0) then
          write(use_unit,*) &
            "FATAL ERROR: embedding multipoles of l>1 is not supported"
        endif
        stop
      endif

!     read multipole data elements (l>0)
      read (8,*,iostat=i_code) desc_str

      if (i_code.ne.0) then
         if (myid.eq.0) then
            write(use_unit,*) "FATAL ERROR: l>0 but no multipole data"
         end if
         stop
      end if

      do while (.not. (flag_eod .or. flag_eof) )
        if (desc_str(1:1).eq."#") then
          continue
        elseif (desc_str.eq."data") then
          backspace(8)
          read(8,*,iostat=i_code) desc_str,  (multipole_data(i_index), &
            i_index=1,num_of_elem,1)
          got_data=.true.
        else
          backspace(8)
          flag_eod =.true.
        endif
        read (8,*,iostat=i_code) desc_str

        if (i_code.ne.0) then
          flag_eod = .true.
          flag_eof = .true.
        end if
      enddo

      backspace(8)

      if(.not. got_data) then
         if (myid.eq.0) then
            write(use_unit,*) "FATAL ERROR: l>0 but no multipole data"
         end if
         stop
      end if

!      write(use_unit,*) "DEBUG: (id,order,data): ", i_multipole,
!     + multipole_order, multipole_data(1:num_of_elem)

!     transform multipole coordinates into atomic units
      do i_index = 1,num_of_elem,1
            multipole_data(i_index) = &
            multipole_data(i_index)/bohr
      enddo
      end

