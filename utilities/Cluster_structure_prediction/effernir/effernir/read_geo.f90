!  Subroutine read_geo provides cluster geometry, and calculates
!  related information

      subroutine read_geo(filename)

      use cluster

      implicit none

      ! fundamental constants

      ! imported variables
      
      ! input
      character*30, intent(in) :: filename

      include "constants.f90"

      ! local variables

      ! eof : end-of-file marker
      ! i_code : iostatus flag
      ! desc_str : line descriptor
      ! species_temp: temporary placeholder for species_name
      ! found : flag indicating whether each species is actually defined.
      
      logical :: eof
      integer :: i_code
      character*20 :: desc_str
      character*20 :: species_temp
      logical :: found

!  counters

      integer :: i_species
      integer :: i_atom

!  begin work

      write(6, *) 
      write(6,'(A)') &
      "------------------------------------------------------------"
      write(6,'(10X,A)') "Reading geometry description geometry.in."
      write(6,'(A)') &
      "------------------------------------------------------------"

      if (.not.associated(initial_coords%coords)) then
         write(*,*) "Memory for arrays not allocated."
         write(*,*) "Aborting."
         stop
      end if
!  initialize

      eof = .false.

      i_atom = 0
    
!  open input file

      open (8, FILE = filename)
!  read first line     

      read (8, *, iostat=i_code) desc_str
      if (i_code .ne. 0) then

        write(*,*) "Invalid file geometry.in, or file not found."
        write(*,*) "Aborting."
        stop

      end if

      do while (.not.eof)

         select case(desc_str) 
         case ("#")
            
            continue
            
         case('atom') 

            backspace(8)
            
            i_atom = i_atom + 1
            
            if (i_atom .gt. size(initial_coords%coords)) then
               write(*,*) "Not enough memory for vector-arrays allocated."
               write(*,*) "Aborting."
               stop
            end if
            
            read(8,*) desc_str, initial_coords%coords(i_atom)%x, &
                 initial_coords%coords(i_atom)%y, initial_coords%coords(i_atom)%z, &
                 species_temp
            
            !  check whether we know this species
             
            found = .false.
            
            do i_species = 1, n_species, 1
               
               if (species_temp .eq. species_name(i_species)) then
                  species(i_atom) = i_species
                  found = .true.
               end if
               
            end do
            
            if (.not.found) then
               
               write(6,'(1X,A,A,A,A)') &
                    "* Species ", species_temp, ", listed in ", &
                    "geometry.in, is not described in input file control.in."
               write(6,*) "* List of ", n_species, " species: "
               do i_species = 1, n_species, 1
                  write(*,*) "* ",i_species,":",species_name(i_species)
               enddo
               write(6,*) "* Aborting."
               stop
               
            end if

         case('initial_charge') 

            backspace(8)
            read(8,*) desc_str, ini_charge_per_atom(i_atom)

         case('initial_moment') 

            backspace(8)
            read(8,*) desc_str, ini_spin_per_atom(i_atom) 

         end select
         
!  next line

        read(8, *, iostat=i_code) desc_str

        if (i_code .ne. 0) then
          eof = .true.
        end if

      end do

!  close geometry.in

      write (*,*) "done."

      close(8)

    end subroutine read_geo
