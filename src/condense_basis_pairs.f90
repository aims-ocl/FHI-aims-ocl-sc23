!****s* FHI-aims/condense_basis_pairs
!  NAME
!    condense_basis_pairs
!  SYNOPSIS

      subroutine condense_basis_pairs()

!  PURPOSE
!   Subroutine "condense_basis_pairs" will find from all the regular
!   basis pairs (n_basis*n_basis) the relevant ones (i.e., those have
!   a non-zero overlap).
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none 
!  OUTPUT
!    none 
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  SOURCE



!   Variables:
!     n_basis_pairs            :
!           number of the nontrival basis pairs
!     n_nghr_basis (n_basis)       :
!           this array tells, for each basis function
!           i_basis=1, ..., n_basis,  how many neighboring basis are relevant.
!     basis_nghr(n_basis, n_basis) :
!           this matrix contains the information of the mapping between the
!           "condensed" basis pairs and the original basis pairs:
!           basis_nghr(i_basis, j_basis) -> i_basis_pair

      use dimensions
      use basis
      use timing ! added for atom bsse      
      use geometry
      use prodbas
      use mpi_tasks
      use localorb_io, only: use_unit
      use runtime_choices, only: calculate_atom_bsse

      implicit none

!    local variables
      real*8   dist
      real*8   radius_sum
      integer  basis_per_atom(n_atoms)

!    counter
      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_atom_1
      integer :: i_atom_2
      integer :: i_fn_1
      integer :: i_fn_2

!    begin to work

!   get maximum number of basis per atom
      basis_per_atom (:) = 0
      do i_basis_1 = 1, n_basis

        i_atom_1 = basis_atom(i_basis_1)
        do i_atom_2 = 1, n_atoms

          if(i_atom_2 .eq. i_atom_1) then
            basis_per_atom(i_atom_2) = &
             basis_per_atom(i_atom_2) + 1
          endif

        enddo
      enddo

      n_max_basis_atom = 0
      do i_atom_2 = 1, n_atoms
         n_max_basis_atom = max(n_max_basis_atom, &
            basis_per_atom(i_atom_2))
      enddo
!      write(use_unit,*)
!      write(use_unit,*) " | n_max_basis_per_atom", n_max_basis_atom

!   get atoms pairs
      do i_atom_1 = 1, n_atoms, 1
        do i_atom_2 = 1, i_atom_1, 1
          atom_pairs(i_atom_2, i_atom_1) = &
             i_atom_1*(i_atom_1-1)/2 + i_atom_2
        enddo
      enddo

      do i_atom_1 = 1, n_atoms, 1
        do i_atom_2 =  i_atom_1+1, n_atoms, 1
          atom_pairs(i_atom_2, i_atom_1) = &
             i_atom_2*(i_atom_2-1)/2 + i_atom_1
        enddo
      enddo

!   get compact basis pairs
      n_basis_pairs = 0
      n_nghr_basis(:) = 0
      basis_nghr(:,:) = 0

      do i_basis_1 = 1, n_basis, 1

         i_atom_1 = basis_atom(i_basis_1)
         i_fn_1 = basis_fn(i_basis_1)

         do i_basis_2 = 1, n_basis, 1

           i_atom_2 = basis_atom(i_basis_2)
           i_fn_2 = basis_fn(i_basis_2)

           dist = (coords(1,i_atom_1) - coords(1,i_atom_2)) * &
                  (coords(1,i_atom_1) - coords(1,i_atom_2)) + &
                  (coords(2,i_atom_1) - coords(2,i_atom_2)) * &
                  (coords(2,i_atom_1) - coords(2,i_atom_2)) + &
                  (coords(3,i_atom_1) - coords(3,i_atom_2)) * &
                  (coords(3,i_atom_1) - coords(3,i_atom_2))

           dist = sqrt(dist)

           radius_sum = outer_radius(i_fn_1) + outer_radius(i_fn_2)

           if(radius_sum .ge. dist) then

             if(i_basis_2 .le. i_basis_1) then
               n_basis_pairs = n_basis_pairs + 1
               n_nghr_basis(i_basis_1) = n_nghr_basis(i_basis_1) + 1
               basis_nghr(i_basis_2, i_basis_1) = n_basis_pairs
               basis_nghr(i_basis_1, i_basis_2) = n_basis_pairs
             endif

           endif

!       loop over i_basis_2
        enddo
!       loop over i_basis_2
      enddo

!       do i_basis_1 = 1, n_basis, 1
!         do i_basis_2 = i_basis_1+1, n_basis, 1
!              basis_nghr(i_basis_2, i_basis_1) = &
!                 basis_nghr(i_basis_1, i_basis_2)
!         enddo
!       enddo


       if(myid .eq. 0) then
         write(use_unit,*)
         write(use_unit,'(A,I10,A, I10)') &
            "  Basis pair condensation :", n_basis*(n_basis+1)/2, &
            "  --> ", n_basis_pairs
       endif
       
!  ATOM BSSE: create a copy of basis_nghr and n_nghr_basis that 
!  will remain as the reference for the subsequent bsse geometries
     if (calculate_atom_bsse) then
       do i_basis_1 = 1, n_basis, 1
         n_nghr_basis_original(i_basis_1) = n_nghr_basis(i_basis_1)
         do i_basis_2 = 1, n_basis, 1
              basis_nghr_original(i_basis_2, i_basis_1) = &
                 basis_nghr(i_basis_2, i_basis_1)
         enddo
       enddo
     endif    


!      deallocate (radial_ovlp)

      end subroutine condense_basis_pairs
!*****
!-----------------------------------------------------------------------------
!****s* FHI-aims/condense_compute_pairs
!  NAME
!   condense_compute_pairs
!  SYNOPSIS

      subroutine condense_compute_pairs &
         ( i_atom, n_compute, i_basis, &
           n_compute_current, i_basis_current, &
           n_compact, compact_pairs, &
           n_compact_current, compact_pairs_current &
          )

!  PURPOSE
!   Now for a sub basis pool 1,2, ..., n_compute (typically n_compute < n_basis),
!   perform roughly the same "condensation" procedure as above. Now for a give
!   basis function i_compute, find the relevant basis functions  from the pool.
!   The information is contained in array "n_compact" and "compact_pairs".
!
!  USES

      use dimensions
      use basis
      use geometry
      use prodbas

!  ARGUMENTS
      integer :: i_atom
      integer :: n_compute
      integer :: i_basis(n_compute)
      integer :: n_compute_current
      integer :: i_basis_current(n_compute)
      integer :: n_compact(n_basis)
      integer :: compact_pairs(n_basis,n_basis)
      integer :: n_compact_current(n_basis)
      integer :: compact_pairs_current(n_max_basis_atom,n_basis)

!  INPUTS
!  o i_atom -- the current atom under consideration
!  o n_compute -- the number of non-zero basis function
!  o i_basis -- array, specifies these non-zero basis function
!  OUTPUT
!   o n_compute_current -- the non-zero basis funciton belonging to the current atom
!   o i_basis_current -- specifies the non-zero basis funciton belonging to the current atom
!   o n_compact -- array, for each basis in the batch, the number of basis
!     functions in the same batch which has zero overlap with the one under
!     consideration
!   o compact_pairs -- two-dimentional array, for each basis in the batch,
!     list the basis functions in the same batch which has zero overlap with
!     the one under consideration
!   o n_compact_current -- array, for each basis in the batch, the number of basis
!     functions in the same batch that belong to the current atom and has
!     zero overlap with the one under consideration
!   o compact_pairs_current -- two-dimentional array, for each basis in the batch,
!     list the basis functions in the same batch which has zero overlap with
!     the one under consideration and belong to the curren atom
!  SOURCE

!    counter
      integer :: i_compute
      integer :: i_compute_1
      integer :: i_basis_1
      integer :: i_basis_2
      integer :: i_atom_1
      integer :: i_atom_2
      integer :: i_index

!  begin to work

      n_compact = 0
      n_compact_current = 0
      compact_pairs(:,:) = 0
      compact_pairs_current(:,:) = 0

      do i_compute = 1, n_compute, 1

         i_basis_1 = i_basis(i_compute)
         i_atom_1 = basis_atom(i_basis_1)
         if(i_atom_1.eq.i_atom) then
             n_compute_current = n_compute_current + 1
             i_basis_current(n_compute_current)=i_compute
         endif

         do i_compute_1 = 1, i_compute, 1

           i_basis_2 = i_basis(i_compute_1)
           i_index = basis_nghr(i_basis_2,i_basis_1)

           if(i_index .ge. 0) then
             n_compact(i_compute) = n_compact(i_compute) + 1
             compact_pairs(n_compact(i_compute),i_compute) = &
               i_compute_1
           endif

         enddo

         do i_compute_1 = 1, n_compute, 1

           i_basis_2 = i_basis(i_compute_1)
           i_atom_2 = basis_atom(i_basis_2)
           i_index = basis_nghr(i_basis_2,i_basis_1)

           if(i_atom_2.eq.i_atom .and. i_index .gt.0) then

               n_compact_current(i_compute) = &
                 n_compact_current(i_compute) + 1
               compact_pairs_current(n_compact_current(i_compute), &
                   i_compute) = i_compute_1

!              write(use_unit,*) i_compute_1, i_compute,
!     +          compact_pairs_current(n_compact_current(i_compute),
!     +               i_compute)

           endif

         enddo
      enddo

      end subroutine condense_compute_pairs
!******
!****s* FHI-aims/bsse_nghr_mapping
!  NAME
!    bsse_nghr_mapping
!  SYNOPSIS

      subroutine bsse_nghr_mapping()

!  PURPOSE
!   Subroutine "bsse_nghr_mapping" will, for the new bsse geometry with 
!   ghost atoms, find from all the regular
!   basis pairs (n_basis*n_basis) the relevant ones (i.e., those have
!   a non-zero overlap).
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none 
!  OUTPUT
!    none 
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  SOURCE



!
      use dimensions
      use timing ! for atom_bsse
      use basis
      use geometry
      use prodbas
      use mpi_tasks
      use localorb_io, only: use_unit

      implicit none

!    local variables

      integer ::      bsse_n_basis_pairs 
      integer, dimension(:), allocatable :: bsse_n_nghr_basis
      integer, dimension(:,:), allocatable :: bsse_basis_nghr    

!    counter
      integer :: i_basis_1
      integer :: i_basis_2


!    begin to work
!    allocate new matrices for bsse geometries
      allocate (bsse_n_nghr_basis(n_basis))
      allocate (bsse_basis_nghr(n_basis,n_basis))

       if (myid.eq.0) then
             write(use_unit,'(1X,A)') &
           "ATOM BSSE : creating basis_nghr for the new structure  "
        endif
!   get compact basis pairs
      bsse_n_basis_pairs = 0
      bsse_n_nghr_basis(:) = 0
      bsse_basis_nghr(:,:) = 0
      
!   assign bsse_n_nghr_basis      
      
      do i_basis_1 = 1, n_basis, 1
        bsse_n_nghr_basis(i_basis_1) = n_nghr_basis_original(basis_mapping(i_basis_1))
      enddo

!   assign bsse_basis_nghr

      do i_basis_1 = 1, n_basis, 1
         do i_basis_2 = 1, i_basis_1, 1
              bsse_basis_nghr(i_basis_2, i_basis_1) = &
                        basis_nghr_original(basis_mapping(i_basis_2), basis_mapping(i_basis_1))
         enddo
     enddo
     ! the lower triangle
      do i_basis_1 = 1, n_basis, 1
         do i_basis_2 = i_basis_1+1, n_basis, 1
              bsse_basis_nghr(i_basis_2, i_basis_1) = &
                 bsse_basis_nghr(i_basis_1, i_basis_2)
         enddo
      enddo
      
      
!    swap old with new indices:  only from second atom onwards
     ! if (current_atom_num_bsse>1) then
       do i_basis_1 = 1, n_basis, 1
         n_nghr_basis(i_basis_1) = bsse_n_nghr_basis(i_basis_1)
         do i_basis_2 = 1, n_basis, 1
              basis_nghr(i_basis_2, i_basis_1) = &
                 bsse_basis_nghr(i_basis_2, i_basis_1)
         enddo
       enddo
     ! endif


      end subroutine bsse_nghr_mapping
