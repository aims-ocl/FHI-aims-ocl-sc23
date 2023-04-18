!****s* FHI-aims/integrate_dipmom_pairstates
!  NAME
!   integrate_dipmom_pairstates
!  SYNOPSIS

      subroutine integrate_dipmom_pairstates &
      (n_occ, n_unocc, n_first_unocc, coord_of_center, &
       basis_l_max, KS_eigenvector, dipole_mom)

!  PURPOSE
!  This subroutine is intended to calculate the dipole moments for
!  every pairs of occupied-virtual product wave function.
!
!  USES

      use dimensions
      use runtime_choices
      use grids
      use geometry
      use basis
      use mpi_utilities
      use synchronize_mpi
      use prodbas
      use constants
      use localorb_io, only: use_unit

      implicit none

!  ARGUMENTS
      integer n_occ
      integer n_unocc
      integer n_first_unocc(n_spin)
      integer basis_l_max(n_species)
      real*8 coord_of_center(3)
      real*8 KS_eigenvector(n_basis, n_states, n_spin)
      real*8 dipole_mom(n_occ,n_unocc,n_spin,3)

!  INPUTS
!  o n_occ -- number of occupied state
!  o n_unocc -- number of unoccupied states
!  o n_first_unocc -- first unoccupied state
!  o basis_l_max -- integer array, the maximal angular momentum number of the basis 
!           functions for each species        
!  o coord_of_center -- coordinate of center of the molecule 
!  o KS_eigenvector -- eigenvector of ground-state DFT/HF calculation
!
!  OUTPUTS
!  o dipole_mom -- dipole moment of each pair product 
!
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

!  local variables

      integer :: l_ylm_max
      integer, dimension(:,:), allocatable :: index_lm
      real*8, dimension(:,:,:), allocatable :: ylm_tab
      real*8, dimension(:,:), allocatable :: tmp_dipole_mom
      real*8, dimension(:,:,:), allocatable :: dipole_mom_basis
      real*8, dimension(:,:,:), allocatable :: wave_times_r(:,:,:)

      real*8 coord_current(3, n_max_angular)
      real*8 coord_wrt_center(3, n_max_angular)
      real*8 dist_tab(n_atoms, n_max_angular)
      real*8 i_r(n_atoms, n_max_angular)
      real*8 dir_tab(3,n_atoms, n_max_angular)
      real*8 trigonom_tab(4,n_atoms, n_max_angular)

      real*8 wave(n_basis, n_max_angular)
      real*8 radial_wave(n_basis)

      integer :: n_compute
      integer :: n_compute_onsite
      integer :: n_max_onsite_basis
      integer :: i_basis(n_basis)
      integer :: i_basis_onsite(n_basis)

!     optimal accounting for matrix multiplications: only use points with nonzero components
      integer :: n_points

!     and condensed version of partition_tabs on angular grids
      real*8 :: partition_tab(n_max_angular)
      real*8 :: partition_tab_2atoms &
                (n_atoms, n_max_angular)

      real time_stamp
      real time_end

!     partition_type_temp makes sure that we use the right partition function for integrals
      integer :: partition_type_temp = 1


!  counters

      integer i_basis_1
      integer i_basis_2
      integer i_atom
      integer i_atom_1
      integer i_radial
      integer i_angular
      integer i_coord
      integer i_spin
      integer i_l
      integer i_m
      integer i_index

      integer i_species

      integer i_compute
      integer i_compute_1

      integer i_task
      integer i_basis_index

!  begin work


      if(myid.eq.0) then
       write(use_unit,'(2X,A,A)') &
        "Integrating the dipole moment assoicated with products of", &
        " KS states ..."
      endif

!     begin with general allocations
        l_ylm_max = l_wave_max

      allocate( ylm_tab( (l_ylm_max+1)**2, n_atoms, &
           n_max_angular) )
      allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max) )

!     initialize index_lm

      i_index = 0
      do i_l = 0, l_wave_max, 1
        do i_m = -i_l, i_l
          i_index = i_index+1
          index_lm(i_m,i_l) = i_index
        enddo
      enddo

      n_max_onsite_basis = 0
      do i_atom = 1, n_atoms, 1 
        n_compute = 0 
        do i_basis_1 = 1, n_basis, 1
          if(basis_atom(i_basis_1).eq.i_atom) then
            n_compute = n_compute + 1
          endif
        enddo
        if (n_max_onsite_basis .lt. n_compute) then
           n_max_onsite_basis = n_compute
        endif
      enddo 
!      write(use_unit,*) "n_max_onsite_basis", n_max_onsite_basis

      allocate( tmp_dipole_mom (n_basis,n_max_onsite_basis))
      allocate( dipole_mom_basis (n_basis,n_basis,3))
      allocate( wave_times_r (n_max_onsite_basis,n_max_angular,3))

      dipole_mom_basis(:,:,:) = 0.d0


      i_task = myid + 1
!     perform partitioned integration, atom by atom, and point by point
!     This will be the outermost loop, to save evaluations of the potential.
!     and the Y_lm functions
      do i_atom = 1, n_atoms, 1

!test
!       if(myid.eq.0) then
!        write(use_unit,*) " | i_atom: ", i_atom
!       endif
!test end

         do i_radial = 1, n_radial(species(i_atom)), 1

!           if (myid.eq.radial_task_list(i_radial,i_atom)) then


!test
!        if(myid.eq.0) then
!          write(use_unit,*) "  i_radial: ", i_radial
!        endif
!test end

           n_compute = 0
           i_basis = 0

           do i_angular = 1, n_angular(i_radial, species(i_atom)), 1


!test
!          if (i_radial.eq.140) then
!            write(use_unit,*) "  i_angular: ", i_angular
!          end if
!test end

!     get current integration point coordinate
              do i_coord = 1, 3, 1
                 coord_current(i_coord,i_angular) = &
                      coords(i_coord,i_atom ) + &
                      r_angular(i_coord, i_angular, i_radial, &
                      species(i_atom)) * &
                      r_radial(i_radial, species(i_atom))

                 coord_wrt_center(i_coord,i_angular) = &
                    coord_current(i_coord,i_angular) - &                    
                    coord_of_center(i_coord) 
              enddo

!     compute atom-centered coordinates of current integration point,
!     as viewed from all atoms
              call tab_atom_centered_coords &
                   ( coord_current(1,i_angular), &
                   dist_tab(1,i_angular), i_r(1,i_angular), &
                   dir_tab(1,1,i_angular) &
                   )


!test
!     if (i_radial.eq.140) then
!     write(use_unit,*) "here 1"
!     end if
!test end

!     evaluate the partition function for two atoms

                  call evaluate_partition_tab_2atoms &
                    ( i_atom, dist_tab(1:n_atoms,i_angular), &
                      i_r(1:n_atoms,i_angular), &
                      w_radial( i_radial, species (i_atom)), &
                      w_angular( i_angular, i_radial, species (i_atom)), &
                      partition_tab_2atoms(1, i_angular), &
                      partition_type_temp )


              partition_tab(i_angular) = 0.d0
              do i_atom_1 =1, n_atoms, 1
               partition_tab(i_angular) = &
                   max( partition_tab(i_angular), &
                        partition_tab_2atoms(i_atom_1, i_angular) )
              enddo

              do i_atom_1 = 1, n_atoms
                 if(dist_tab(i_atom_1,i_angular) .lt. 1.e-15) then
                   partition_tab(i_angular) = 0.d0
                   exit
                 endif
              enddo

!     determine which basis functions are relevant at current integration point,
!     and tabulate their indices
              if (partition_tab(i_angular).gt.0.d0) &
                   then

                 call prune_basis_v1(dist_tab(1,i_angular), n_compute, &
                      i_basis)

              end if
           enddo

           i_basis_onsite(:) = 0
           n_compute_onsite = 0
           do i_compute = 1, n_compute, 1
            if (i_atom .eq. basis_atom(i_basis(i_compute))) then
               n_compute_onsite = n_compute_onsite + 1   
               i_basis_onsite(n_compute_onsite) = i_compute
            endif
           enddo

           if (n_compute.gt.0) then

              n_points = 0
              do i_angular = 1, n_angular(i_radial, species(i_atom)), 1

                 if (partition_tab(i_angular).gt.0.d0) &
                      then
!     execute only if partition_tab.gt.0 here, i.e. if the integration point
!     makes sense
                    n_points = n_points + 1
!                    partition_tab_2atoms(:, n_points) = &
!                         partition_tab_2atoms(:, i_angular)
!                    partition(:, n_points) = &
!                         partition_tab(:, i_angular)
!     compute trigonometric functions of spherical coordinate angles
!     of current integration point, viewed from all atoms
                    call tab_trigonom &
                         ( dir_tab(1,1,i_angular), &
                         trigonom_tab(1,1,i_angular) &
                         )

!     tabulate distance and Ylm's w.r.t. other atoms
                    call tab_wave_ylm &
                         ( trigonom_tab(1,1,i_angular), basis_l_max, &
                         l_ylm_max, &
                         ylm_tab(1,1,i_angular) )

!           tabulate total wave function value for each basis function
                    call evaluate_waves_v0 &
                         (i_r(1,i_angular), l_ylm_max, &
                         ylm_tab(1,1,i_angular), &
                         dist_tab(1,i_angular), index_lm, n_compute, &
                         i_basis, radial_wave(1), &
                         wave(1,n_points))

                     do i_coord = 1, 3, 1
                      do i_compute_1 = 1, n_compute_onsite, 1

                         i_compute = i_basis_onsite(i_compute_1)
                         wave_times_r(i_compute_1,n_points,i_coord) = &
                          wave(i_compute,n_points)  * &
                          coord_wrt_center(i_coord, i_angular) 

                      enddo
                     enddo

                     do i_compute = 1, n_compute, 1
                       wave(i_compute,n_points)= wave(i_compute,n_points)  &
                        * partition_tab_2atoms(basis_atom(i_basis(i_compute)), &
                                                i_angular)
                     enddo

!     end if (partition_tab.gt.0)
             end if

!     end angular integration loop
          enddo
!     add the contribution from the current shell

          do i_coord=1, 3, 1
            tmp_dipole_mom(:,:) = 0.d0
           call dgemm('N','T',n_compute,n_compute_onsite,n_points,1.d0, &
                      wave(1,1), n_basis, &
                      wave_times_r(1,1,i_coord), n_max_onsite_basis, 0.d0, &
                      tmp_dipole_mom(1,1), n_basis )
!           tmp_dipole_mom(:,:) = 0.d0
!           do i_point = 1, n_points, 1
!             do i_compute_1 = 1, n_compute_onsite, 1
!               do i_compute = 1, n_compute, 1
!                   tmp_dipole_mom(i_compute, i_compute_1) = &
!                     tmp_dipole_mom(i_compute, i_compute_1) + &
!                     wave(i_compute,i_point) * &
!                     wave_times_r(i_compute_1,i_point,i_coord)
!                enddo
!              enddo 
!            enddo
!            write(use_unit,*) i_radial, n_points, tmp_dipole_mom(1,1), tmp_dipole_mom(1,2)
     
           do i_compute = 1, n_compute_onsite,  1
             i_basis_1 = i_basis(i_basis_onsite(i_compute))
             do i_compute_1 = 1, n_compute,  1
               i_basis_2 = i_basis(i_compute_1)
               dipole_mom_basis(i_basis_2,i_basis_1, i_coord) = &
                   dipole_mom_basis(i_basis_2, i_basis_1, i_coord) + &
                   tmp_dipole_mom(i_compute_1,i_compute)
             enddo
           enddo

          enddo
!     end if (n_compute.gt.0) then
       end if

!test
!            if (i_radial.eq.140) then
!            write(use_unit,*) "here 3"
!            end if
!test end

!       end radial integration loop
      enddo

!     end integration loop over atoms
      enddo

      do i_basis_1 = 1, n_basis, 1
        do i_basis_2 = 1, i_basis_1-1, 1
         if(basis_atom(i_basis_2) .ne. basis_atom(i_basis_1)) then
          dipole_mom_basis(i_basis_2,i_basis_1,:) = &
            dipole_mom_basis(i_basis_2,i_basis_1,:) + &
            dipole_mom_basis(i_basis_1,i_basis_2,:) 

           dipole_mom_basis(i_basis_1,i_basis_2,:) = &
             dipole_mom_basis(i_basis_2,i_basis_1,:) 
          endif
         enddo
       enddo
         
!         do i_basis_1 =1, n_basis, 1
!          do i_basis_2 = 1, i_basis_1, 1
!             write(use_unit, '(4X,4I5,2X,3F16.10)') &
!                 i_basis_1, i_basis_2, basis_l(i_basis_1), &
!                 basis_l(i_basis_2), &
!                 dipole_mom_basis(i_basis_2,i_basis_1,:)
!         enddo
!        enddo

      deallocate(tmp_dipole_mom)
      allocate(tmp_dipole_mom(n_basis,n_occ))

      do i_coord = 1, 3, 1
        do  i_spin = 1, n_spin, 1
          call dgemm('N', 'N', n_basis, n_occ, n_basis, 1.d0, &
                     dipole_mom_basis(1,1,i_coord), n_basis, &
                     KS_eigenvector(1,1,i_spin), n_basis, 0.d0, &      
                     tmp_dipole_mom(1,1), n_basis)
                     
          call dgemm('T', 'N', n_occ, n_unocc, n_basis, 1.d0, &
                     tmp_dipole_mom(1,1), n_basis, &
                     KS_eigenvector(1,n_states-n_unocc+1,i_spin), n_basis, &
                     0.d0, dipole_mom(1,1,i_spin,i_coord), n_occ)

!         tmp_dipole_mom(:,:) = 0.d0
!         do i_basis_1 = 1, n_occ, 1
!           do i_basis_2 = 1, n_basis, 1
!              do i_compute = 1, n_basis, 1
!                tmp_dipole_mom(i_basis_2, i_basis_1) = &
!                    tmp_dipole_mom(i_basis_2, i_basis_1) +  &
!                    dipole_mom_basis(i_basis_2, i_compute, i_coord) * &
!                    KS_eigenvector(i_compute, i_basis_1, i_spin)
!              enddo
!           enddo
!         enddo
!        
!          write(use_unit,*)"tmp_dipole_mom", tmp_dipole_mom(:,:)
!
!         dipole_mom(:,:,i_spin,i_coord) = 0.d0
!         do i_basis_1 = 1, n_unocc, 1
!           do i_basis_2 = 1, n_occ, 1
!              do i_compute = 1, n_basis, 1
!                 dipole_mom(i_basis_2, i_basis_1,i_spin,i_coord) = &
!                     dipole_mom(i_basis_2, i_basis_1,i_spin,i_coord) + &
!                    tmp_dipole_mom(i_compute, i_basis_2) *  &
!                    KS_eigenvector(i_compute, i_basis_1+n_states-n_unocc, i_spin)
!              enddo
!           enddo
!         enddo
        
!         do i_basis_1 =1, n_unocc, 1
!          do i_basis_2 = 1, n_occ, 1
!             write(use_unit, '(4X,4I5,2X,3F16.10)') &
!                 i_coord, i_spin, i_basis_1+n_states-n_unocc, i_basis_2, &
!                 dipole_mom(i_basis_2,i_basis_1,i_spin,i_coord)
!         enddo
!        enddo
!
        enddo
      enddo

      if (allocated(ylm_tab)) then
        deallocate(ylm_tab)
      end if
      if (allocated(index_lm)) then
        deallocate(index_lm)
      end if
      if (allocated(tmp_dipole_mom)) then
        deallocate(tmp_dipole_mom)
      end if
      if (allocated(dipole_mom_basis)) then
        deallocate(dipole_mom_basis)
      end if
      if (allocated(wave_times_r)) then
        deallocate(wave_times_r)
      end if
      return
      end subroutine integrate_dipmom_pairstates

!----------------------------------------------------------------------
!******
