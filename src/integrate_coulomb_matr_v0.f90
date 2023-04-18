!****s* FHI-aims/integrate_coulomb_matr_v0
!  NAME
!   integrate_coulomb_matr_v0
!  SYNOPSIS

      subroutine integrate_coulomb_matr_v0 &
      (basis_l_max, coulomb_matr)

!  PURPOSE
!  this subroutine is intended to calculate the coulomb interaction matrix
!  elements between every two basis functions, revised from the subroutine
!  "integrate_real_hamiltonian_matrix".
!  <i| v(r-r') |j>
!
!  USES

      use dimensions
      use runtime_choices
      use grids
      use geometry
      use basis
      use prodbas
      use mpi_tasks
      use synchronize_mpi
      use constants
      use localorb_io, only: use_unit

      implicit none

!  ARGUMENTS

      integer basis_l_max (n_species)
      real*8 coulomb_matr( n_basbas,n_loc_prodbas )

!  INPUTS
!  o basis_l_max -- integer array, the maximal angular momentum number of the basis
!           functions for each species
!  OUTPUTS
!  o coulomb_matr -- real array, the Coulomb interaction matrix within the auxiliary
!            basis         
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
!    v_times_wave:
!    ovlp_prodbas: overlap matrix for the product basis functions

      integer :: l_ylm_max
      integer, dimension(:,:), allocatable :: index_lm
      real*8, dimension(:,:,:), allocatable :: ylm_tab
      real*8, dimension(:,:,:), allocatable :: &
                      v_times_radialwaves_spl
      real*8, dimension(:,:), allocatable :: v_times_waves
      real*8, dimension(:,:), allocatable :: ovlp_prodbas
      real*8, dimension(:,:), allocatable :: wave


      real*8 coord_current(3, n_max_angular)
      real*8 dist_tab(n_atoms, n_max_angular)
      real*8 i_r(n_atoms, n_max_angular)
      real*8 dir_tab(3,n_atoms, n_max_angular)
      real*8 trigonom_tab(4,n_atoms, n_max_angular)


      integer :: n_compute
      integer :: n_compute_current
      integer :: n_prodbas_current
      integer :: i_basis(n_basbas)
      integer :: i_basis_current(n_basbas)
      integer :: i_prodbas(n_loc_prodbas)
      integer :: i_prodbas_current(n_loc_prodbas)
      integer :: basbas_l_max(n_species)

!     optimal accounting for matrix multiplications: only use points with nonzero components
      integer :: n_points

!     and condensed version of partition_tabs on angular grids
      real*8 :: partition_tab(n_max_angular)
      real*8 :: partition(n_max_angular)
      real*8 :: partition_tab_2atoms &
                (n_atoms,n_max_angular)

!      real*8 :: x_matrix(n_basbas,n_basbas)
!      real*8, dimension(:,:,:), allocatable :: gradient_basis_wave

!      real*8, dimension(n_max_angular) :: zora_operator
!      logical, dimension(n_max_angular) :: t_zora
!      real*8, dimension(n_basis,3,n_max_angular) :: zora_vector

      real time_stamp
      real time_end

!     partition_type_temp makes sure that we use the right partition function for integrals
      integer :: partition_type_temp = 1


!  counters

      integer i_basis_1
      integer i_basis_2
      integer i_basbas
      integer i_atom
      integer i_atom_1
      integer i_radial
      integer i_angular
      integer i_grid
      integer i_index, i_l, i_m
      integer i_coord

      integer i_species

      integer i_compute
      integer i_compute_2
      integer i_point

      integer i_task
!  begin work

      if(myid.eq.0) then
       write(use_unit,*)
       write(use_unit,'(2X,A,A)') &
       "Integrating the coulomb interaction matrix for basis functions" &
           , " ..."
       write(use_unit,*)
      endif

!     partition_type_temp makes sure that we use the right partition function for integrals
      if (partition_type.eq.6) then
        partition_type_temp = 6
      else
        ! just use type 1 anyway
        partition_type_temp = 1
      end if

!     begin with general allocations
      l_ylm_max = 2* l_wave_max

      allocate( v_times_radialwaves_spl ( &
                n_max_spline, n_hartree_grid, n_loc_prodbas) )

      allocate( v_times_waves (n_loc_prodbas, n_max_angular) )

      allocate( ylm_tab( (l_ylm_max+1)**2, n_atoms, &
                n_max_angular) )
      allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max) )
      allocate( wave(n_basbas, n_max_angular))
!      allocate (ovlp_prodbas(n_basbas, n_basbas))

!     initialize

      i_index = 0
!      do i_basis_1 = 1, n_basis, 1
!          do i_basis_2 = 1, i_basis_1, 1
!            i_index = i_index+1

!           coulomb_matr(i_index) = 0.
!        enddo
!      enddo

      call integrate_v_times_radialwaves &
            (v_times_radialwaves_spl)


       coulomb_matr=0.d0


!     initialize index_lm

      i_index = 0
      do i_l = 0, l_ylm_max, 1
        do i_m = -i_l, i_l
          i_index = i_index+1
          index_lm(i_m,i_l) = i_index
        enddo
      enddo

      basbas_l_max = 2*basis_l_max

      i_task = myid + 1

!     perform partitioned integration, atom by atom, and point by point
!     This will be the outermost loop, to save evaluations of the potential.
!     and the Y_lm functions
      do i_atom = 1, n_atoms, 1

!test
!        write(use_unit,*) "i_atom: ", i_atom
!test end

         do i_radial = 1, n_radial(species(i_atom)), 1

!           if(myid.eq.radial_task_list(i_radial,i_atom)) then
!test
!          write(use_unit,*) "  i_radial: ", i_radial
!test end

           n_compute = 0
           n_compute_current = 0
           i_basis = 0
           i_basis_current = 0

           n_prodbas_current = 0
           i_prodbas_current = 0


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
              enddo

!     compute atom-centered coordinates of current integration point,
!     as viewed from all atoms
!          if(i_atom.eq.8.and.i_radial.eq.1) then
              call tab_atom_centered_coords &
                   ( coord_current(1,i_angular), &
                   dist_tab(1,i_angular), i_r(1,i_angular), &
                   dir_tab(1,1,i_angular) &
                   )
!              write(use_unit,*) "corrd_current", coord_current(:,i_angular)
!              write(use_unit,*) "corrd_atom", coords(:,:)
!          endif

!              call evaluate_partition_tab
!     +             (i_atom, dist_tab(1:n_atoms,i_angular),
!     +              i_r(1:n_atoms,i_angular),
!     +              w_radial( i_radial, species (i_atom)),
!     +              w_angular( i_angular, i_radial, species (i_atom)),
!     +              partition_tab(i_angular),
!     +              partition_type_temp )

!     evaluate the partition function for three atoms
                  call evaluate_partition_tab_2atoms &
                    ( i_atom, dist_tab(1:n_atoms,i_angular), &
                      i_r(1:n_atoms,i_angular), &
                      w_radial( i_radial, species (i_atom)), &
                      w_angular( i_angular, i_radial, species (i_atom)), &
                      partition_tab_2atoms(1, i_angular), &
                      partition_type_temp )

                   partition_tab(i_angular) = 0.d0
                   do i_atom_1 =1, n_atoms
                      partition_tab(i_angular) = &
                       max( partition_tab(i_angular), &
                            partition_tab_2atoms(i_atom_1,i_angular) )
                   enddo
                   do i_atom_1 =1, n_atoms
                    if( dist_tab(i_atom_1,i_angular) .lt. 1.e-15) then
                       partition_tab(i_angular) = 0.d0
                       exit
                    endif
                   enddo

!test
!     if (i_radial.eq.140) then
!     write(use_unit,*) "here 1"
!     end if
!test end

!     determine which basis functions are relevant at current integration point,
!     and tabulate their indices
!               if (partition_tab(i_angular).gt.0.d0)
!     +             then
!                 call prune_basbas_v1(dist_tab(1,i_angular), n_compute,
!     +                i_basis)
!               end if

              enddo

!     from the relevant basis function set i_basis find those which belong to
!     the atom i.
              n_compute = n_basbas
              do i_compute = 1, n_basbas

                i_basis(i_compute) = i_compute

                if(basbas_atom(i_basis(i_compute)).eq.i_atom) then
                    n_compute_current = n_compute_current + 1
                    i_basis_current(n_compute_current) = &
                         i_compute
                endif
             enddo

             do i_basis_1 = 1, n_loc_prodbas

                i_prodbas(i_basis_1) = i_basis_1

                i_basbas = map_prodbas(i_basis_1, i_task)
                if(i_basbas.gt.0) then
                   if(basbas_atom(i_basbas) .eq. i_atom) then
                       n_prodbas_current = n_prodbas_current + 1
                       i_prodbas_current ( n_prodbas_current) = &
                          i_basis_1
                   endif
                endif
             enddo

!test
!     if (i_radial.eq.140) then
!     write(use_unit,*) "here 2"
!     end if
!test end

            if (n_compute_current.gt.0) then

              n_points = 0
              v_times_waves =0.d0

              do i_angular = 1, n_angular(i_radial, species(i_atom)), 1

                 if (partition_tab(i_angular).gt.0.d0) &
                      then
!     execute only if partition_tab.gt.0 here, i.e. if the integration point
!     makes sense
                    n_points = n_points + 1
                    partition_tab_2atoms(:,n_points) = &
                         partition_tab_2atoms(:,i_angular)
!     compute trigonometric functions of spherical coordinate angles
!     of current integration point, viewed from all atoms
                    call tab_trigonom &
                         ( dir_tab(1,1,i_angular), &
                         trigonom_tab(1,1,i_angular) &
                         )

!     tabulate distance and Ylm's w.r.t. other atoms
                    call tab_wave_ylm &
                         ( trigonom_tab(1,1,i_angular), basbas_l_max, &
                         l_ylm_max, &
                         ylm_tab(1,1,i_angular) )

!          if(i_atom.eq.8.and.i_radial.eq.1) then
!                    write(use_unit,*)"dist", n_atoms, dist_tab(:,i_angular)
!                    write(use_unit,*)"partition", partition_tab(i_angular)
!           tabulate total wave function value for each basis function
                    call evaluate_prod_waves &
                         (i_r(1,i_angular), l_ylm_max, &
                         ylm_tab(1,1,i_angular), &
                         dist_tab(1,i_angular), index_lm, &
                         n_compute, &
                         i_basis, &
                         wave(1,n_points))


!     evaluate v * psi at this particular point

                    call evaluate_v_times_waves &
                         (l_ylm_max, ylm_tab(1,1,i_angular), &
                         dist_tab(1,i_angular), index_lm, &
                         n_loc_prodbas, i_prodbas, &
                         v_times_radialwaves_spl(1,1,1), &
                         v_times_waves(1:n_loc_prodbas,n_points))

!           if(i_radial .eq. 5) then 
!            write(use_unit,*) "i_r", i_r(:,i_angular)
!            write(use_unit,'(I6,5f13.8)') n_points, wave(1,n_points), wave(18,n_points), v_times_waves(35,n_points), v_times_waves(52,n_points)
!            write(use_unit,'(4f13.8)') partition_tab_2atoms(:,n_points) 
!            write(use_unit,'(A,4f13.8)') "i_r",i_r(:,i_angular)
!            write(use_unit,*)"wave",wave(1,n_points),wave(1,n_points)
!            write(use_unit,*)"vwave",v_times_waves(32,n_points),v_times_waves(34,n_points)
!            write(use_unit,*)"partition",partition(n_points)
!            write(use_unit,*) "radialwave",v_times_radialwaves_spl(1,:,134)
!            write(use_unit,*) "ylm",ylm_tab(:,:,i_angular)
!           endif

!     end if (partition_tab.gt.0)
             end if

!     end angular integration loop
          enddo
!     add non-relativistic contributions to the Hamiltonian matrix elements

          call evaluate_coulomb_matr_shell_v0 &
               ( n_points,i_task,i_atom, &
                 partition_tab_2atoms(1,1), &
                 n_compute, i_basis(1), &
                 n_compute_current, i_basis_current(1), &
                 n_loc_prodbas, i_prodbas(1), &
                 n_prodbas_current, i_prodbas_current(1), &
                 v_times_waves(1,1), &
                 wave(1,1), &
                 coulomb_matr(1,1))


!     end if (n_compute.gt.0) then
       end if
!       write(use_unit,'(2I6,3f13.8)') i_atom, i_radial, coulomb_matr(17,35),coulomb_matr(18,35),coulomb_matr(19,35)

!       write(110, '(100f16.8)') r_radial(i_radial, species(i_atom)),
!     +      ( x_matrix(i_basis_1, i_basis_1)/
!     +        w_radial(i_radial,species(i_atom)),
!     +        i_basis_1=1, n_basis)

!test
!            if (i_radial.eq.140) then
!            write(use_unit,*) "here 3"
!            end if
!test end
!    MPI task distribution
!       endif

!       end radial integration loop
      enddo

!     end integration loop over atoms
      enddo

!      call sync_matrix(coulomb_matr,n_basbas,n_basbas)

!   calculate O^-1 V O^-1
!$$ but, this is the old version, in the new "V" version, this is
!$$ not needed any more

!      do i_basis_1 =1, n_basbas
!        do i_basis_2 = i_basis_1 , n_basbas
!               coulomb_matr(i_basis_1,i_basis_2) =
!     +         0.5d0* ( coulomb_matr(i_basis_1,i_basis_2) +
!     +                  coulomb_matr(i_basis_2,i_basis_1) )

!               coulomb_matr(i_basis_2,i_basis_1) =
!     +           coulomb_matr(i_basis_1,i_basis_2)

!        enddo
!      enddo

!       open (100,file='coulomb.dat')
!       read(100,*) coulomb_matr(:,:)
!       close(100)
!       call coulomb_matr_transform (n_basbas, coulomb_matr,
!     +                             ovlp_prodbas)

!      do i_basis_1=1, n_basbas
!       write(use_unit,'(I6,f18.6)') i_basis_1,
!     +          coulomb_matr(i_basis_1,i_basis_1)
!      enddo

!      if(myid.eq.0) then
!      write(use_unit,'(2X,A,2X,A,7X,A)')
!     +          "i_basis_1", "i_basis_2", "Coulomb interaction matrix"

!       if(myid.eq.0) then
!       do i_basis_1 =1, 20, 1
!        do i_basis_2 =1, 20, 1
!          write(use_unit, '(4X,I5,5X,I5,9X,2f16.10)') &
!                     i_basis_1, i_basis_2, &
!                     coulomb_matr(i_basis_2,i_basis_1), &
!                     coulomb_matr(i_basis_1,i_basis_2)
!         enddo
!        enddo
!        do i_basis_1 = n_loc_prodbas/2+1, n_loc_prodbas, 1
!         do i_basis_2 = 1, n_basbas/2, 1
!        do i_basis_1 = 1, 20, 1
!         do i_basis_2 = 1, 20, 1
!          write(use_unit, '(4X,I4,5X,I4,2X,2I4,2X,3f16.10)') &
!                   i_basis_1-n_loc_prodbas/2, i_basis_2, &
!                   basbas_atom(i_basis_2),basbas_atom(i_basis_1), &
!!                   basbas_atom(i_basis_2+n_loc_prodbas/2),basbas_atom(i_basis_1-n_loc_prodbas/2), &
!                   coulomb_matr(i_basis_2,i_basis_1), &
!                   coulomb_matr(i_basis_1-n_loc_prodbas/2,i_basis_2+n_loc_prodbas/2), &
!                   -coulomb_matr(i_basis_2,i_basis_1)+ &
!                   coulomb_matr(i_basis_1-n_loc_prodbas/2,i_basis_2+n_loc_prodbas/2)
!         enddo
!        enddo
!       endif
!test


      if (allocated(ylm_tab)) then
        deallocate(ylm_tab)
      endif
      if (allocated(index_lm)) then
        deallocate(index_lm)
      endif
      if (allocated(v_times_radialwaves_spl)) then
        deallocate( v_times_radialwaves_spl )
      endif
      if (allocated( v_times_waves) ) then
        deallocate( v_times_waves )
      endif
      if (allocated( wave) ) then
        deallocate( wave )
      endif

      return
      end subroutine integrate_coulomb_matr_v0

!----------------------------------------------------------------------
!******
