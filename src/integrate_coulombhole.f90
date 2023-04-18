!****s* FHI-aims/integrate_coulombhole
!  NAME
!   integrate_coulombhole
!  SYNOPSIS

      subroutine integrate_coulombhole &
      ( n_KS_states, &
        partition_tab, basis_l_max, &
        screened_coulomb_matr, KS_eigenvector, &
        coulomb_hole)

!  PURPOSE
!  Subroutine evaluate_coulombhole_shell evaluates the static coulomb
!  hole part of the self energy:
!     sig_coh = < i | W(r,r) |i >, |i> being the KS state,
!
!  USES

      use dimensions
      use runtime_choices
      use grids
      use geometry
      use basis
      use prodbas
      use constants
      use mpi_tasks
      use localorb_io, only : use_unit

      implicit none

!  ARGUMENTS

      integer n_KS_states
      real*8, dimension(n_full_points) :: partition_tab
      integer basis_l_max (n_species)
      real*8, dimension(n_basbas,n_basbas) :: screened_coulomb_matr
      real*8, dimension(n_basis,n_states, n_spin) :: KS_eigenvector

      real*8, dimension(n_KS_states, n_spin) :: coulomb_hole

!  INPUTS
!  o n_KS_states -- number of single particle states for which the static Coulomb hole
!        contribution to the self-energy nedds to be calculated
!  o partition_tab -- the partition function for evaluating the real space integrals
!  o basis_l_max -- the maximal angular momentum components for the basis functions
!          of each species
!  o screened_coulomb_matr -- the screend Coulomb interacition matrix (with auxiliary
!          basis) at the frequency zero
!  o KS_eigenvector -- the single particle (KS or HF) eigenvector
!
!  OUTPUTS
!  o coulomb_hole -- the static Coulomb hole contribution of the self-energy
!          for each KS/HF state. 
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


      real*8 coord_current(3, n_max_angular)
      real*8 dist_tab(n_atoms, n_max_angular)
      real*8 i_r(n_atoms, n_max_angular)
      real*8 dir_tab(3,n_atoms, n_max_angular)
      real*8 trigonom_tab(4,n_atoms, n_max_angular)

      real*8 phi_times_psi &
                (n_basis*n_KS_states, n_max_angular)
      real*8 radial_wave(n_basis)
      real*8 wave(n_basis, n_max_angular)
      real*8 KS_wave(n_KS_states,n_max_angular,n_spin)

      real*8 wave_prod(n_basbas, n_max_angular)

      integer :: n_compute
      integer :: n_prod_compute
      integer :: i_basis(n_basis)
      integer :: i_basbas(n_basbas)
      integer :: basbas_l_max(n_species)

!     optimal accounting for matrix multiplications: only use points with nonzero components
      integer :: n_points

!     and condensed version of partition_tabs on angular grids
       real*8 :: partition(n_max_angular)
!      real*8 :: partition_tab_3atoms
!     +          (n_atoms*(n_atoms+1)/2,n_max_angular)

!      real*8, dimension(:,:,:), allocatable :: gradient_basis_wave

!      real*8, dimension(n_max_angular) :: zora_operator
!      logical, dimension(n_max_angular) :: t_zora
!      real*8, dimension(n_basis,3,n_max_angular) :: zora_vector

      real time_stamp
      real time_end

!  counters

      integer i_basis_1
      integer i_basis_2
      integer i_basis_3
      integer i_basbas_1
      integer i_state
      integer i_atom
      integer i_atom_1
      integer i_atom_2
      integer i_radial
      integer i_angular
      integer i_grid
      integer i_index, i_l, i_m
      integer i_coord
      integer i_spin

      integer i_species

      integer i_compute
      integer i_compute_2

      integer i_full_points
      integer i_full_points_2

!  begin work

      write(use_unit,*)
      write(use_unit,'(2X,A,A)') &
        "Integrating the static Coulomb hole part of the self energy ", &
          "..."

!     begin with general allocations
        l_ylm_max = 2*l_wave_max

      allocate( ylm_tab( (l_ylm_max+1)**2, n_atoms, &
           n_max_angular) )
      allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max) )

!     initialize

      coulomb_hole(:,:)=0

!     initialize index_lm

      i_index = 0
      do i_l = 0, 2*l_wave_max, 1
        do i_m = -i_l, i_l
          i_index = i_index+1
          index_lm(i_m,i_l) = i_index
        enddo
      enddo

      basbas_l_max = 2*basis_l_max

      i_full_points = 0
      i_full_points_2 = 0
!     perform partitioned integration, atom by atom, and point by point
!     This will be the outermost loop, to save evaluations of the potential.
!     and the Y_lm functions
      do i_atom = 1, n_atoms, 1

!test
!        write(use_unit,*) "i_atom: ", i_atom
!test end

         do i_radial = 1, n_radial(species(i_atom)), 1


!test
!          write(use_unit,*) "  i_radial: ", i_radial
!test end

           n_compute = 0
           i_basis = 0
           n_prod_compute = 0
           i_basbas = 0


           do i_angular = 1, n_angular(i_radial, species(i_atom)), 1

            i_full_points_2 = i_full_points_2 + 1
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

!     determine which basis functions are relevant at current integration point,
!     and tabulate their indices
              if (partition_tab(i_full_points_2).gt.0.d0) &
                   then

                 call prune_basis_v1(dist_tab(1,i_angular), n_compute, &
                      i_basis)
                 call prune_prodbas_v1(myid+1, dist_tab(1,i_angular), &
                             n_prod_compute, i_basbas)

              end if

!     evaluate the partition function for three atoms
!     now used for this moment, it is not obvious how to exploit this property
!     for product basis


           enddo

!     from the relevant basis function set i_basis find those which belong to
!     the atom i.

!           i_basbas = 0
!           n_prod_compute = 0

!           do i_basbas_1 = 1, n_basbas, 1

!             if (basbas_atom(i_basbas_1).eq.i_atom) then

!               n_prod_compute = n_prod_compute +1
!               i_basbas (n_prod_compute)= i_basbas_1
!             endif

!           enddo
!              write(use_unit,*)"2",n_prod_compute
!              write(use_unit,*)i_basbas(:)
!test
!     if (i_radial.eq.140) then
!     write(use_unit,*) "here 2"
!     end if
!test end

           if (n_prod_compute.gt.0.and.n_compute.gt.0) then

              n_points = 0
              KS_wave(:,:,:) = 0.d0
              do i_angular = 1, n_angular(i_radial, species(i_atom)), 1

                 i_full_points = i_full_points + 1
                 if (partition_tab(i_full_points).gt.0.d0) &
                      then
!     execute only if partition_tab.gt.0 here, i.e. if the integration point
!     makes sense
                    n_points = n_points + 1
                    partition(n_points) = &
                         partition_tab(i_full_points)
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

!           tabulate total wave function value for each basis function
                    call evaluate_waves_v0 &
                         (i_r(1,i_angular), l_ylm_max, &
                         ylm_tab(1,1,i_angular), &
                         dist_tab(1,i_angular), index_lm, n_compute, &
                         i_basis, radial_wave(1), &
                         wave(1,n_points))

                    call evaluate_prod_waves &
                         (i_r(1,i_angular), l_ylm_max, &
                         ylm_tab(1,1,i_angular), &
                         dist_tab(1,i_angular), index_lm, &
                         n_prod_compute, &
                         i_basbas, &
                         wave_prod(1,n_points))


!   evaluate KS wave function at this point
                    do i_spin = 1, n_spin
                     do i_state = 1, n_KS_states, 1
                       do  i_compute = 1, n_compute, 1
                         KS_wave (i_state,n_points,i_spin) = &
                          KS_wave(i_state,n_points,i_spin) + &
                          wave(i_compute,n_points) * &
                          KS_eigenvector( i_basis(i_compute), &
                                          i_state,i_spin )
                       enddo
                     enddo
                    enddo

!     end if (partition_tab.gt.0)
             end if

!     end angular integration loop
          enddo
!     integrate for this shell
          do i_spin = 1, n_spin
            call evaluate_coulombhole_shell &
               ( n_points, n_KS_states, &
                 partition(1:n_points), &
                 n_prod_compute, &
                 i_basbas(1:n_prod_compute), &
                 KS_wave(1:n_KS_states,1:n_points,i_spin), &
                 wave_prod(1:n_prod_compute,1:n_points), &
                 screened_coulomb_matr, coulomb_hole(1,i_spin))
         enddo

      else

        i_full_points = i_full_points + &
                     n_angular(i_radial, species(i_atom))
!     end if (n_compute.gt.0) then
       end if

!test
!            if (i_radial.eq.140) then
!            write(use_unit,*) "here 3"
!            end if
!test end

!test
!            if ( (i_atom.eq.1) .and.
!     +           ((i_radial.eq.47).or.(i_radial.eq.48)) .and.
!     +           (i_angular.eq.1) ) then
!              write(use_unit,*) "      n_compute: ", n_compute
!              write(use_unit,*) "      i_basis  : ",
!     +          (i_basis(i_compute), i_compute = 1,n_compute,1)
!            end if
!test end
!       end radial integration loop
      enddo

!     end integration loop over atoms
      enddo

!test


      write(use_unit,'(5X,A,7X,A)') &
               "i_state", "Sig_COH"
      do i_spin = 1, n_spin
        do i_state = 1, n_KS_states

           write(use_unit, '(7X,I5,5X,F16.10)') &
           i_state,  coulomb_hole(i_state,i_spin)

        enddo
      enddo
!ctest

      if (allocated(ylm_tab)) then
        deallocate(ylm_tab)
      end if
      if (allocated(index_lm)) then
        deallocate(index_lm)
      end if
      return
      end subroutine integrate_coulombhole

!----------------------------------------------------------------------
!***********
