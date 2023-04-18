!****s* FHI-aims/integrate_ovlp3fn
!  NAME
!   integrate_ovlp3fn
!  SYNOPSIS

subroutine integrate_ovlp3fn(basis_l_max,ext_l_max, ovlp_3fn, ovlp_type)

  !  PURPOSE
  !
  !    This subroutine is intended to calculate the integration over three
  !    basis functions, revised from the subroutine
  !    "integrate_real_hamiltonian_matrix".
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
  use species_data, only: species_pseudoized
  implicit none

  !  ARGUMENTS
  integer, intent(IN) :: basis_l_max (n_species)
  integer, intent(IN) :: ext_l_max (n_species)
  real*8, intent(OUT) :: ovlp_3fn(n_basis_pairs, n_loc_prodbas)
  integer, intent(IN) :: ovlp_type

  !  INPUTS
  !  o ext_l_max -- integer array, the maximal angular momentum number of the auxiliary  
  !           basis functions for each species 
  !  o ovlp_type -- either OVLP_TYPE_OVERLAP, OVLP_TYPE_COULOMB, or
  !                 OVLP_TYPE_HSE
  !  OUTPUTS
  !  o ovlp_3fn -- real array, the three-center overlap integral (O integral) over 
  !           two NAO basis functions and one auxiliary basis function. This is the
  !           central quantity of the entire formalism of post-Hartree-Fock calculation.
  !           Later on, the is a transformation from ovlp_3fn to ovlp_3KS, the latter
  !           being the the integral over two single-particle orbital (KS or HF) and
  !           one auxiliary basis. 
  !
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE

  !  local variables

  integer :: l_ylm_max
  integer :: l_ylm_max_basbas
  integer, dimension(:,:), allocatable :: index_lm
  integer, dimension(:,:), allocatable :: index_lm_basbas
  ! we allocate different index_lm and ylm_tab for normal basis and
  ! auxiliary basis functions as they might differ
  real*8, dimension(:,:,:), allocatable :: ylm_tab
  real*8, dimension(:,:,:), allocatable :: ylm_tab_basbas
  real*8, dimension(:,:,:), allocatable :: v_times_radialwaves_spl
  real*8, dimension(:,:), allocatable :: v_times_waves

  real*8 coord_current(3, n_max_angular)
  real*8 dist_tab(n_atoms, n_max_angular)
  real*8 i_r(n_atoms, n_max_angular)
  real*8 dir_tab(3,n_atoms, n_max_angular)
  real*8 trigonom_tab(4,n_atoms, n_max_angular)

  real*8 wave(n_basis, n_max_angular)
  real*8 radial_wave(n_basis)

  integer :: n_compute
  integer :: n_compute_current
  integer :: n_prodbas_current
  integer :: n_prodbas_other
  integer :: i_basis(n_basis)
  integer :: i_basis_current(n_basis)
  integer :: i_prodbas(n_loc_prodbas)
  integer :: i_basbas(n_loc_prodbas)
  integer :: i_prodbas_current(n_loc_prodbas)
  integer :: i_prodbas_other(n_loc_prodbas)

  integer :: n_prod_compute

  integer :: n_compact(n_basis)
  integer :: compact_pairs(n_basis,n_basis) 
  integer :: n_compact_current(n_basis)
  integer :: compact_pairs_current(n_max_basis_atom, n_basis)

  integer :: basbas_l_max(n_species)
  integer :: normal_l_max(n_species)

  !     optimal accounting for matrix multiplications: only use points with nonzero components
  integer :: n_points
  integer :: n_zeros

  !     and condensed version of partition_tabs on angular grids
  real*8 :: partition_tab(n_max_angular)
  real*8 :: partition(n_max_angular)
  !      real*8 :: partition_tab_3atoms
  !     +          (n_atoms*(n_atoms+1)/2,n_max_angular)
  real*8 :: partition_tab_3atoms &
  (n_max_angular,n_atoms*(n_atoms+1)/2)

  !      real*8, dimension(:,:,:), allocatable :: gradient_basis_wave

  !      real*8, dimension(n_max_angular) :: zora_operator
  !      logical, dimension(n_max_angular) :: t_zora
  !      real*8, dimension(n_basis,3,n_max_angular) :: zora_vector

  real time_stamp
  real time_end

  real*8 ovlp_sq

  !     partition_type_temp makes sure that we use the right partition function for integrals
  integer :: partition_type_temp = 1


  !  counters

  integer i_basis_1
  integer i_basis_2
  integer i_basis_3
  integer i_basbas_1
  integer i_prodbas_1
  integer i_atom
  integer i_atom_1
  integer i_atom_2
  integer i_radial
  integer i_angular
  integer i_grid
  integer i_index, i_l, i_m
  integer i_coord

  integer i_species

  integer i_compute
  integer i_compute_1
  integer i_compute_2

  integer i_task
  integer i_basis_index
  logical :: do_coulomb
  character(*), parameter :: func = 'integrate_ovlp3fn'

  !  begin work

!  if(myid.eq.0)then
!     write(use_unit,*) 'ovlp_type', ovlp_type
!     write(use_unit,*) 'OVLP_TYPE_COULOMB',OVLP_TYPE_COULOMB 
!     write(use_unit,*) 'OVLP_TYPE_HSE',OVLP_TYPE_HSE 
!     write(use_unit,*) 'OVLP_TYPE_OVERLAP',OVLP_TYPE_OVERLAP 
!  endif

  select case (ovlp_type)
  case(OVLP_TYPE_COULOMB)
     do_coulomb = .true.
     if (use_hse.and.(hse_omega_hf.ne.0.0d0).and. .not. use_gw_and_hse &
         .and. .not. use_dftpt2_and_hse) then
        call aims_stop('Cannot calculate bare Coulomb matrix within HSE', func)
     end if
  case(OVLP_TYPE_HSE)
     do_coulomb = .true.
     if (.not. (use_hse.and.(hse_omega_hf.ne.0.0d0)).or.  use_gw_and_hse) then
        call aims_stop('Cannot calculate screened Coulomb matrix without HSE', func)
     end if
  case(OVLP_TYPE_LR)
     do_coulomb = .true.
     if ((.not. (use_hse.and.(hse_omega_hf.ne.0.0d0)) .or.  use_gw_and_hse) &
         .and. .not. lrc_pt2_started) then
        call aims_stop('Cannot calculate long-range Coulomb matrix without HSE or LRC-PT2', func)
     end if
  case(OVLP_TYPE_OVERLAP)
     do_coulomb = .false.
  case default
     write(use_unit,*) "ovlp_type: ", ovlp_type
     call aims_stop('Invalid ovlp_type')
     stop   ! Satisfy compiler
  end select


  if(myid.eq.0) then
     write(use_unit,*)
     write(use_unit,*)"----------------------------------------------------"
     write (use_unit,'(2X,A,A,A)') &
     "Integrating the 3-basis-function ", &
     & trim(OVLP_TYPE_NAMES(ovlp_type)), " matrix ..."
  endif

  !     partition_type_temp makes sure that we use the right partition function for integrals
  if (partition_type.eq.6) then
     partition_type_temp = 6
  else
     ! just use type 1 anyway
     partition_type_temp = 1
  end if

  !     begin with general allocations
  l_ylm_max = 2*l_wave_max
  allocate( ylm_tab( (l_ylm_max+1)**2, n_atoms, &
  n_max_angular) )
  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max) )
  
  if (flag_auxil_basis == PRODBAS_OPT) then
    l_ylm_max_basbas = 2*maxval(basbasfn_l)
  else
    l_ylm_max_basbas = 2*l_ext_max  
  endif

  allocate( ylm_tab_basbas( (l_ylm_max_basbas+1)**2, n_atoms, &
  n_max_angular) )
  allocate( index_lm_basbas( -l_ylm_max_basbas:l_ylm_max_basbas, 0:l_ylm_max_basbas) )
  allocate( v_times_waves (n_loc_prodbas, n_max_angular) )

  !     first integrate over v(r,r')*P_u(r') on the logrithmic radial gird

  if (do_coulomb) then
     allocate( v_times_radialwaves_spl ( &
     n_max_spline, n_hartree_grid, n_loc_prodbas) )
     call integrate_v_times_radialwaves &
     (v_times_radialwaves_spl)
  end if

  !     initialize

  ovlp_3fn(:,:) = 0.d0

  !     initialize index_lm

  i_index = 0
  do i_l = 0, l_ylm_max, 1
     do i_m = -i_l, i_l
        i_index = i_index+1
        index_lm(i_m,i_l) = i_index
     enddo
  enddo

  i_index = 0
  do i_l = 0, l_ylm_max_basbas, 1
     do i_m = -i_l, i_l
        i_index = i_index+1
        index_lm_basbas(i_m,i_l) = i_index
     enddo
  enddo

  normal_l_max = 2*basis_l_max

  if (flag_auxil_basis == PRODBAS_OPT) then
    basbas_l_max = 2*maxval(basbasfn_l)
  else
    basbas_l_max = 2*ext_l_max  
  endif

  i_task = myid + 1
  ! perform partitioned integration, atom by atom, and point by point
  ! This will be the outermost loop, to save evaluations of the potential.
  ! and the Y_lm functions
  do i_atom = 1, n_atoms, 1
     if(species_pseudoized(species(i_atom))) cycle
     !test
     if(myid.eq.0) then
        write(use_unit,*) " | i_atom: ", i_atom
     endif
     !test end

     do i_radial = 1, n_radial(species(i_atom)), 1

        ! if (myid.eq.radial_task_list(i_radial,i_atom)) then

        n_compute = 0
        i_basis = 0

        n_compute_current = 0
        i_basis_current = 0

        i_prodbas = 0
        i_basbas = 0

        n_prodbas_current = 0
        i_prodbas_current = 0

        n_prodbas_other = 0
        i_prodbas_other = 0

        n_prod_compute = 0
        do i_angular = 1, n_angular(i_radial, species(i_atom)), 1

           ! get current integration point coordinate
           do i_coord = 1, 3, 1
              coord_current(i_coord,i_angular) = &
              coords(i_coord,i_atom )  + &
              r_angular(i_coord, i_angular, i_radial, &
              species(i_atom)) * &
              r_radial(i_radial, species(i_atom))
           enddo

           ! compute atom-centered coordinates of current integration point,
           ! as viewed from all atoms
           call tab_atom_centered_coords &
           ( coord_current(1,i_angular), &
           dist_tab(1,i_angular), i_r(1,i_angular), &
           dir_tab(1,1,i_angular) &
           )

           ! evaluate the partition function for three atoms
           call evaluate_partition_tab_3atoms &
           (i_atom, dist_tab(1:n_atoms,i_angular), &
           i_r(1:n_atoms,i_angular), &
           w_radial( i_radial, species (i_atom)), &
           w_angular( i_angular, i_radial, species (i_atom)), &
           partition_tab_3atoms(i_angular, &
           1:n_atoms*(n_atoms+1)/2), &
           partition_type_temp )

           partition_tab(i_angular) = 0.d0
           do i_atom_1 =1, n_atoms*(n_atoms+1)/2, 1
              partition_tab(i_angular) = &
              max( partition_tab(i_angular), &
              partition_tab_3atoms(i_angular,i_atom_1) )
           enddo
           do i_atom_1 = 1, n_atoms
              if(dist_tab(i_atom_1,i_angular) .lt. 1.e-15) then
                 partition_tab(i_angular) = 0.d0
                 exit
              endif
           enddo

           !     determine which basis functions are relevant at current integration point,
           !     and tabulate their indices
           if (partition_tab(i_angular).gt.0.d0) then

              call prune_basis_v1(dist_tab(1,i_angular), n_compute, i_basis)
              call prune_prodbas_v1 &
              ( i_task, dist_tab(1,i_angular), n_prod_compute, i_prodbas )

              do i_compute = 1, n_prod_compute
                 i_prodbas_1 = i_prodbas(i_compute)
                 i_basbas(i_compute) = map_prodbas(i_prodbas_1,myid+1)
              enddo

           end if
        enddo   ! i_angular

        ! from the relevant basis function set i_basis find those which belong to
        ! the atom i.
        call condense_compute_pairs &
        ( i_atom, n_compute, i_basis, &
        n_compute_current, i_basis_current, &
        n_compact, compact_pairs, &
        n_compact_current, compact_pairs_current &
        )

        ! Fixme: in future this has to changed to exploit some kind of locality
        !        for auxiliary basis.
        do i_compute_1 = 1, n_prod_compute
           i_basis_1 = i_prodbas(i_compute_1)
           i_basbas_1 = map_prodbas(i_basis_1, i_task)
           if(i_basbas_1.gt.0) then
              if(basbas_atom(i_basbas_1).eq.i_atom) then
                 n_prodbas_current = n_prodbas_current + 1
                 i_prodbas_current(n_prodbas_current) = i_compute_1
              else
                 n_prodbas_other = n_prodbas_other + 1
                 i_prodbas_other(n_prodbas_other) = i_compute_1
              endif
           endif
        enddo

        if (n_compute.gt.0) then

           n_points = 0
           partition(:) = 0
           do i_angular = 1, n_angular(i_radial, species(i_atom)), 1

              if (partition_tab(i_angular).gt.0.d0) &
              then
                 ! execute only if partition_tab.gt.0 here, i.e. if the integration point
                 ! makes sense
                 n_points = n_points + 1
                 partition_tab_3atoms(n_points,1:n_atoms*(n_atoms+1)/2) = &
                 partition_tab_3atoms(i_angular,1:n_atoms*(n_atoms+1)/2)
                 ! compute trigonometric functions of spherical coordinate angles
                 ! of current integration point, viewed from all atoms
                 call tab_trigonom &
                 ( dir_tab(1,1,i_angular), &
                 trigonom_tab(1,1,i_angular) &
                 )

                 ! tabulate distance and Ylm's w.r.t. other atoms
                 call tab_wave_ylm &
                 ( trigonom_tab(1,1,i_angular), normal_l_max, &
                 l_ylm_max, &
                 ylm_tab(1,1,i_angular) )
                 
                 ! tabulate distance and Ylm's w.r.t. other atoms
                 call tab_wave_ylm &
                 ( trigonom_tab(1,1,i_angular), basbas_l_max, &
                 l_ylm_max_basbas, &
                 ylm_tab_basbas(1,1,i_angular) )                 


                 ! tabulate total wave function value for each basis function
                 call evaluate_waves_v0 &
                 (i_r(1,i_angular), l_ylm_max, &
                 ylm_tab(1,1,i_angular), &
                 dist_tab(1,i_angular), index_lm, n_compute, &
                 i_basis, radial_wave, &
                 wave(1,n_points))

                 if (do_coulomb) then
                    call evaluate_v_times_waves &
                    (l_ylm_max_basbas, ylm_tab_basbas(1,1,i_angular), &
                    dist_tab(1,i_angular), index_lm_basbas, &
                    n_prod_compute, i_prodbas, &
                    v_times_radialwaves_spl, &
                    v_times_waves(1:n_loc_prodbas,n_points))
                 else
                    call evaluate_prod_waves &
                    ( i_r(1,i_angular), l_ylm_max_basbas, ylm_tab_basbas(1,1,i_angular), &
                    dist_tab(1,i_angular), index_lm_basbas, &
                    n_prod_compute, i_basbas, &
                    v_times_waves(1:n_loc_prodbas,n_points) &
                    )
                 end if


              end if   ! if (partition_tab.gt.0)
           enddo   ! angular integration loop

           !     add the contribution from the current shell
           call evaluate_ovlp3fn_shell_v0 &
           ( n_points,i_task,i_atom, &
           partition_tab_3atoms(1:n_points, &
           1:n_atoms*(n_atoms+1)/2), &
           n_prod_compute, i_prodbas, &
           n_prodbas_current, i_prodbas_current, &
           n_prodbas_other, i_prodbas_other, &
           n_compute, i_basis(1:n_compute), &
           n_compute_current, i_basis_current, &
           n_compact, n_compact_current, &
           compact_pairs, compact_pairs_current, &
           wave(1:n_compute,1:n_points), &
           v_times_waves, ovlp_3fn )

        end if  ! if (n_compute.gt.0)
     enddo   ! radial integration
  enddo   ! atom integration

  ! call symmetrize_ovlp3basis_matrix(ovlp_3fn)

  if (allocated(ylm_tab)) deallocate(ylm_tab)
  if (allocated(index_lm)) deallocate(index_lm)
  if (allocated(index_lm_basbas)) deallocate(index_lm_basbas)
  if (allocated(v_times_radialwaves_spl)) then
     deallocate(v_times_radialwaves_spl)
  endif
  if (allocated(v_times_waves)) deallocate(v_times_waves)

!  do i_basis_1 = 1, n_basis, 1
!    do i_basis_2 = 1, i_basis_1, 1
!      i_index = basis_nghr(i_basis_2,i_basis_1)
!      do i_prodbas_1 = 1, n_loc_prodbas, 1
!        if(basis_atom(i_basis_1).eq.2 .and. basis_atom(i_basis_2).eq.1 .and. basbas_atom(i_prodbas_1) .eq. 3) then
!         write(use_unit,'(9I4,f18.10)') basis_atom(i_basis_1), basis_atom(i_basis_2), basbas_atom(i_prodbas_1), &
!                basis_l(i_basis_1), basis_l(i_basis_2), basbas_l(i_prodbas_1), &
!                basis_m(i_basis_1), basis_m(i_basis_2), basbas_m(i_prodbas_1), ovlp_3fn(i_index, i_prodbas_1)
!        endif
!      enddo
!    enddo
!  enddo
end subroutine integrate_ovlp3fn
!******
