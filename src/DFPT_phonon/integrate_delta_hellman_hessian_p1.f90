!------------------------------------------------------------------------------------------------------------

!****s* FHI-aims/integrate_delta_hellman_hessian_p1
!  NAME
!   integrate_delta_hellman_hessian_p1
!  SYNOPSIS

subroutine integrate_delta_hellman_hessian_p1 &
     ( hellman_feynman_term)


!  PURPOSE
!    The subroutine sums up the _delta_ hartree potential from the multipole components
!    which have been calculated before hand. The actual summation is done here.
!
!
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use free_atoms
  use spline
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use pbc_lists
  use load_balancing
  use hartree_potential_real_p0
  ! rho_multipole from hartree_potential_storage conflicts with the local variable here,
  ! so just only use get_rho_multipole_spl() from this module.
  use hartree_potential_storage, only : get_rho_multipole_supercell_spl , & 
                                 delta_v_hartree_part_at_zero_supercell, &
                                 delta_v_hartree_deriv_l0_at_zero_supercell



  implicit none

!  ARGUMENTS



  real*8, target, dimension(3,n_centers_in_sc_DFPT)        ::  hellman_feynman_term

!  INPUTS
!  OUTPUT
! o hellman_feynman_term
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

  integer index_lm(-l_pot_max:l_pot_max, 0:l_pot_max )

  real*8 coord_current(3)
  real*8 dist_tab_sq
  real*8 dist_tab_in
  real*8 dist_tab_out
  real*8 dir_tab(3)
  real*8 dir_tab_in(3)
  real*8 dir_tab_out(3)
  real*8 log_weight
  real*8 radial_weight
  real*8 trigonom_tab(4)
  real*8 i_r
  real*8 i_r_log
  real*8 ylm_tab((l_pot_max+1)**2)


  integer, parameter :: n_coeff_hartree = 2 ! Number of spline coeffs for current_delta_v_hart_part_spl

  integer l_h_dim

  
  real*8, dimension(:,:,:), allocatable :: current_delta_rho_multipole_spl
  real*8, dimension(:,:,:), allocatable :: current_delta_v_hart_part_spl
  real*8, dimension((l_pot_max+1)**2) :: delta_v_hartree_multipole_component
  real*8, dimension(:), allocatable :: delta_v_hartree_multipole_deriv
  real*8 v_hartree_gradient_temp(3), d_v_hartree_free_d_r
  
  !-------------------add for hellman_feynman gradient------------------
  real*8, dimension(:),allocatable :: adap_outer_radius_sq
  real*8, dimension(:,:), allocatable :: dylm_dtheta_tab
  real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab
  !-------------------end add for hellman_feynman gradient------------------







  integer :: current_spl_atom
  integer :: current_center
  integer :: atom_of_splines

  integer :: l_atom_max

  !  counters

  integer i_atom_2
  integer i_center, i_center_multipole
  integer i_batch
  integer i_l
  integer i_m
  integer i_coord
  integer i_index
  integer i_lm


  !  external functions
  real*8, external :: ddot
  character*100 :: info_str

  integer :: mpierr
  integer :: info



  ! Load balancing stuff

  integer n_my_batches_work  ! Number of batches actually used
  integer n_full_points_work ! Number of integration points actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

  ! Pointers to the actually used array


  integer n_bp

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all

  if(use_batch_permutation > 0) then
    call localorb_info("integrate_delta_hellman_hessian  load balancing.", use_unit,'(2X,A)', OL_norm )
  else
    call localorb_info("integrate_delta_hellman_hessian.", use_unit,'(2X,A)', OL_norm )
  endif


  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)

  n_bp = use_batch_permutation
!  if(use_batch_permutation > 0) then
!
!    n_my_batches_work = batch_perm(n_bp)%n_my_batches
!    n_full_points_work = batch_perm(n_bp)%n_full_points
!
!    batches_work => batch_perm(n_bp)%batches
!    partition_tab => batch_perm(n_bp)%partition_tab
!
!    allocate(rho(n_full_points_work))
!    call permute_point_array(n_bp,delta_rho_std,rho)
!!
!
!
!    allocate(potential(n_full_points_work))

!  else

    n_my_batches_work = n_my_batches
    n_full_points_work = n_full_points
    batches_work => batches

!  endif

  if(get_batch_weights) then
    allocate(batch_times(n_my_batches_work))
    batch_times(:) = 0
  endif

  !-----------------------------------------------------------------------------





  allocate(adap_outer_radius_sq(n_atoms),stat=info)
  call check_allocation(info, 'adap_outer_radius_sq          ')


  if (.not.allocated(current_delta_rho_multipole_spl)) then
     allocate(current_delta_rho_multipole_spl &
          ((l_pot_max+1)**2, n_max_spline, n_max_radial+2),stat=info)
     call check_allocation(info, 'current_delta_rho_multipole_spl     ')
  end if

  if (.not.allocated(current_delta_v_hart_part_spl)) then
     allocate(current_delta_v_hart_part_spl &
          ((l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid),stat=info)
     call check_allocation(info, 'current_delta_v_hart_part_spl ')
  end if

  if (.not.allocated(delta_v_hartree_multipole_deriv)) then
        allocate(delta_v_hartree_multipole_deriv((l_pot_max+1)**2),stat=info)
        call check_allocation(info, 'delta_v_hartree_multipole_deri')
  end if



  if (.not.allocated(dylm_dtheta_tab)) then
        allocate(dylm_dtheta_tab((l_pot_max+1)**2, n_centers_hartree_potential),stat=info)
        call check_allocation(info, 'dylm_dtheta_tab               ')

  end if

  if (.not.allocated(scaled_dylm_dphi_tab)) then
        allocate(scaled_dylm_dphi_tab((l_pot_max+1)**2,n_centers_hartree_potential),stat=info)
        call check_allocation(info, 'scaled_dylm_dphi_tab          ')

  end if




!--------------------------(1) near field : rho_multipole------------------

  call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()


    ! First loop over the grid: We run over the Hartree potential center by center, adding the Hartree
    ! contribution of that atom at each point of the integration grid

    do i_center = 1, n_centers_in_sc_DFPT, 1    ! ===> R1 in my equation

       i_atom_2 =  center_in_sc_DFPT_to_atom(i_center)
       hellman_feynman_term(1:3,i_center) = 0.0d0 
       !atom_of_splines = 0

    do i_center_multipole = 1, n_centers_in_sc_DFPT, 1  ! ===> R2 in my equation

       current_center   = center_in_sc_DFPT_to_center(i_center_multipole)
       current_spl_atom = center_to_atom(current_center)


            ! Tabulate distances and directions to all atoms -
            ! including relative positions on logarithmic and
            ! radial integration grids.

            call tab_two_center_coords_PBC & 
                 ( i_center_multipole, i_center,  & 
                 dist_tab_sq,  &
                 dir_tab ) ! = i_center-i_center_multipole = R1- R2 

            if (i_center_multipole.eq. i_center ) then
              dist_tab_sq = r_grid_min(species(current_spl_atom))**2 + 1e-15
            end if
          
              
            ! At each integration point, the Hartree potential components coming from
            ! different atoms are split into two groups:
            ! Any components coming from atoms close by are evaluated by explicit
            ! numerical splines; far away atoms are simply represented by
            ! an analytical long-distance multipole potential.
            !if (dist_tab_sq.lt.multipole_radius_sq(current_spl_atom)) then
               ! begin with everything related to n_atoms_in ...

               !if(i_center_multipole /= atom_of_splines) then
                call get_rho_multipole_supercell_spl(current_delta_rho_multipole_spl, i_center_multipole, & 
                                                     current_spl_atom)
                call integrate_delta_v_hartree(current_delta_rho_multipole_spl, current_delta_v_hart_part_spl, &
                                               n_coeff_hartree, current_spl_atom )
               !  atom_of_splines = i_center_multipole
               !endif


               call tab_single_atom_centered_coords_radial_log_p0 &
                    ( current_center, dist_tab_sq, dir_tab,  &
                    dist_tab_in, i_r, i_r_log, dir_tab_in )

               ! for all inner atoms, we need ylm functions and their gradients explicitly
               call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)

               call tab_single_gradient_ylm_p0 &
                       ( trigonom_tab, l_hartree, l_pot_max,  &
                       current_center,  &
                       ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab)

               ! Now loop over all inner atoms (those which are close enough
               ! to the current integration point so we need explicit numerical splines)
               !     partitioned
               !     hartree potentials need to be summed up
               !     according to Delley (eq. 12c)

               l_h_dim = (l_hartree(species_center(current_center))+1)**2

               ! obtain spline-interpolated values of the multipole components
               ! of the partitioned Hartree potential, splined on the logarithmic
               ! integration grid

               call spline_vector_v2 &
                    ( i_r_log, &
                    current_delta_v_hart_part_spl, &
                    (l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid, &
                    n_grid(species_center(current_center)), &
                    l_h_dim, &
                    delta_v_hartree_multipole_component)

!--------------------for delta_V(0)-------------------------
               if (i_center_multipole.eq. i_center ) then
                  delta_v_hartree_multipole_component(1) =    &  
                       delta_v_hartree_part_at_zero_supercell(current_spl_atom)

                  do i_lm = 2, l_h_dim

                     delta_v_hartree_multipole_component(i_lm) = 2* current_delta_v_hart_part_spl(i_lm,1,1) &
                          - current_delta_v_hart_part_spl(i_lm,1,2)

                  end do
               end if
!--------------------end for delta_V(0)-------------------------


                  ! call spline vector derivative
                  ! dot priduct.
                  ! abs(V_radial_deriv) * dir(:,i_center_L)

                  call tab_single_radial_weights_v2 &
                       ( current_spl_atom, dist_tab_in, i_r, &
                       log_weight, radial_weight )

                  ! splines are now derivatives df/di where i is the grid point index
                  ! must convert to df/dr = df/di * di/dr

                  call spline_deriv_vector_v2 &
                       ( i_r_log, current_delta_v_hart_part_spl, &
                       (l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid,  &
                       n_grid(species_center(current_center)), &
                       l_h_dim, &
                       delta_v_hartree_multipole_deriv)

                  delta_v_hartree_multipole_deriv(:) = &
                       delta_v_hartree_multipole_deriv(:) * log_weight


!--------------------for delta_V(0)-------------------------
                  if(i_center_multipole == i_center) then

                     ! l=1 components

      delta_v_hartree_multipole_deriv(2) = delta_v_hartree_deriv_l0_at_zero_supercell(1, i_center_multipole)
      delta_v_hartree_multipole_deriv(3) = delta_v_hartree_deriv_l0_at_zero_supercell(2, i_center_multipole)
      delta_v_hartree_multipole_deriv(4) = delta_v_hartree_deriv_l0_at_zero_supercell(3, i_center_multipole)

                  endif

!--------------------end for delta_V(0)-------------------------


                  !--------shanghui add for only delta_hellman_feynman_hessian---------------
                   v_hartree_gradient_temp = 0.0d0
                   d_v_hartree_free_d_r = 0.0d0
                  !--------shanghui end add only delta_hellman_feynman_hessian---------
                 

                  ! obtain gradients for present atom ...
                  call evaluate_v_hartree_gradient &
                       ( dist_tab_in, dir_tab_in,  &
                       trigonom_tab,  &
                       ylm_tab,  &
                       dylm_dtheta_tab,  &
                       scaled_dylm_dphi_tab,  &
                       delta_v_hartree_multipole_component,  &
                       delta_v_hartree_multipole_deriv, &
                       d_v_hartree_free_d_r, &
                       l_h_dim, &
                       v_hartree_gradient_temp)


              !-------shanghui add for DFPT_all_q:(1) near-field-----------------
               v_hartree_gradient_temp(1:3)=  &
               species_z(species(i_atom_2)) * v_hartree_gradient_temp(1:3)
                !write(use_unit,*) i_center_multipole,i_center,dir_tab 

                !dir_tab=i_center(R1) -i_center_multipole(R2) 

                !write(use_unit,*) 'print dV^R2(R1-R2)/dR1:'- v_hartree_gradient_temp(1:3)

                hellman_feynman_term(1:3, i_center)  = hellman_feynman_term(1:3, i_center) &
                - v_hartree_gradient_temp(1:3) 
              !-------shanghui end add for DFPT_all_q-----------------


            !end if ! multipole radius



    end do  ! end loop over i_center_multipole
    
    end do ! loop over i_center







  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_world,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for potential: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)


  if(get_batch_weights) then
    call set_batch_weights(n_bp,batch_times)
    deallocate(batch_times)
  endif





!  if(use_batch_permutation > 0) then
!
!    call permute_point_array_back(n_bp,1,potential,potential_std)
!    deallocate(potential)
!
!    deallocate(rho)
!
!  endif

  ! Write one last line to bound output
  write(info_str,*) ' '
  call localorb_info(info_str, use_unit,'(A)',OL_norm)


  if (allocated(delta_v_hartree_multipole_deriv)) then
     deallocate(delta_v_hartree_multipole_deriv)
  end if

  if (allocated(current_delta_rho_multipole_spl)) then
     deallocate(current_delta_rho_multipole_spl)
  end if
  if (allocated(current_delta_v_hart_part_spl)) then
     deallocate(current_delta_v_hart_part_spl)
  end if

  if (allocated(dylm_dtheta_tab)) then
     deallocate(dylm_dtheta_tab)
  end if
  if (allocated(scaled_dylm_dphi_tab)) then
     deallocate(scaled_dylm_dphi_tab)
  end if



end subroutine integrate_delta_hellman_hessian_p1
!******
!------------------------------------------------------------------------------------------------------------
