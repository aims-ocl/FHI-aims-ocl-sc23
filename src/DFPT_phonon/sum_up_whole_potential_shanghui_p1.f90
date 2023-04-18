!------------------------------------------------------------------------------------------------------------

!****s* FHI-aims/sum_up_whole_potential_shanghui_p1
!  NAME
!   sum_up_whole_potential_shanghui_p1
!  SYNOPSIS

subroutine sum_up_whole_potential_shanghui_p1 &
     (partition_tab_std, delta_rho_std, delta_potential_std)


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
  use hartree_potential_storage, only : get_rho_multipole_supercell_spl !,get_rho_multipole_spl 

  implicit none

!  ARGUMENTS


  integer  i_atom
  real*8, target, dimension(n_full_points)        :: partition_tab_std
  real*8, target, dimension(n_full_points)        :: delta_rho_std
  real*8, target, dimension(n_full_points)        :: delta_potential_std

!  INPUTS
! o partition_tab_std -- values of partition function
! o delta_rho_std -- delta_electron density
!
!  OUTPUT
! o delta_potential_std -- delta_ Hartree potential
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

  logical forces_on 

  real*8 :: hartree_multipole_error

  ! work arrays
  real*8, allocatable :: delta_rho_multipole(:)

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


  real*8, dimension(:,:,:), allocatable :: current_delta_rho_multipole_spl
  real*8, dimension(:,:,:), allocatable :: current_delta_v_hart_part_spl

  !     for spline_vector
  real*8, dimension((l_pot_max+1)**2) :: delta_rho_multipole_component
  real*8, dimension((l_pot_max+1)**2) :: delta_v_hartree_multipole_component
  integer l_h_dim

  real*8 delta_rho_multipole_aux
  real*8 delta_v_hartree_aux

  real*8, dimension(:),allocatable :: adap_outer_radius_sq

  integer :: current_spl_atom
  integer :: current_center
  integer :: atom_of_splines

  integer :: l_atom_max

  !  counters

  integer i_center_multipole, i_center_multipole_in_hamiltonian
  integer i_batch
  integer i_l
  integer i_m
  integer i_coord
  integer i_index,i_index_inside
  integer i_lm
  integer :: i_spin

  integer :: i_lat
  integer :: i_full_points

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

  real*8, pointer :: partition_tab(:)
  real*8, pointer :: delta_rho(:)
  real*8, pointer :: delta_potential(:)

  integer n_bp

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all


  if(use_batch_permutation > 0) then
    call localorb_info("Summing up the Hartree potential with load balancing.", use_unit,'(2X,A)', OL_norm )
  else
    call localorb_info("Summing up the Hartree potential.", use_unit,'(2X,A)', OL_norm )
  endif


    forces_on=.false.
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
    partition_tab => partition_tab_std

    delta_rho => delta_rho_std
    delta_potential => delta_potential_std

!  endif

  if(get_batch_weights) then
    allocate(batch_times(n_my_batches_work))
    batch_times(:) = 0
  endif

  !-----------------------------------------------------------------------------




  allocate(delta_rho_multipole(n_full_points_work),stat=info)
  call check_allocation(info, 'delta_rho_multipole')
  

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




!----------------------------(2) far field multipole moment------------
  !  initialize index_lm
  i_index = 0
  do i_l = 0, l_pot_max, 1
     do i_m = -i_l, i_l
        i_index = i_index + 1
        index_lm(i_m, i_l) = i_index
     enddo
  enddo






!--------------------------(1) near field : rho_multipole------------------

  ! Initialize potential (which is an output variable and runs over all grid points
  !delta_potential   = 0.d0

  delta_rho_multipole = 0.d0

  call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()


    ! First loop over the grid: We run over the Hartree potential center by center, adding the Hartree
    ! contribution of that atom at each point of the integration grid

!    do i_center = 1, n_centers_hartree_potential, 1
!      current_center   = centers_hartree_potential(i_center)
!      current_spl_atom = center_to_atom(current_center)



     do i_center_multipole = 1,n_centers_in_sc_DFPT

        current_center = center_in_sc_DFPT_to_center(i_center_multipole)
        current_spl_atom = center_to_atom(current_center)
        atom_of_splines = 0
        !if(myid.eq.0) write(use_unit,*) 'i_center_multipole:',i_center_multipole
      !
      ! Now follows the real work: Summing the multipole potential, density, and their
      ! derivatives on the integration grid
      !

      ! Reset grid counter for current Hartree potential center
      i_full_points = 0

      do i_batch = 1, n_my_batches_work

        if(get_batch_weights) time_start = mpi_wtime()

        ! loop over one batch
        do i_index = 1, batches_work(i_batch)%size, 1

          ! i_full_points is the index that indicates where we are in the entire grid (for external quanities like rho, potential, ...)
          i_full_points = i_full_points + 1

          ! Only execute if partition_tab is .gt. zero, else
          ! we can run into 1/0 !!!
          if (partition_tab(i_full_points).gt.0.d0) then

            ! get current integration point coordinate
            coord_current(:) = batches_work(i_batch) % points(i_index) % coords(:)

            ! Tabulate distances and directions to a single atom -
            ! including relative positions on logarithmic and
            ! radial integration grids.

            call tab_single_atom_centered_coords_p0 &
                 ( current_center, &
                 coord_current,  &
                 dist_tab_sq,  &
                 dir_tab )

            ! VB: For uniformity, determine the maximum angular momentum required for the present atom 
            ! right here

            l_atom_max = l_hartree(species(current_spl_atom))
!            do while ( (outer_potential_radius(l_atom_max, current_spl_atom) .lt. dist_tab_sq ) & 
!                 .and. (l_atom_max.gt.0) ) 
!              l_atom_max = l_atom_max - 1
!            enddo


            ! At each integration point, the Hartree potential components coming from
            ! different atoms are split into two groups:
            ! Any components coming from atoms close by are evaluated by explicit
            ! numerical splines; far away atoms are simply represented by
            ! an analytical long-distance multipole potential.

            ! We treat first the atoms close by
 
!            if (dist_tab_sq.lt.multipole_radius_sq(current_spl_atom) ) then
          !  if (dist_tab_sq.lt.131.648909d0 ) then

            !----shanghui treat everything using MP--------------
            !if (dist_tab_sq.lt.max(adap_outer_radius_sq(current_spl_atom),multipole_radius_sq_c(current_spl_atom)) ) then 

              if(i_center_multipole /= atom_of_splines) then
                !call get_rho_multipole_spl(current_delta_rho_multipole_spl, current_spl_atom)
                call get_rho_multipole_supercell_spl(current_delta_rho_multipole_spl, i_center_multipole, & 
                                                     current_spl_atom)
                call integrate_delta_v_hartree(current_delta_rho_multipole_spl, current_delta_v_hart_part_spl, &
                                               n_coeff_hartree, current_spl_atom )
                atom_of_splines = i_center_multipole 
               ! to avoid spl at every i_full_points, which could be costy
              endif
              

              call tab_single_atom_centered_coords_radial_log_p0 &
                   ( current_center, dist_tab_sq, dir_tab,  &
                   dist_tab_in, i_r, i_r_log, dir_tab_in )

              ! for an inner atom we need ylm functions and their gradients explicitly

              call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)

              call tab_single_wave_ylm_p2 &
                     ( trigonom_tab, l_atom_max,  &
                     l_pot_max, ylm_tab)


              ! For an inner atoms (those which are close enough
              ! to the current integration point so we need explicit numerical splines)
              !     partitioned
              !     hartree potentials need to be summed up
              !     according to Delley (eq. 12c)
              l_h_dim = (l_atom_max + 1)**2

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

                ! sum up the Hartree potential contribution from the present inner atom
!                delta_v_hartree_aux = &
!                       ddot ( l_h_dim, delta_v_hartree_multipole_component, 1, ylm_tab, 1 )
                delta_v_hartree_aux = 0.0d0
 
                do i_index_inside =1, l_h_dim
                delta_v_hartree_aux = delta_v_hartree_aux + & 
                delta_v_hartree_multipole_component(i_index_inside)* & 
                ylm_tab(i_index_inside)
                enddo


              !-------shanghui add for DFPT_all_q:(1) near-field-----------------
                delta_potential( i_full_points) = &
                delta_potential( i_full_points) + &
                delta_v_hartree_aux
              !-------shanghui end add for DFPT_all_q-----------------


              ! Obtain spline-interpolated values of the multipole density,
              ! this time splined on the radial integration grid
              call spline_vector_v2 &
                   ( i_r+1, current_delta_rho_multipole_spl, &
                   (l_pot_max+1)**2, n_max_spline, n_max_radial+2,  &
                   n_radial(species_center(current_center))+2,  &
                   l_h_dim, &
                   delta_rho_multipole_component)

                ! sum up the multipole density contribution from the present inner atom
                !delta_rho_multipole_aux = &
                !       ddot ( l_h_dim, delta_rho_multipole_component, 1, ylm_tab, 1)
                delta_rho_multipole_aux = 0.0d0
 
                do i_index_inside =1, l_h_dim
                delta_rho_multipole_aux = delta_rho_multipole_aux + & 
                delta_rho_multipole_component(i_index_inside)* & 
                ylm_tab(i_index_inside)
                enddo

                delta_rho_multipole(i_full_points) = & 
                delta_rho_multipole(i_full_points) + & 
                delta_rho_multipole_aux


!            else if (dist_tab_sq .lt. adap_outer_radius_sq(current_spl_atom)) then
!              ! the current center is in the far-distance part-----------------
!
!



           !end if  ! end if for separating current_center either far-distance or near part




          end if ! if (partition_tab(i_full_points).gt.0)

        end do ! end loop over a batch
        if(get_batch_weights) batch_times(i_batch) = batch_times(i_batch) + mpi_wtime() - time_start
      end do ! end loop over batches
      
 
     

    enddo ! i_center_multipole 


    i_full_points = 0 
    hartree_multipole_error = 0.d0
    do i_batch = 1, n_my_batches_work
    do i_index = 1, batches_work(i_batch)%size, 1

       i_full_points = i_full_points + 1

       if (partition_tab(i_full_points).gt.0.d0) then

            hartree_multipole_error = &
            hartree_multipole_error  + &
          (abs( delta_rho(i_full_points) - delta_rho_multipole(i_full_points) ))**2 &
                * partition_tab(i_full_points)
        endif 
     enddo 
     enddo 
 


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



 !-------shanghui begin parallel------  
  call sync_sum_up_whole_potential_shanghui(  &
       hartree_multipole_error)
 !-------shanghui end parallel------  



  ! write root-mean square error of Hartree potential
  hartree_multipole_error =   &
       sqrt(hartree_multipole_error)

     write(info_str,'(2X,A,1X,E14.6)')  &
          "| RMS delta_rho error from multipole expansion :",   &
          hartree_multipole_error
     call localorb_info(info_str, use_unit,'(A)', OL_norm )



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

  if(allocated(delta_rho_multipole)) deallocate(delta_rho_multipole)

  if (allocated(current_delta_rho_multipole_spl)) then
     deallocate(current_delta_rho_multipole_spl)
  end if
  if (allocated(current_delta_v_hart_part_spl)) then
     deallocate(current_delta_v_hart_part_spl)
  end if


end subroutine sum_up_whole_potential_shanghui_p1
!******
!------------------------------------------------------------------------------------------------------------
