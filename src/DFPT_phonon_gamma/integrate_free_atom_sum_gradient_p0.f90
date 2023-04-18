!------------------------------------------------------------------------------------------------------------

!****s* FHI-aims/integrate_free_atom_sum_gradient_p0
!  NAME
!   integrate_free_atom_sum_gradient_p0
!  SYNOPSIS

subroutine integrate_free_atom_sum_gradient_p0 &
     (partition_tab_std, rho_free_gradient,v_free_gradient)

!  PURPOSE
!    The subroutine get the gradient of total free atom density ( rho_free_gradient )
!    and gradient of total free atom potential ( v_free_gradient ) for phonon_gamma
!
!    shanghui 2013.12.30
!
!  USES

  use constants, only: pi4_inv
  use dimensions
  use runtime_choices
  use grids
  use geometry   !species(n_atom)
  use species_data
  use free_atoms
  use spline
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use pbc_lists
  use load_balancing

  implicit none

!  ARGUMENTS


  !-------------------shanghui add here------------------------------------
  real*8, target, dimension(n_full_points)        :: partition_tab_std
  real*8, target, dimension(3,n_atoms,n_full_points) :: rho_free_gradient
  real*8, target, dimension(3,n_atoms,n_full_points) :: v_free_gradient
  !------------------shanghui end add------------------------------------- 

!  INPUTS
! o partition_tab_std -- values of partition function
!
!  OUTPUT
! o rho_free_gradient -- gradient of total free atom density
! o v_free_gradient -- gradient of total free atom potential
!  AUTHOR
!    Honghui Shang @ FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
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


  ! Local variables

!------------------real calc-----------
  integer i_center
  integer i_batch
  integer i_coord
  integer i_index
  integer :: i_full_points
  integer :: current_spl_atom
  integer :: current_center
  integer  i_atom, j_atom


  real*8 coord_current(3)

!-------(1)---------
  real*8 dist_tab_sq
  real*8 dir_tab(3)

!-------(2)---------
  real*8 dist_tab_in
  real*8 i_r
  real*8 i_r_log
  real*8 dir_tab_in(3)

!-------(3)--------
  real*8 log_weight
  real*8 radial_weight


  real*8 :: d_v_hartree_free_d_r
  real*8 :: d_rho_free_d_r

  real*8 :: rho_free_gradient_temp(3)
  real*8 :: v_free_gradient_temp(3)

  real*8, pointer :: partition_tab(:)
!-----------------end real calc-----------


!-------------mpi---------------------------
  integer :: info
  character*200 :: info_str

  ! Load balancing stuff
  integer n_my_batches_work  ! Number of batches actually used
  integer n_full_points_work ! Number of integration points actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

  integer n_bp

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all
!-------------end mpi---------------------------


  if(use_batch_permutation > 0) then
    call localorb_info("Summing up the Hartree potential with load balancing.", use_unit,'(2X,A)', OL_norm )
  else
    call localorb_info("Summing up the Hartree potential.", use_unit,'(2X,A)', OL_norm )
  endif

  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)

  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then

    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    n_full_points_work = batch_perm(n_bp)%n_full_points

    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab

  else
    n_my_batches_work = n_my_batches
    n_full_points_work = n_full_points
    batches_work => batches
    partition_tab => partition_tab_std

  endif

  if(get_batch_weights) then
    allocate(batch_times(n_my_batches_work))
    batch_times(:) = 0
  endif

  !-----------------------------------------------------------------------------




  call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()


    do i_center = 1, n_centers_hartree_potential, 1

      current_center   = centers_hartree_potential(i_center)
      current_spl_atom = center_to_atom(current_center)

      ! Now follows the real work: Summing the multipole potential, density, and their
      ! derivatives on the integration grid

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
            !------(1) get dist_tab_sq,dir_tab-----
            call tab_single_atom_centered_coords_p0 &
                 ( current_center, &
                 coord_current,  &
                 dist_tab_sq,  &
                 dir_tab )


            rho_free_gradient_temp(1:3)=0.0d0
            v_free_gradient_temp(1:3)=0.0d0


            !-----near field-----
            !if (dist_tab_sq.lt.multipole_radius_sq(current_spl_atom) ) then
            

            !-----keep the same as free_rho_superpos-----
            if (dist_tab_sq.lt.multipole_radius_free_sq(species(current_spl_atom))) then

             !------(2) get dist_tab_in, i_r, i_r_log, dir_tab_in-----
             call tab_single_atom_centered_coords_radial_log_p0 &
                   ( current_center, dist_tab_sq, dir_tab,  &
                   dist_tab_in, i_r, i_r_log, dir_tab_in )

                ! need radial derivatives di/dr for the benefit of
                ! spline derivative evaluation later
             !------(3) get log_weight, radial_weight
             call tab_single_radial_weights_v2 &
                     ( current_spl_atom, dist_tab_in, i_r, &
                     log_weight, radial_weight )


                d_v_hartree_free_d_r =  &
                     val_spline_deriv( i_r_log, &
                     free_pot_es_spl(1,1,species(current_spl_atom)),  &
                     n_grid(species(current_spl_atom))) * log_weight  
                       
                d_rho_free_d_r = pi4_inv *  &
                     val_spline_deriv( i_r_log,  &
                     free_rho_spl(1,1,species(current_spl_atom)),  &
                     n_grid(species(current_spl_atom)))  * log_weight


                 rho_free_gradient_temp(:)=dir_tab_in(:) * d_rho_free_d_r
                 v_free_gradient_temp(:)=dir_tab_in(:) * d_v_hartree_free_d_r

             endif 
            !endif ! only near part



           do i_coord = 1, 3

           rho_free_gradient(i_coord,center_to_atom(current_center),i_full_points)= &
           rho_free_gradient(i_coord,center_to_atom(current_center),i_full_points)+ &
           rho_free_gradient_temp(i_coord)

           v_free_gradient(i_coord,center_to_atom(current_center),i_full_points)= &
           v_free_gradient(i_coord,center_to_atom(current_center),i_full_points)+ &
           v_free_gradient_temp(i_coord)
             
           enddo


          end if ! end if (partition_tab.gt.0.d0)
        end do  ! end loop over points in a batch
        if(get_batch_weights) batch_times(i_batch) = batch_times(i_batch) + mpi_wtime() - time_start
      end do ! end loop over batches
    end do  ! end loop over source atoms



         


  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_world,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for free_atom_sum_gradient: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)


  if(get_batch_weights) then
    call set_batch_weights(n_bp,batch_times)
    deallocate(batch_times)
  endif


  ! Write one last line to bound output
  write(info_str,*) ' '
  call localorb_info(info_str, use_unit,'(A)',OL_norm)



end subroutine integrate_free_atom_sum_gradient_p0
!******
!------------------------------------------------------------------------------------------------------------
