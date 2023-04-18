!****s* FHI-aims/get_angular_grid
!  NAME
!    get_angular_grid
!  SYNOPSIS

subroutine get_absorption &
     (KS_eigen, KS_vec, KS_vec_complex, occ_numbers, chemical_potential, &
      partition_tab, l_shell_max, ep1_in)

!  PURPOSE
!
!  Wrapper function for calculating and outputting of the  absorption
!
!  0. ep1_in from character to integer (x->1,y->2,z->3);
!     allocations 
!  1.function 'calculate_mommat_p0', real space integral of atom centered basis 
!    function i with the gradient of function j for in each cell (cell distance)
!  2.function 'construct_overlap', 'fourier transform' of real space integrals 
!    by summing up the phases resulting from cell distances
!    funtion is called once for the upper triangle and once for the lower
!    triangle of the sparse matrix from first step
!  3.function 'calc_moment_p0', basis transformation from atom centered basis to 
!    KS-orbitals
!  4.function 'calc_dielectric_function', momentum matrix elements are summed up
!    to yield the absorption
!  (5.) function 'out_absorption', output of result to file
!  USES
  use calculate_mommat_base
  use dimensions
  use runtime_choices
  use localorb_io
  use mpi_utilities
  use synchronize_mpi_basic, only: sync_vector, sync_real_number
  use pbc_lists
!  ARGUMENTS

  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) :: KS_eigen
  complex*16, dimension(n_basis, n_states, n_spin, n_k_points), INTENT(IN)::  &
                                                                 KS_vec_complex
  real*8, dimension(n_basis, n_states, n_spin, n_k_points), INTENT(IN)::  KS_vec
  real*8, dimension(n_states, n_spin,n_k_points), INTENT(IN) :: occ_numbers
  real*8, INTENT(IN) :: chemical_potential
  real*8, target, dimension(n_full_points), INTENT(IN) :: partition_tab
  integer, dimension(n_species), INTENT(IN) :: l_shell_max 
  CHARACTER(len=1), INTENT(IN):: ep1_in



!  INPUTS
!   o KS_eigen
!   o KS_vec/KS_vec_complex
!   o occ_numbers
!   o chemical_potential
!   o partition_tab
!   o l_shell_max
!  OUTPUT
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

  real*8, dimension(:), allocatable :: dielectric_function
  real*8:: omega_plasma

  character*128 :: info_str
  integer :: info
  integer :: n_state_min_in
  integer :: n_state_max_in
  logical :: use_absorption
  !  counters

  integer :: i_k
  integer :: new_k_point
  integer :: coord_1

  !  begin work
    write(info_str,'(6X,A,1X,I4)') "Absorption post processing starts"
    call localorb_info ( info_str )
    call get_state_minmax_K(KS_eigen, n_state_min_in, n_state_max_in)
    if (.not.allocated( dielectric_function)) then
      allocate( dielectric_function ( n_omega),stat=info)
      call check_allocation(info, 'dielectric_function')
    end if
    if (ep1_in=='x')then
      coord_1=1
    elseif(ep1_in=='y')then
      coord_1=2
    elseif(ep1_in=='z')then
      coord_1=3
    else
      coord_1=1
    endif
    call allocate_mommat()
    call calculate_mommat_p0 ( partition_tab, l_shell_max, &
                              mommat_full_oned_up, mommat_full_oned_low,coord_1)
    i_k = 0
    omega_plasma=0
    dielectric_function=0
    do new_k_point = 1,n_k_points
      if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points ) then
	! k-point stored on current mpi task
	i_k = i_k+1 !new_k_point
	call allocate_mommat_k()
	call allocate_moment_one(n_state_min_in, n_state_max_in)
	call construct_overlap( mommat_full_oned_up, mommat_full_w_up,&
                                        mommat_full_w_complex_up, new_k_point,&
                                        work_ovl_mom )
	call construct_overlap( mommat_full_oned_low, &
                     mommat_full_w_low, mommat_full_w_complex_low, new_k_point,&
                     work_ovl_mom )
	call calc_moment_p0(moment_one,mommat_full_w_up, &
                     mommat_full_w_low, mommat_full_w_complex_up, &
                     mommat_full_w_complex_low, KS_vec(:,:,:,i_k) ,&
                     KS_vec_complex(:,:,:,i_k) , new_k_point, coord_1, &
                     n_state_min_in, n_state_max_in)
	call calc_dielectric_function(moment_one,moment_one,&
                     dielectric_function,omega_plasma,chemical_potential,&
                     KS_eigen(:,:,new_k_point),k_weights(new_k_point),&
                     occ_numbers(:,:,new_k_point),.True.,widthone, &
                     widthtwo, n_state_min_in, n_state_max_in)
	call clean_mommat()
      endif
    enddo
    if(.not. use_local_index) call sync_vector(dielectric_function, n_omega )
    if(.not. use_local_index) call sync_real_number(omega_plasma)
    if (myid == 0) then
    call out_absorption(dielectric_function,ep1_in)
    endif           
    call clean_mommat_final()
    if (allocated(dielectric_function))  deallocate(dielectric_function)
    write(info_str,'(6X,A,1X,I4)')"Absorption post processing finished"

end subroutine get_absorption
!******	
