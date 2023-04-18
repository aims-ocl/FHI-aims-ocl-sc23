!****s* FHI-aims/output_cube_files_p2
!  NAME
!  output_cube_files_p2
!  SYNOPSIS
subroutine output_cube_files_p2_cpscf (first_order_density_matrix)
!TODO: Timing information
!TODO : Check uses in all subfunctions
  !  PURPOSE
  !  Plot charge density for periodic systems. 
  !  This routine is thought as an improvement
  !  over the old output_cube_files_p1
  !
  !  USES
  !
  use dimensions
  use plot
  use runtime_choices
  use physics
  use synchronize_mpi
  use mpi_tasks
  use mpi_utilities
  use timing
  use hartree_potential_recip !for potential output..
  use constants
  use free_atoms,   only:average_free_es_pot 
  use lpb_solver_utilities, only: dielecfunc_from_density_points, dielecfunc_gradient_from_density_gradient,&
    phi_zero_mpb, kBT_mpb,z_mpb,kappainf_mpb,kappafunc_from_density_points,dielecfunc_gradient_from_density_gradient_points,&
	rhomin_kappa,rhomax_kappa,rhomin_kappa_anion,rhomax_kappa_anion, use_separate_alpha
  use esp_grids,only:  n_full_points_esp, esp_get_n_compute_dens
  use esp_charges,only: multipole_expantion_rho,sum_up_potential,esp_collect_pot,&
                        get_vxc_esp
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  use species_data, only: l_shell_max
    !  INPUTS
  !   none
  !  OUTPUT
  !   Writes charge density to file
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

 implicit none
  real*8 first_order_density_matrix(n_hamiltonian_matrix_size)

  !local variables
  real*8,     dimension(:,:), allocatable :: densmat
  real*8,     dimension(:),   allocatable :: densmat_sparse
  real*8,     dimension(:),   allocatable :: KS_vec
  complex*16, dimension(:),   allocatable :: KS_vec_complex
  real*8, dimension(:),   allocatable :: KS_vec_complex_im
  real*8, dimension(:),   allocatable :: KS_vec_complex_re
  complex*16, dimension(:),   allocatable :: KS_vec_soc_up
  complex*16, dimension(:),   allocatable :: KS_vec_soc_dn

  real*8 :: cube_units(3) !Step size along each edge  
  real*8 :: offset_coord(3)
  character(LEN=15):: file_format
  character(LEN=50) :: stm_filename

  integer :: n_points !Number of points of the cube file
  !real*8, dimension(:,:), allocatable :: points !actual cartesian coordiantes of the grid points
  real*8, dimension(:,:), allocatable :: minicube_points !cartesian coordinates of the grid points of each minicube
  real*8 en_lower_limit, en_upper_limit

  real*8, dimension (:), allocatable :: cube_output !stores the output of the cube
  real*8, dimension (:), allocatable :: minicube_output !stores the output of the minicube
  real*8, dimension (:,:), allocatable :: minicube_kindens !kinetic on each minicube 
  complex*16, dimension (:), allocatable :: minicube_wavefunc_complex
  complex*16, dimension (:), allocatable :: minicube_wavefunc_complex_temp
  !real*8, dimension (:), allocatable :: minicube_tmp !tempoary information for the minicube 
  real*8, dimension (:), allocatable :: minicube_density !Density on each minicube 

  !SR: MPB solvation
  real*8, dimension (:), allocatable :: dielec_func  
  real*8, dimension (:), allocatable :: ion_dens, kappa_func,&
delta_v_cube,delta_v_gradient_abs,kappa_func_anion
  real*8 :: current_potential
  integer :: i_cube_point
  real*8, dimension (:), allocatable :: my_free_hartree_superpos, my_free_rho_superpos
  real*8, dimension(:,:), allocatable :: my_rho
  real*8, dimension (:,:,:), allocatable :: my_rho_gradient
  real*8, dimension(:,:), allocatable :: my_free_rho_gradient_superpos
  real*8, dimension (:,:), allocatable :: dielec_func_grad
  real*8, dimension (:,:), allocatable :: kappa_func_gradient,&
kappa_func_gradient_anion !Density on each minicube 
  real*8, dimension (:,:), allocatable :: minicube_density_gradient_tp !Density on each minicube 
  !end MPB solvation
  real*8, dimension (:,:), allocatable :: minicube_density_gradient !Density on each minicube 
  real*8, dimension (:), allocatable :: minicube_initial_density
  real*8, dimension (:), allocatable :: minicube_wavefunc
  real*8, dimension (:), allocatable :: minicube_D
  real*8, dimension (:), allocatable :: minicube_D0
  real*8, dimension (:), allocatable :: minicube_embedding_potential
  real*8, dimension (:,:), allocatable :: elf_dens ! for elf, need both spin-up and spin-down densities
  real*8, dimension (:,:,:), allocatable :: elf_dens_grad ! density gradient for both spin channels
  real*8, dimension (:,:), allocatable :: elf_kindens ! kinetic energy density for both spin channels
  real*8 :: tp1, tp2

  ! ESP infrastructure
  real*8, dimension(:),allocatable                       :: &
                                            delta_v_hartree_part_at_zero_trans
  real*8, dimension(:,:),allocatable                     :: &
                                        delta_v_hartree_deriv_l0_at_zero_trans
  real*8, dimension( :, :),allocatable :: multipole_moments_trans
  real*8, dimension(:),allocatable                     :: &
                                                     multipole_radius_sq_trans
  integer, dimension(:),allocatable                      :: &
                                              l_hartree_max_far_distance_trans
  real*8, dimension(:, :),allocatable          :: outer_potential_radius_trans
  real*8, target, dimension(:) ,allocatable        :: partition_tab_esp
  real*8, target, dimension(:,:) ,allocatable        :: potential_trans
  real*8, dimension(:),allocatable :: free_hartree_superpos_trans
  real*8, dimension(:) ,allocatable        :: radius_esp_min
  real*8, dimension(:) ,allocatable        :: radius_esp_max
  real*8, dimension(:,:), allocatable :: hamiltonian_two

  !helper
  integer :: info !Gets the return from allocate statements
  real*8 :: inv_bohr_3 
  real*8 :: sqrt_inv_bohr_3 
  character(LEN=100) :: info_str
  real*8 :: current_coord

  !Variables governing minicubes
  integer :: n_minicubes_x
  integer :: n_minicubes_y
  integer :: n_minicubes_z
  integer :: n_minicubes
  integer :: minicube_npoints
  integer :: n_minicubes_per_task
  integer :: i_minipoint
  integer :: max_divisor
  logical :: dumpme
  logical :: dumpcube

  
  !counters
  integer :: i_cube !Counter for the i-th Cube-File; the total number of requested cube files is given by n_cube
  integer :: i_edge1, i_edge2, i_edge3, i_edge
  integer :: i_point
  integer :: i_coord ! 1-3 .. x-z
  integer :: i_k_point
  integer :: i_state
  integer :: i_spin
  integer :: eig_dens_count
  integer :: i_minicube_x
  integer :: i_minicube_y
  integer :: i_minicube_z
  integer :: i_minicube
  integer :: i_minicube_point
  integer :: i_task

 integer :: blocksize = 0
  
  integer:: n_compute_basis_local, n_compute_basis
  integer :: n_compute, n_compute_a
  integer :: i_kpoint

  integer :: i_basis

  !For potential:
  real*8  dip_length, dip_origin, dip_gradient
  real*8, dimension(3) ::  dummy


  !More variables
  real*8, allocatable :: average(:)
  character*100 :: average_filename
  integer :: FirstPoint, LastPoint, points_in_plane, i_plane, i_x, i_y !Counters
  integer :: i_current

  character(*), parameter :: func = 'output_cube_files_p2'

  !begin work
  !check for units and issue a warning if legacy option is enabled
  i_current = 0

  write(info_str, '(1X, A)') ''
  call localorb_info(info_str)
  if(trim(cube_content_unit).eq.'legacy') then
    write(info_str, '(2X, A)') '* WARNING: The output  of  the cube file is in 1/Ang^3  or'
    call localorb_info(info_str)
    write(info_str, '(2X, A)') '* WARNING: 1/Ang^3/2,  but  the cube voxel  definition  is'
    call localorb_info(info_str)
    write(info_str, '(2X, A)') '* WARNING: in atomic units. Be careful when postprocessing' 
    call localorb_info(info_str)
    write(info_str, '(2X, A)') '* WARNING: Consistent units can be chosen in control.in by'
    call localorb_info(info_str)
    write(info_str, '(2X, A)') '*               cube_content_unit  bohr                 '
    call localorb_info(info_str)
    write(info_str, '(2X, A)') ''
    call localorb_info(info_str)
    warn_cube_file_format = .true.
  endif

   ! Set things that should be constant
   !FIXME: This is not exactly clean programming, but changing this once instead 
   !FIXME: of writing half a dozen of if-statements save time and makes the code
   !FIXME: easier to read, in my opinion. Feel free to change.
   if(cube_content_unit.eq.'legacy') then
        inv_bohr_3 = 1.0d0/(bohr**3)             
        sqrt_inv_bohr_3 = ( inv_bohr_3 )**(0.5d0)
   else
        inv_bohr_3 = 1.0d0
         sqrt_inv_bohr_3 = 1.0d0
   endif
 
   do i_cube=1,n_cube,1 !Iterate over all cubes
        

   !Defaults
   en_upper_limit= 1d10
   en_lower_limit=-1d10
   n_points = (cube_edge_steps(1,i_cube))*(cube_edge_steps(2,i_cube))&
             *(cube_edge_steps(3,i_cube))
  !Safeguard against very, very small cube files
   cube_divisor(i_cube)=min(min(min(cube_edge_steps(1,i_cube),cube_edge_steps(2,i_cube)),& 
                                    cube_edge_steps(3,i_cube)),cube_divisor(i_cube))
!Simplifications and allocation of minicube-arrays
   call aims_allocate( minicube_output, cube_divisor(i_cube)**3,    "minicube_output" )
   call aims_allocate( minicube_points, 3, cube_divisor(i_cube)**3, "minicube_points" )
   
   select case (cube_type(i_cube))
   case('total_density')
           allocate (minicube_density(cube_divisor(i_cube)**3),stat=info) !allocate density
           call check_allocation(info, 'minicube_density', func)
           minicube_density=0.0
   case('delta_density')
           allocate (minicube_density(cube_divisor(i_cube)**3),stat=info) !allocate density
           call check_allocation(info, 'minicube_density', func)
           allocate (minicube_initial_density(cube_divisor(i_cube)**3),stat=info) !allocate density
           call check_allocation(info, 'minicube_initial_density', func)
           minicube_density=0.0	   
   case('delta_v')
           allocate (delta_v_cube(cube_divisor(i_cube)**3),stat=info)
           call check_allocation(info, 'delta_v_cube', func)
           allocate (delta_v_gradient_abs(cube_divisor(i_cube)**3),stat=info)
           call check_allocation(info, 'delta_v_gradient_abs', func)   
   case('ion_dens')
           allocate (minicube_density(cube_divisor(i_cube)**3),stat=info) !allocate density
           call check_allocation(info, 'minicube_density', func)   
           minicube_density=0.0	
           allocate  (minicube_density_gradient(cube_divisor(i_cube)**3,3),stat=info) ! allocate minicube_gradient
           call check_allocation(info,'minicube_density_gradient', func)
           minicube_density_gradient=0.0
           allocate  (minicube_density_gradient_tp(cube_divisor(i_cube)**3,3),stat=info) ! allocate minicube_gradient
           call check_allocation(info,'minicube_density_gradient_tp', func)
           minicube_density_gradient_tp=0.0
           allocate (ion_dens(cube_divisor(i_cube)**3),stat=info)
           call check_allocation(info, 'ion_dens', func)
           allocate (delta_v_cube(cube_divisor(i_cube)**3),stat=info)
           call check_allocation(info, 'delta_v_cube', func)
           allocate (delta_v_gradient_abs(cube_divisor(i_cube)**3),stat=info)
           call check_allocation(info, 'delta_v_gradient_abs', func)
           ion_dens = 0d0
	   allocate (dielec_func(cube_divisor(i_cube)**3),stat=info)
           call check_allocation(info, 'dielec_func', func)
	   allocate (kappa_func(cube_divisor(i_cube)**3),stat=info)
           call check_allocation(info, 'kappa_func', func)
	   if (use_separate_alpha) then
		allocate (kappa_func_anion(cube_divisor(i_cube)**3),stat=info)
			   call check_allocation(info, 'kappa_func_anion', func)
	   allocate (kappa_func_gradient_anion(3,cube_divisor(i_cube)**3),stat=info)
           call check_allocation(info, 'kappa_func_gradient_anion', func)           
	   end if
	   allocate (kappa_func_gradient(3,cube_divisor(i_cube)**3),stat=info)
           call check_allocation(info, 'kappa_func_gradient', func)           
	   allocate (my_free_hartree_superpos(cube_divisor(i_cube)**3),stat=info)
	   call check_allocation(info, 'my_free_hartree_superpos', func)
	   allocate (my_rho(n_spin,cube_divisor(i_cube)**3),stat=info)
	   call check_allocation(info, 'my_rho', func)	 
	   allocate (my_free_rho_superpos(cube_divisor(i_cube)**3),stat=info)
	   call check_allocation(info, 'my_free_rho_superpos', func)	  
	   allocate (my_rho_gradient(3,n_spin,cube_divisor(i_cube)**3),stat=info)
	   call check_allocation(info, 'my_rho_gradient', func)	
	   allocate (my_free_rho_gradient_superpos(3,cube_divisor(i_cube)**3),stat=info)
	   call check_allocation(info, 'my_free_rho_gradient_superpos', func)	
           allocate (dielec_func_grad(3,cube_divisor(i_cube)**3),stat=info) !allocate density
           call check_allocation(info, 'dielec_func_grad', func)

   case('dielec_func')
!            allocate (minicube_initial_density(cube_divisor(i_cube)**3),stat=info) !allocate density
!            call check_allocation(info, 'minicube_initial_density', func)
           allocate (dielec_func(cube_divisor(i_cube)**3),stat=info) !allocate density
           call check_allocation(info, 'dielec_func', func)
	   dielec_func=0.0
!            allocate (dielec_func_grad(3,cube_divisor(i_cube)**3),stat=info) !allocate density
!            call check_allocation(info, 'dielec_func_grad', func)
! 	   dielec_func_grad=0.0
!            allocate  (minicube_density_gradient(cube_divisor(i_cube)**3,3),stat=info) ! fallocate minicube_gradient
!            call check_allocation(info,'minicube_density_gradient', func)
!            minicube_density_gradient=0.0
           allocate (minicube_density(cube_divisor(i_cube)**3),stat=info) !allocate density
           call check_allocation(info, 'minicube_density', func)
!            allocate (minicube_initial_density(cube_divisor(i_cube)**3),stat=info) !allocate density
!            call check_allocation(info, 'minicube_initial_density', func)
           minicube_density=0.0	              
   case('stm')
     if(cube_stm(3,i_cube).gt.0) then 
           allocate (minicube_density(cube_divisor(i_cube)**3),stat=info) !allocate density
           call check_allocation(info, 'minicube_density', func)
           minicube_density=0.0
           cube_stm(1,i_cube)=chemical_potential
           if(cube_stm(2,i_cube).le.0)then
                  en_upper_limit=cube_stm(1,i_cube)
                  en_lower_limit=cube_stm(1,i_cube)+cube_stm(2,i_cube)
           endif
           if(cube_stm(2,i_cube).ge.0)then
                  en_upper_limit=cube_stm(1,i_cube)+cube_stm(2,i_cube)
                  en_lower_limit=cube_stm(1,i_cube)
            endif
     endif
   case('spin_density')  !spin density is simply the difference of up and down
           allocate (minicube_density(cube_divisor(i_cube)**3),stat=info) !allocate density
           call check_allocation(info, 'minicube_density', func)
           minicube_density=0.0
           cube_type(i_cube)='total_density'
           cube_spin(1,i_cube)=1
           cube_spin(2,i_cube)=-1
    case('elf')
           allocate (elf_dens(cube_divisor(i_cube)**3,n_spin),stat=info)
           call check_allocation(info, 'elf_dens', func)
           allocate (elf_dens_grad(3,cube_divisor(i_cube)**3,n_spin),stat=info)
           call check_allocation(info, 'elf_dens_grad', func)
           allocate (elf_kindens(cube_divisor(i_cube)**3,n_spin),stat=info)
           call check_allocation(info, 'elf_kindens', func)
    case('eigenstate')
        if (real_eigenvectors.or.out_cube_soc) then
            call aims_allocate(minicube_wavefunc, cube_divisor(i_cube)**3, "minicube_wavefunc") ! allocate wave function
            minicube_wavefunc=0.0
        else
            call aims_allocate(minicube_wavefunc_complex, &
                 cube_divisor(i_cube)**3, "minicube_wavefunc_complex") ! allocate wave function
            minicube_wavefunc_complex=0.0
        endif
! eigenstate of non-soc when soc triggered
    case('eigenstate_non_soc')
        if (real_eigenvectors) then
            call aims_allocate(minicube_wavefunc, cube_divisor(i_cube)**3, "minicube_wavefunc") ! allocate wave function
            minicube_wavefunc=0.0
        else
            call aims_allocate(minicube_wavefunc_complex, &
                 cube_divisor(i_cube)**3, "minicube_wavefunc_complex") ! allocate wave function
            minicube_wavefunc_complex=0.0
        endif     
    case('eigenstate_imag')
        if (real_eigenvectors.and..not.out_cube_soc) &
             call aims_stop('Output of imiaginary part of eigenstate for system with real eigenvectors requested')
        call aims_allocate ( minicube_wavefunc, cube_divisor(i_cube)**3,"minicube_wavefunc" )
        minicube_wavefunc=0.0
!        call aims_allocate( minicube_wavefunc_complex, cube_divisor(i_cube)**3, "minicube_wavefunc_complex" )
!        minicube_wavefunc_complex=0.0
    case('eigenstate_density')
        if (real_eigenvectors.and..not.out_cube_soc) then
            call aims_allocate ( minicube_wavefunc, cube_divisor(i_cube)**3, "minicube_wavefunc" )
            minicube_wavefunc=0.0
        else
            call aims_allocate( minicube_wavefunc_complex, cube_divisor(i_cube)**3, "minicube_wavefunc_complex" )
            minicube_wavefunc_complex = (0.0d0,0.0d0)
            if (out_cube_soc) then
                 call aims_allocate( minicube_wavefunc_complex_temp, cube_divisor(i_cube)**3, "minicube_wavefunc_complex_temp" )
                 minicube_wavefunc_complex_temp = (0.0d0,0.0d0)
            end if
        endif
    case('long_range_potential')
            allocate(minicube_embedding_potential(cube_divisor(i_cube)**3),stat=info)
            call check_allocation(info,'mincube_embedding_potential')
            minicube_embedding_potential=0.d0
    case('potential')
        continue
    case('hartree_potential')
        continue
    case('xc_potential')
        continue
    case default
        call aims_stop('Unkown cube type in Simplification block')
    end select
   
     !allocations of general arrays
     if (any(cube_type_needs_densmat)) then
          ! WPH: These matrices should be the part of the memory bottleneck, unless the user picks a comically
          !      dense cube file or there are unnecessarily large work matrices in called subroutines.  (There
          !      will likely be work matrices needed with similar dimensions, at least in the LAPACK case.)

          ! Incompatible options should have already been spotted when reading control.in, but to be on the safe
          ! side, do it again
          if (out_cube_soc) then
              call localorb_info("")
              call aims_stop_coll( "* You have selected a cube output that requires the density matrix be computed. &
                                    & This is currently not supported when using spin-orbit coupling.  Exiting.", func )
          end if
          if (use_local_index) then
              call localorb_info("")
              call aims_stop_coll( "* You have selected a cube output that requires the density matrix be computed. &
                                    & This is currently not supported with the use_local_index keyword.  Exiting.", func )
          end if
          if (use_load_balancing) then
              call localorb_info("")
              call aims_stop_coll( "* You have selected a cube output that requires the density matrix be computed. &
                                    & This is currently not supported with the load_balancing keyword.  Exiting.", func )
          end if
          if (packed_matrix_format == PM_none) then
              allocate(densmat(n_centers_basis_T,n_centers_basis_T),stat=info) 
              call check_allocation(info,'densmat',func)
              allocate(densmat_sparse(1),stat=info) !Dummy allocation
              call check_allocation(info,'densmat_sparse',func)
              densmat(1:n_centers_basis_T,1:n_centers_basis_T)=0.0 
              densmat_sparse(1)=0.0
          elseif (packed_matrix_format == PM_index) then
              allocate(densmat(1,1),stat=info) !Dummy allocation
              call check_allocation(info,'densmat',func)
              allocate(densmat_sparse(n_hamiltonian_matrix_size),stat=info) 
              call check_allocation(info,'denmat_sparse',func)
              densmat(1,1)=0.0 
              densmat_sparse(1:n_hamiltonian_matrix_size)=0.0
         endif
    end if
    if (any(cube_type_needs_eigenvectors)) then
         ! WPH: If statement not necessary, because these arrays are relatively
         !      small and could be blanket allocated, but I'm a sucker for symmetry.
         if (out_cube_soc .and. cube_type_needs_soc_eigenvectors(i_cube)) then
! out_cube_soc and we really plot the soc eigenvectors
                   allocate(KS_vec_complex(n_basis),stat=info)
                   call check_allocation(info, 'KS_vec_complex', func)
                   KS_vec_complex(1:n_basis)=0.0
!              allocate( KS_vec_soc_up(n_basis),stat=info)
!              call check_allocation(info, 'KS_vec_soc_up', func)
!              KS_vec_soc_up(1:n_basis)=0.0
!              allocate( KS_vec_soc_dn(n_basis),stat=info)
!              call check_allocation(info, 'KS_vec_soc_dn', func)
!              KS_vec_soc_dn(1:n_basis)=0.0
         else 
              if (real_eigenvectors) then
                   allocate( KS_vec(n_basis),stat=info)
                   call check_allocation(info, 'KS_vec', func)
                   KS_vec(1:n_basis)=0.0
              else 
                   allocate(KS_vec_complex(n_basis),stat=info)
                   call check_allocation(info, 'KS_vec_complex', func)
                   KS_vec_complex(1:n_basis)=0.0
              end if
         endif
    endif
    !General allocations - 
    allocate (cube_output(n_points), stat=info) !allocate output array
    call check_allocation(info,'cube_output',func)
    cube_output = 0.0
    !allocate (points(3,n_points),stat=info) !allocate array
    !call check_allocation(info, 'points', func)

    !Minicube_allocations - allocate all minicubes just to be sure.
!    allocate(minicube_output(cube_divisor(i_cube)**3),stat=info)
!    call check_allocation(info, 'minicube_output', func)
!    allocate(minicube_tmp(cube_divisor(i_cube)**3),stat=info)
!    call check_allocation(info, 'minicube_tmp', func)
!    allocate (minicube_density(cube_divisor(i_cube)**3),stat=info) !allocate density
!    call check_allocation(info, 'minicube_density', func)
!    allocate  (minicube_wavefunc(cube_divisor(i_cube)**3),stat=info) ! allocate wave function
!    call check_allocation(info,'minicube_wavefunc', func)
!    allocate  (minicube_density_gradient(cube_divisor(i_cube)**3,3),stat=info) ! allocate minicube_gradient
!    call check_allocation(info,'minicube_density_gradient', func)
!    allocate  (minicube_wavefunc_complex(cube_divisor(i_cube)**3),stat=info) ! allocate wave function
!    call check_allocation(info,'minicube_wavefunc_complex', func)
!    allocate (minicube_points(3,cube_divisor(i_cube)**3),stat=info) !allocate array
!    call check_allocation(info, 'minicube_points', func)
!    allocate (minicube_initial_density(cube_divisor(i_cube)**3),stat=info) !allocate density
!    call check_allocation(info, 'minicube_initial_density', func)
!    allocate (minicube_kindens(2, cube_divisor(i_cube)**3),stat=info) !allocate density
!    call check_allocation(info, 'minicube_kindens', func)


    !Basic preparations
    if(cube_type_needs_densmat(i_cube)) then
        if (packed_matrix_format==PM_none)  call prepare_densmat(i_cube,en_lower_limit,en_upper_limit,densmat,densmat_sparse)
        if (packed_matrix_format==PM_index) call prepare_densmat_sparse(i_cube,en_lower_limit,en_upper_limit,densmat,densmat_sparse)
    endif

    densmat_sparse = first_order_density_matrix

    if (cube_type_needs_eigenvectors(i_cube)) then
        if(out_cube_soc .and. cube_type_needs_soc_eigenvectors(i_cube)) then
             call prepare_soc_eigenvector( i_cube, KS_vec_complex)
        else 
             if (real_eigenvectors)      call prepare_eigenvector(i_cube,KS_vec)
             if (.not.real_eigenvectors) call prepare_eigenvector_complex(i_cube,KS_vec_complex)
        end if
    end if

    if (n_periodic.eq.3) then
         !prepate basic variables
          dip_length =  abs((maxval(lattice_vector(3,:)) - minval(lattice_vector(3,:))))
          dip_origin = vacuum_z_level + dip_length * 0.5d0
          dip_gradient = -(pot_jump/dip_length)/0.5d0 
    endif

!   Initializations, just to be sure
!   minicube_points(1:3,cube_divisor(i_cube)**3)=0.0
!   minicube_wavefunc_complex(1:cube_divisor(i_cube)**3) = 0.0
!   minicube_kindens(1:2,1:cube_divisor(i_cube)**3)=0.0
!   !minicube_tmp(1:cube_divisor(i_cube)**3)=0.0
!   points(1:3,1:n_points)=0.0
!   cube_output(1:n_points)= 0.0
!   minicube_density(i_minicube_point)=0.0
!   minicube_density_gradient(i_minicube_point,1:3)=0.0
!   minicube_initial_density(i_minicube_point)=0.0
!   minicube_wavefunc(i_minicube_point) = 0.0


     !Determine in how many batches ("minicubes") the cube file will be divided
     n_minicubes_x =  int(ceiling(cube_edge_steps(1,i_cube)/dble(cube_divisor(i_cube))))
     n_minicubes_y =  int(ceiling(cube_edge_steps(2,i_cube)/dble(cube_divisor(i_cube))))
     n_minicubes_z =  int(ceiling(cube_edge_steps(3,i_cube)/dble(cube_divisor(i_cube))))
     n_minicubes =  n_minicubes_x*n_minicubes_y*n_minicubes_z  
     if (myid.eq.0) then
        write(use_unit,*) ''
        write(info_str,*) '  Starting cube ', i_cube
        call localorb_info(info_str)
        call localorb_info('   Cube will be divided in ')
        write(info_str,fmt='(A,I0,A)') '     ', n_minicubes_x,  ' minicubes along the first, ' 
        call localorb_info(info_str)
        write(info_str,fmt='(A,I0,A)') '     ', n_minicubes_y,  ' minicubes along the second, and ' 
        call localorb_info(info_str)
        write(info_str,fmt='(A,I0,A)') '     ', n_minicubes_z,  ' minicubes along the third cube edge ' 
        call localorb_info(info_str)
    endif
  

   !Work out and store all the points
   !Calculate points
   !Origin must be in center - OTH update: We want to calculate everything as if
   !the point was in the central unit cell. However, we still want to allow the
   !user to choose the origin of his output file whereever he wants. Hence, the
   !following mapping should not be done. Rather, the map to center cell occurs
   !in the subroutine where the coordinates of each cube voxel is actually
   !calculated
   !if (n_periodic.gt.0)  call map_to_center_cell(cube_origin(1:3,i_cube))

   do i_coord = 1,3,1 !calculate step size
      cube_units(i_coord) = cube_edge_unit(i_coord,1,i_cube)+ &
                            cube_edge_unit(i_coord,2,i_cube)+ &
                            cube_edge_unit(i_coord,3,i_cube)
      offset_coord = real(cube_edge_steps(1:3,i_cube)-1)/2.d0* cube_units
   enddo

   !OTH: I commented this out to get rid of the points array, which 
   !     is stored on each CPU and can become huge for many systems.
   !     Instead, the coordinate of the hypothetical point is now
   !     calculated on the fly for each minicube in a seperate 
   !     subroutine.
   !     JWs concerns still apply.
   ! i_point = 1
   ! call localorb_info('   Setting up cube grid points')
   ! do i_edge1=1,(cube_edge_steps(1,i_cube)),1
   !      do i_edge2=1,(cube_edge_steps(2,i_cube)),1
   !         do i_edge3=1,(cube_edge_steps(3,i_cube)),1
   !          do i_coord = 1,3,1 !for all 3 spatial coordiantes
   !             points(i_coord,i_point)=(i_edge1-1)*cube_edge_unit(i_coord,1,i_cube)+ &
   !                               (i_edge2-1)*cube_edge_unit(i_coord,2,i_cube)+ &
   !                               (i_edge3-1)*cube_edge_unit(i_coord,3,i_cube)  &
   !                                 -offset_coord(i_coord)+cube_origin(i_coord,i_cube)
   !          enddo
             ! JW: FIXME: Two issues:
             ! (1) For eigenstate output (densities are fine), mapping back to
             !     the central unit cell means that the phase factor is going
             !     to be wrong.  We should either avoid to map back (e.g. by
             !     providing our own centers list in evaluate_wavefunc_*()) or
             !     record the corresponding cell_indices to reapply the needed
             !     phase.

             ! (2) For densities (everything using
             !     evaluate_rho_from_densmat()), mapping back might be
             !     extremely expensive in terms of memory.  We prune the
             !     Cbasis functions for all n_compute functions which are
             !     nonzero at least on /one/ grid point.  Mapping might end up
             !     having up to eight distant clusters of points within one
             !     minicube.  As the pruned density matrix is of size
             !     n_compute**2, mapping might cause an increase of memory
             !     footprint (and evaluation time, though less critical) of
             !     64.
   !          if (n_periodic.gt.0) then
   !                call map_to_center_cell(points(1:3,i_point))
   !          endif
   !          i_point=i_point+1 !Proceed to next point
   !        enddo
   !     enddo
   ! enddo

    !Now that we have all the points, we can distribute the minicubes..
if ((cube_type(i_cube).ne.'potential').and.(cube_type(i_cube).ne.'hartree_potential').and.(cube_type(i_cube).ne.'xc_potential')) then !...unless we want to calculate the potential, in which case the philosophy of the two routines conflicts

   !Find out how many minicubes per task
   if (n_tasks.gt.1) call localorb_info('   Dividing grid points on processors')
   if (n_tasks==1)   call localorb_info('   Calculating all minicubes on main processor')
   n_minicubes_per_task = ceiling(n_minicubes/dble(n_tasks))

   do  i_task = 0,(n_minicubes_per_task-1),1
          dumpcube =.false.
          i_minicube_x=modulo(modulo((myid+i_task*n_tasks),(n_minicubes_x*n_minicubes_y)),n_minicubes_x)+1
          i_minicube_y=floor(modulo((myid+i_task*n_tasks),(n_minicubes_x*n_minicubes_y))/dble(n_minicubes_x))+1
          i_minicube_z=floor((myid+i_task*n_tasks)/dble(n_minicubes_x*n_minicubes_y))+1
	!status output
 	if (myid==0 .and. int((1d0*i_task)/(n_minicubes_per_task-1)*100)==i_current .and.&
                solvent_method .eq. SOLVENT_MPB) then
 	! since the output can get quite slow, I add here a status output
	  write(use_unit,*) '>>> ', i_current, '%'
	  i_current = i_current + 1
	end if
   
   !Safeguard: Inelegant, but 'cycle' does not seem to work, nor a if-statement
   !which does nothing over the whole do-loop. Both break mpi runs for some reason. 
   !So this is the best I could come up with. At least it works and doesn't cost time anyway. 
        if (i_minicube_x.gt.n_minicubes_x) then
               i_minicube_x=n_minicubes_x
               dumpcube =.true.
        endif
        if (i_minicube_y.gt.n_minicubes_y) then
                i_minicube_y=n_minicubes_y
                dumpcube=.true.
        endif
        if (i_minicube_z.gt.n_minicubes_z) then
                  dumpcube =.true.
                  i_minicube_z=n_minicubes_z
        endif
   
   
      !Assign the global points to minicubes
        do i_minicube_point = 1,cube_divisor(i_cube)**3, 1
                 !This rather complicated-looking block makes sure that all points of the
                 !minicube are close together (i.e., in a cubic block)
                 i_point = floor((i_minicube_point-1)/dble(cube_divisor(i_cube)**2))*& 
                           cube_edge_steps(1,i_cube)*cube_edge_steps(2,i_cube)+ &
                           floor((modulo(i_minicube_point-1,cube_divisor(i_cube)**2))/& 
                           dble(cube_divisor(i_cube)))*cube_edge_steps(1,i_cube) + &
                           modulo(((modulo(i_minicube_point-1,cube_divisor(i_cube)**2)+1)-1),cube_divisor(i_cube))+1  + &
                            (i_minicube_x-1)*cube_divisor(i_cube) + & !offset in x
                            (i_minicube_y-1)*cube_divisor(i_cube)*cube_edge_steps(1,i_cube) + & !offset in y
                            (i_minicube_z-1)*cube_divisor(i_cube)*cube_edge_steps(1,i_cube)*cube_edge_steps(2,i_cube) !offset in z
               if (i_point.gt.n_points) i_point = n_points !ugly but efficient safeguard against points outside of array
               !minicube_points(1:3,i_minicube_point)=points(1:3,i_point)
               call obtain_coordinate_of_point(i_point,i_cube,minicube_points(1:3,i_minicube_point))
               if (n_periodic.gt.0)  call map_to_center_cell(minicube_points(1:3,i_minicube_point))
        enddo !i_minicube_point
   

   !Data mining - the core of the routine
   select case (cube_type(i_cube))
   case ('total_density') 
        if (packed_matrix_format==PM_none) then 
            call evaluate_rho_from_densmat(cube_divisor(i_cube)**3, minicube_points, .false., densmat, minicube_density)
        else 
            call evaluate_rho_from_densmat(cube_divisor(i_cube)**3, minicube_points, .true., densmat_sparse,  minicube_density)
        endif
        minicube_output=minicube_density*inv_bohr_3
   case('stm') 
        if (packed_matrix_format==PM_none) then                       
            call evaluate_rho_from_densmat(cube_divisor(i_cube)**3, minicube_points, .false., densmat, minicube_density)
        else            
            call evaluate_rho_from_densmat(cube_divisor(i_cube)**3, minicube_points, .true., densmat_sparse,  minicube_density)
        endif
        minicube_output=minicube_density*hartree*cube_stm(2,i_cube)*inv_bohr_3
   case('delta_density') 
        if (packed_matrix_format==PM_none) then
            call evaluate_rho_from_densmat(cube_divisor(i_cube)**3, minicube_points, .false., densmat, minicube_density)
        else
            call evaluate_rho_from_densmat(cube_divisor(i_cube)**3, minicube_points, .true., densmat_sparse,  minicube_density)
        endif
        call output_cube_initial_density(cube_divisor(i_cube)**3,minicube_points,minicube_initial_density)
        minicube_output=(minicube_density-minicube_initial_density)*inv_bohr_3
    case('delta_v') 
    !attention this method works only for cluster. with and without MPBE calculation
    !no far field terms implemented, yet and no periodic conditions. no guarantee
    !for correct results. tested for cluster systems.
    !for output of spline coefficients, needed for coupling to FEM code KARDOS
    !call write_free_spline_coeff()    
        call evaluate_potential_cube(cube_divisor(i_cube)**3, minicube_points, .false., delta_v_cube, delta_v_gradient_abs)
        minicube_output = delta_v_cube !delta_v_gradient_abs
    case('ion_dens') 
!	the same valid as for delta_v. no far field terms, yet and as also the MPBE solver only for non-periodic systems
! 	1.) evaluate all free atom properties, we need only my_free_hartree_superpos (vfree)
	call  initialize_cube_grid_storage( &
	      cube_divisor(i_cube)**3, minicube_points, &
	      my_free_hartree_superpos, my_free_rho_superpos,  &
	      my_free_rho_gradient_superpos, &
	      my_rho, my_rho_gradient)   
! 	2.) evaluate rho
!          !Calculate spin density
	if (packed_matrix_format==PM_none) then
	      call evaluate_rho_from_densmat(cube_divisor(i_cube)**3, minicube_points, .false., densmat, minicube_density)
	else
	      call evaluate_rho_from_densmat(cube_divisor(i_cube)**3, minicube_points, .true., densmat_sparse,  minicube_density)
	end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	if someone needs this, this includes the calculation of density gradients, dielectric function and its gradient
!	for checks
! !         !Calculate spin density gradient 
! !         !FIXME: This routine uses finite differences and is therefore really, really slow. 
! !         !It should be straightforward to just adjust evaluate_rho_from_densmat to work on the gradient instead.
!         if (packed_matrix_format==PM_none) then
!             call evaluate_rho_gradient_from_densmat(cube_divisor(i_cube)**3, &
!                     minicube_points, .false., densmat, &
!                     minicube_density_gradient)
!         else
!             call evaluate_rho_gradient_from_densmat(cube_divisor(i_cube)**3, &
!                     minicube_points, .true., densmat_sparse, &
!                     minicube_density_gradient)
!         endif
! !         !evaluate eps on cubic grid
! !        call output_cube_initial_density(cube_divisor(i_cube)**3,minicube_points,minicube_initial_density)
!         call dielecfunc_from_density_points( minicube_density, cube_divisor(i_cube)**3, dielec_func  )!         !evaluate grad eps (not necessary, only to have it if one wants to see it)
!           call dielecfunc_gradient_from_density_gradient_points(dielec_func, dielec_func_grad, minicube_density,&
!   	  minicube_density_gradient, cube_divisor(i_cube)**3)
! !         !evaluate alpha function
!          do i_coord=1, 3, 1
! 	  minicube_density_gradient_tp(i_coord,:) = minicube_density_gradient(:,i_coord)
! 	 end do
	minicube_density_gradient_tp = 0d0 !we do not need this, if needed uncomment the above lines
	dielec_func = 1d0 ! we need this only for alphakind = 0 in kappafunc_from_density but usually alphakind = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	3.)    	evaluate delta_v from multipole moments on cubic grid
	  !evaluate delta_v potential
        call evaluate_potential_cube( cube_divisor(i_cube)**3, minicube_points, .false., delta_v_cube, delta_v_gradient_abs  )

! 	4.)	evaluate kappa function
	call kappafunc_from_density_points( minicube_density, minicube_density_gradient_tp, &
	  dielec_func, cube_divisor(i_cube)**3, kappa_func,kappa_func_gradient,rhomin_kappa,rhomax_kappa) 
	if (use_separate_alpha) then
		call kappafunc_from_density_points( minicube_density, minicube_density_gradient_tp, &
		dielec_func, cube_divisor(i_cube)**3, kappa_func_anion,kappa_func_gradient_anion,rhomin_kappa_anion,rhomax_kappa_anion) 
	end if
! 	5.) 	evaluate ion density
	ion_dens = 0d0
	do i_cube_point = 1, cube_divisor(i_cube)**3, 1
	  if (kappa_func(i_cube_point) > 0d0) then
 	    current_potential = delta_v_cube(i_cube_point) + my_free_hartree_superpos(i_cube_point)
		if (.not.use_separate_alpha) then
			ion_dens(i_cube_point) =  pi4_inv*kappa_func(i_cube_point)/kappainf_mpb*kappainf_mpb**2*&
		      sinh(z_mpb*current_potential/kBT_mpb)/(1d0-phi_zero_mpb+phi_zero_mpb*&
		      kappa_func(i_cube_point)/kappainf_mpb*cosh(z_mpb*current_potential/kBT_mpb))
		else
			ion_dens(i_cube_point) =  pi4_inv/2d0*kappainf_mpb**2*1d0/kappainf_mpb*&
		      (kappa_func(i_cube_point)*exp(-z_mpb*current_potential/kBT_mpb)/(1d0-phi_zero_mpb+phi_zero_mpb*kappa_func(i_cube_point)/kappainf_mpb*cosh(z_mpb*current_potential/kBT_mpb))+&
			  kappa_func_anion(i_cube_point)*exp(z_mpb*current_potential/kBT_mpb)/(1d0-phi_zero_mpb+phi_zero_mpb*kappa_func_anion(i_cube_point)/kappainf_mpb*cosh(z_mpb*current_potential/kBT_mpb)))
	    end if
      end if
	end do
        minicube_output = ion_dens 
   case('dielec_func') 
	
        if (packed_matrix_format==PM_none) then
            call evaluate_rho_from_densmat(cube_divisor(i_cube)**3, minicube_points, .false., densmat, minicube_density)
        else
            call evaluate_rho_from_densmat(cube_divisor(i_cube)**3, minicube_points, .true., densmat_sparse,  minicube_density)
        endif
!        call output_cube_initial_density(cube_divisor(i_cube)**3,minicube_points,minicube_initial_density)
        call dielecfunc_from_density_points( minicube_density, cube_divisor(i_cube)**3, dielec_func  )
! 	write(use_unit,*) 'shape', shape(minicube_density_gradient)
! 	call evaluate_rho_gradient_from_densmat(cube_divisor(i_cube)**3, minicube_points, .true., densmat_sparse,  minicube_density_gradient)
! 	call dielecfunc_gradient_from_density_gradient(dielec_func, dielec_func_grad, minicube_initial_density,&
! 	  minicube_density_gradient,cube_divisor(i_cube)**3)

! 	do i_point = 1, cube_divisor(i_cube)**3, 1
! 	  minicube_output(i_point) = sqrt(sum(dielec_func_grad(:,i_point)**2))
! 	  write(use_unit,*) minicube_output(i_point)
! 	end do

	minicube_output = dielec_func
	
   case('elf')
      call evaluate_elf(i_cube, cube_divisor, minicube_points, elf_dens, elf_dens_grad, elf_kindens)
      if (n_spin.eq.1) then
         ! Savin et al. version
         do i_point = 1, cube_divisor(i_cube)**3
            tp1 = 0.3*pisq3**(2*third)*(elf_dens(i_point,1))**(8*third)
            tp2 = elf_kindens(i_point,1)*elf_dens(i_point,1) - &
                 0.125*dot_product(elf_dens_grad(:,i_point,1),elf_dens_grad(:,i_point,1))
            minicube_output(i_point) = tp1**2/(tp1**2+tp2**2+1d-6)
         enddo! i_point
      else ! spin-polarized case
         if (cube_elf(i_cube).eq.0) then
            ! Savin et al. version
            elf_dens(:,1) = elf_dens(:,1) + elf_dens(:,2)
            elf_kindens(:,1) = elf_kindens(:,1) + elf_kindens(:,2)
            elf_dens_grad(:,:,1) = elf_dens_grad(:,:,1) + elf_dens_grad(:,:,2)
            do i_point = 1, cube_divisor(i_cube)**3
               tp1 = 0.3*pisq3**(2*third)*(elf_dens(i_point,1))**(8*third)
               tp2 = elf_kindens(i_point,1)*elf_dens(i_point,1) - &
                    0.125*dot_product(elf_dens_grad(:,i_point,1),elf_dens_grad(:,i_point,1))
               minicube_output(i_point) = tp1**2/(tp1**2+tp2**2+1d-6)
            enddo! i_point
         elseif (cube_elf(i_cube).eq.1) then
            ! Becke-Edgecombe ELF, spin channel 1
            do i_point = 1, cube_divisor(i_cube)**3
               tp1 = 0.6*(2*pisq3)**(2*third)*(elf_dens(i_point,1))**(8*third)
               tp2 = elf_kindens(i_point,1)*elf_dens(i_point,1) - &
                    0.25*dot_product(elf_dens_grad(:,i_point,1),elf_dens_grad(:,i_point,1))
               minicube_output(i_point) = tp1**2/(tp1**2+tp2**2+1d-6)
            enddo! i_point
         elseif (cube_elf(i_cube).eq.2) then
            ! Becke-Edgecombe ELF, spin channel 2
            do i_point = 1, cube_divisor(i_cube)**3
               tp1 = 0.6*(2*pisq3)**(2*third)*(elf_dens(i_point,2))**(8*third)
               tp2 = elf_kindens(i_point,2)*elf_dens(i_point,2) - &
                    0.25*dot_product(elf_dens_grad(:,i_point,2),elf_dens_grad(:,i_point,2))
               minicube_output(i_point) = tp1**2/(tp1**2+tp2**2+1d-6)
            enddo! i_point
         else
            ! Kohout-Savin ELF
            elf_kindens(:,1) = elf_kindens(:,1) + elf_kindens(:,2)
            do i_point = 1, cube_divisor(i_cube)**3
               tp1 = 0.3*(2*pisq3)**(2*third)*&
                    (elf_dens(i_point,1)**(5*third) + elf_dens(i_point,2))*&
                    elf_dens(i_point,1)*elf_dens(i_point,2)
               tp2 = elf_kindens(i_point,1)*elf_dens(i_point,1)*elf_dens(i_point,2) - &
           0.125*dot_product(elf_dens_grad(:,i_point,1),elf_dens_grad(:,i_point,1))*elf_dens(i_point,2) - &
           0.125*dot_product(elf_dens_grad(:,i_point,2),elf_dens_grad(:,i_point,2))*elf_dens(i_point,1)
               minicube_output(i_point) = tp1**2/(tp1**2+tp2**2+1d-6)
            enddo
         endif ! if (cube_state(1,i_cube).eq.0)
      endif ! if (n_spin.eq.1)


   ! WPH: For output using SOC-perturbed eigenvectors, we calculate the minicubes for the spin-up spinor 
   !      half and the spin-down spinor half individually, then sum them together.  This allows us to use
   !      the scalar-relativistic cube subroutines for each half.
   case('eigenstate')
         if (out_cube_soc .and. cube_type_needs_soc_eigenvectors(i_cube)) then
                   allocate(KS_vec_complex_re(n_basis),stat=info)
                   call check_allocation(info, 'KS_vec_complex_re', func)
                   do i_basis = 1, n_basis
                      KS_vec_complex_re(i_basis) = real(KS_vec_complex(i_basis))
!                      write(info_str, '(I3, F10.6)') i_basis, KS_vec_complex_re(i_basis)
!                      call localorb_info(info_str)
                   end do
                   call evaluate_wavefunc_real(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube), &
                                              .false., n_basis, 1, KS_vec_complex_re, minicube_wavefunc)
                  minicube_output=minicube_wavefunc
                  ! AJL, Oct2018: Moved this to final deallocations during debugging
                  ! deallocate(KS_vec_complex_re)
!                  call evaluate_wavefunc_cmplx(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube),  &
!                                            .false., n_basis, 1, KS_vec_complex, minicube_wavefunc_complex)
!                  minicube_output=real(minicube_wavefunc_complex)
              ! Spin-up spinors
!              call evaluate_wavefunc_cmplx(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube),      &
!                                        .false., n_basis, 1, KS_vec_soc_up, minicube_wavefunc_complex)
!              minicube_output = real(minicube_wavefunc_complex)
!              ! Spin-down spinors
!              minicube_wavefunc_complex = (0.0d0,0.0d0)
!              call evaluate_wavefunc_cmplx(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube),      &
!                                        .false., n_basis, 1, KS_vec_soc_dn, minicube_wavefunc_complex)
!              minicube_output = minicube_output + real(minicube_wavefunc_complex)
         else
              if (real_eigenvectors) then
                  call evaluate_wavefunc_real(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube),   &
                                             .false., n_basis, 1, KS_vec, minicube_wavefunc)
                  minicube_output=minicube_wavefunc
              else
                  call evaluate_wavefunc_cmplx(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube),  &
                                            .false., n_basis, 1, KS_vec_complex, minicube_wavefunc_complex)
                  minicube_output=real(minicube_wavefunc_complex)
              endif
         end if
         minicube_output=minicube_output*sqrt_inv_bohr_3
    case('eigenstate_non_soc')
              if (real_eigenvectors) then
                  call evaluate_wavefunc_real(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube),   &
                                             .false., n_basis, 1, KS_vec, minicube_wavefunc)
                  minicube_output=minicube_wavefunc
              else
                  call evaluate_wavefunc_cmplx(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube),  &
                                            .false., n_basis, 1, KS_vec_complex, minicube_wavefunc_complex)
                  minicube_output=real(minicube_wavefunc_complex)
              endif
    case('eigenstate_imag')
         if (out_cube_soc .and. cube_type_needs_soc_eigenvectors(i_cube)) then
                   allocate(KS_vec_complex_im(n_basis),stat=info)
                   call check_allocation(info, 'KS_vec_complex_im', func)
                   do i_basis = 1, n_basis
                      KS_vec_complex_im(i_basis) = aimag(KS_vec_complex(i_basis))
                   end do
!                   KS_vec_complex = KS_vec_complex_im
!                   deallocate(KS_vec_complex_im) 
                   call evaluate_wavefunc_real(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube), &
                                              .false., n_basis, 1, KS_vec_complex_im, minicube_wavefunc)
                  minicube_output=minicube_wavefunc
                  ! AJL, Oct2018: Moved this to final deallocations during debugging
                  ! deallocate(KS_vec_complex_im)
!                   minicube_output=real(minicube_wavefunc_complex)
              ! Spin-up spinors
!              call evaluate_wavefunc_cmplx(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube),      &
!                                        .false., n_basis, 1, KS_vec_soc_up, minicube_wavefunc_complex)
!              minicube_output = aimag(minicube_wavefunc_complex)
!              ! Spin-down spinors
!              minicube_wavefunc_complex = (0.0d0,0.0d0)
!              call evaluate_wavefunc_cmplx(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube),      &
!                                        .false., n_basis, 1, KS_vec_soc_dn, minicube_wavefunc_complex)
!              minicube_output = minicube_output + aimag(minicube_wavefunc_complex)
        else
              if (real_eigenvectors) then
                    minicube_output=0.0d0
              else
                   call evaluate_wavefunc_cmplx(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube), &
                                              .false., n_basis, 1, KS_vec_complex, minicube_wavefunc_complex)
                   minicube_output=aimag(minicube_wavefunc_complex)
              endif
        end if
        minicube_output=minicube_output*sqrt_inv_bohr_3
     case('eigenstate_density')
        if (out_cube_soc) then
             ! Spin-up spinors
!             call evaluate_wavefunc_cmplx(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube),    &
!                                         .false., n_basis, 1, KS_vec_soc_up, minicube_wavefunc_complex)
!             minicube_wavefunc_complex_temp = minicube_wavefunc_complex
!              ! Spin-down spinors
!             minicube_wavefunc_complex = (0.0d0,0.0d0)
!             call evaluate_wavefunc_cmplx(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube),    &
!                                         .false., n_basis, 1, KS_vec_soc_dn, minicube_wavefunc_complex)
!             minicube_wavefunc_complex = minicube_wavefunc_complex + minicube_wavefunc_complex_temp
!             do i_minicube_point=1,cube_divisor(i_cube)**3,1
!                  minicube_output(i_minicube_point) = real(minicube_wavefunc_complex(i_minicube_point))**2 + &
!                                              aimag(minicube_wavefunc_complex(i_minicube_point))**2
!             end do
        else
             if (real_eigenvectors) then
                 call evaluate_wavefunc_real(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube),     &
                                            .false., n_basis, 1, KS_vec, minicube_wavefunc)
                 minicube_output=minicube_wavefunc**2
             else
                 call evaluate_wavefunc_cmplx(cube_divisor(i_cube)**3, minicube_points, cube_state(2,i_cube),    &
                                             .false., n_basis, 1, KS_vec_complex, minicube_wavefunc_complex)
                 do i_minicube_point=1,cube_divisor(i_cube)**3,1
                      minicube_output(i_minicube_point)= real(minicube_wavefunc_complex(i_minicube_point))**2+  &
                                                  aimag(minicube_wavefunc_complex(i_minicube_point))**2
                 enddo
             end if
        end if
        minicube_output=minicube_output*inv_bohr_3
   case('long_range_potential')
        do i_minicube_point=1,cube_divisor(i_cube)**3,1
                call  evaluate_potential(minicube_points(:,i_minicube_point), &
                                         dip_gradient,dip_origin,dip_length,  &
                                         pot_jump, minicube_output(i_minicube_point))
                call embedding_potential(minicube_points(:,i_minicube_point), &
                                         minicube_embedding_potential(i_minicube_point),dummy)
                minicube_output(i_minicube_point)=minicube_output(i_minicube_point)+minicube_embedding_potential(i_minicube_point)
         enddo
         minicube_output=(minicube_output-average_potential)*hartree
         if (myid.eq.0) then
             ! write(use_unit,*) 'Debug info: potential shift [H] ', average_potential
             ! write(use_unit,*) 'Debug info: potential shift [eV] ', average_potential*hartree
         endif


   case('potential')
        call aims_stop('Current implementation of output cube potential not compatible with distribution in minicubes')
   case('hartree_potential')
        call aims_stop('Current implementation of output cube potential not compatible with distribution in minicubes')
   case('xc_potential')
        call aims_stop('Current implementation of output cube potential not compatible with distribution in minicubes')
                 
   case default 
        call aims_stop('Unknown cube type during data mining') 
   end select
   
   !Now re-collect data in a single large array
                  do i_minicube_point = 1, cube_divisor(i_cube)**3, 1
                    dumpme=.false.
                    i_point = floor((i_minicube_point-1)/dble(cube_divisor(i_cube)**2))& 
                              *cube_edge_steps(1,i_cube)*cube_edge_steps(2,i_cube)+ &
                              floor((modulo(i_minicube_point-1,cube_divisor(i_cube)**2))/& 
                              dble(cube_divisor(i_cube)))*cube_edge_steps(1,i_cube) + &
                              modulo(((modulo(i_minicube_point-1,cube_divisor(i_cube)**2)+1)-1),cube_divisor(i_cube))+1  + &
                               (i_minicube_x-1)*cube_divisor(i_cube) + & !offset in x
                               (i_minicube_y-1)*cube_divisor(i_cube)*cube_edge_steps(1,i_cube) + & !offset in y
                               (i_minicube_z-1)*cube_divisor(i_cube)*cube_edge_steps(1,i_cube)*cube_edge_steps(2,i_cube) !offset in z
                 if (i_minicube_x.eq.n_minicubes_x .or. i_minicube_y.eq.n_minicubes_y .or. i_minicube_z.eq.n_minicubes_z) then !Borderline cube, dump unneccesary points
                      if (i_minicube_z.eq.n_minicubes_z) then
                         max_divisor = modulo(cube_edge_steps(3,i_cube)-1,cube_divisor(i_cube))+1
                         if (i_minicube_point.gt.max_divisor*cube_divisor(i_cube)**2) dumpme=.true. !dump point
                      endif
                      if (i_minicube_y.eq.n_minicubes_y) then
                        max_divisor=modulo(cube_edge_steps(2,i_cube)-1,cube_divisor(i_cube))+1
                        if (modulo(i_minicube_point-1,cube_divisor(i_cube)**2)+1.gt.max_divisor*cube_divisor(i_cube)) dumpme=.true.
                      endif
                      if (i_minicube_x.eq.n_minicubes_x) then
                          max_divisor=modulo(cube_edge_steps(1,i_cube)-1,cube_divisor(i_cube))+1
                          if (modulo(modulo(i_minicube_point-1,cube_divisor(i_cube)**2),cube_divisor(i_cube))+1.gt.max_divisor) & 
                             dumpme=.true. !dump
                      endif
                       
                 endif
                     if(.not.dumpme.and. .not.dumpcube) then
                       cube_output(i_point) = minicube_output(i_minicube_point)
                       ! Safeguard against too small values, they don't necessarily write
                       ! correctly to the cubefile due to lacking space for the exponent.
                       ! It seems that only two digits are safe here.
                       ! BL: 1e-99 is zero when converted to real(kind=4)
                       ! Is this "accuracy" necessary? 1e-16 should be sufficient.
                       if (ABS(cube_output(i_point)) < 1.0e-16) cube_output(i_point) = 0.0
                     endif
                  enddo
                                 
   enddo !i_task - over tasks
   call sync_vector(cube_output,n_points)
! elseif(cube_type(i_cube).eq.'ion_dens') then
!        if (myid==0) write(use_unit,*) 'catching it before'
!        n_minicubes=(ceiling(n_points/dble(cube_divisor(i_cube)**3)))
!        do i_minicube=1,n_minicubes,1
!           minicube_points=0.0 !Reset 
!           minicube_npoints=cube_divisor(i_cube)**3
!           if (i_minicube.eq.n_minicubes) minicube_npoints=mod(n_points-1,cube_divisor(i_cube)**3)+1 !Number of points for the minicube
!           do i_point = 1,minicube_npoints,1
!              call obtain_coordinate_of_point((i_minicube-1)*cube_divisor(i_cube)**3+i_point,i_cube,minicube_points(:,i_point))
!              !minicube_points(:,i_point)=points(:,(i_minicube-1)*cube_divisor(i_cube)**3+i_point) !Divide the points
!           enddo
! 	  call evaluate_potential_cube( minicube_npoints, minicube_points,  minicube_output  )
!         !  call output_whole_potential_on_arbitrary_points(minicube_npoints,minicube_points,minicube_output)
!           do i_point=1,minicube_npoints,1
!             cube_output((i_minicube-1)*cube_divisor(i_cube)**3+i_point)=minicube_output(i_point) !Collect in output array
!           enddo
!        enddo

elseif(cube_type(i_cube).eq.'potential') then
!  !    TODO: Make somehow compatible with actual division in minicubes and its distribution amongst the CPUs
!     if (myid.eq.0) then
!       write(info_str,*) ' Calculating potential for output on cube grid. This can take a while...'
!       call localorb_info(info_str)
!       write(info_str,*) '* This is an EXPERIMENTAL feature. Make sure to doublecheck the results before using for anything.'
!       call localorb_info(info_str)
!      write(info_str,*) '* The correctness of the results are not guaranteed, so proceed at your own risk.'
!       call localorb_info(info_str)
!      if (cube_divisor(i_cube).eq.10) then
!         write(info_str,*) 'Cube divisor was left at its default value of 10.'
!         call  localorb_info(info_str)
!         write(info_str,*) 'This will usually result in very slow processing.'
!         call  localorb_info(info_str)
!         write(info_str,*) 'For better results, try increasing cube divisor to 45.'
!         call  localorb_info(info_str)
!      endif
!     endif
!  !    To circumvent problem with too large arrays, split them into 90k points maximum size
!  !    FIXME: This is a very inelegant AND slow solution which should be improved
!  !    cube_divisor is only taken here as limit because this is the size of the minicube output...
!  !    TODO: Make somehow compatible with actual division in minicubes and its distribtion amongst the CPUs
!       n_minicubes=(ceiling(n_points/dble(cube_divisor(i_cube)**3)))
!       do i_minicube=1,n_minicubes,1
!          minicube_points=0.0 !Reset 
!          minicube_npoints=cube_divisor(i_cube)**3
!          if (i_minicube.eq.n_minicubes) minicube_npoints=mod(n_points-1,cube_divisor(i_cube)**3)+1 !Number of points for the minicube
!          do i_point = 1,minicube_npoints,1
!             call obtain_coordinate_of_point((i_minicube-1)*cube_divisor(i_cube)**3+i_point,i_cube,minicube_points(:,i_point))
!             !minicube_points(:,i_point)=points(:,(i_minicube-1)*cube_divisor(i_cube)**3+i_point) !Divide the points
!          enddo
!          call output_whole_potential_on_arbitrary_points(minicube_npoints,minicube_points,minicube_output)
!          do i_point=1,minicube_npoints,1
!            cube_output((i_minicube-1)*cube_divisor(i_cube)**3+i_point)=minicube_output(i_point) !Collect in output array
!          enddo
!       enddo
!       cube_output(:)=(cube_output(:)-average_potential)*hartree !Provide unit
!       !Shift by average potential
!       write(info_str,*) ' Attention: The potential is NOT shifted to the average potential being zero and '
!       call localorb_info(info_str)
!       write(info_str,*) ' therefore does not correspond directly to the numbers in the output '
!       call localorb_info(info_str)
!
!     !write(info_str,*) ' Calculating potential for output on cube grid. This can take a while...'
!     !call localorb_info(info_str)
!     !call output_whole_potential_on_arbitrary_points(n_points,points,cube_output)
!     !cube_output(:)=cube_output(:)*hartree

elseif((cube_type(i_cube).eq.'hartree_potential').or.(cube_type(i_cube).eq.'xc_potential')) then
    if(.not. allocated( radius_esp_min))then
      allocate( radius_esp_min(n_species),stat=info)
      if(info/=0)then
        write(use_unit,*)'Error in allocation: radius_esp_min'
	stop
      end if
    end if
    if(.not. allocated( radius_esp_max))then
      allocate( radius_esp_max(n_species),stat=info)
      if(info/=0)then
	write(use_unit,*)'Error in allocation: radius_esp_max'
	stop
      end if
    end if
    radius_esp_min = 0d0
    radius_esp_max = 10d0
    call esp_partition_grid( radius_esp_min, radius_esp_max, .False.,&
                               1,.true.,.false.,&
    cube_edge_steps(1,i_cube), cube_edge_steps(2,i_cube), cube_edge_steps(3,i_cube),&
           .true.,cube_edge_unit(1:3,1:3,i_cube), offset_coord - cube_origin(1:3,i_cube))

    if(.not. allocated( partition_tab_esp))then
      allocate( partition_tab_esp(n_full_points_esp),stat=info)
      if(info/=0)then
	write(use_unit,*)'Error in allocation: partition_tab_esp'
	stop
      end if
    end if

    if(n_periodic.eq.0)then
       partition_tab_esp = 1d0
    else
      partition_tab_esp(:) = (1d0/dble(cube_edge_steps(1,i_cube)* cube_edge_steps(2,i_cube)* &
                                cube_edge_steps(3,i_cube)))*cell_volume
    endif


    if (cube_type(i_cube).eq.'hartree_potential') then
    
      if (allocated(potential_trans)) then
	deallocate(potential_trans)
      end if   
      allocate( potential_trans(n_full_points_esp,1),stat=info)
      if(info/=0)then
	write(use_unit,*)'Error in allocation: potential_trans'
	stop
      end if
      
      if(.not. allocated( delta_v_hartree_part_at_zero_trans))then
	allocate( delta_v_hartree_part_at_zero_trans(n_atoms),stat=info)
	if(info/=0)then
	    write(use_unit,*)'Error in allocation: delta_v_hartree_part_at_zero_trans'
	    stop
	end if
      end if
      if(.not. allocated( delta_v_hartree_deriv_l0_at_zero_trans))then
	allocate( delta_v_hartree_deriv_l0_at_zero_trans(3,n_atoms),stat=info)
	if(info/=0)then
	    write(use_unit,*)'Error in allocation: delta_v_hartree_deriv_l0_at_zero_trans'
	    stop
	end if
      end if
      if(.not. allocated( multipole_moments_trans))then
	allocate( multipole_moments_trans(( l_pot_max + 1)**2, n_atoms),stat=info)
	if(info/=0)then
	    write(use_unit,*)'Error in allocation: multipole_moments_trans'
	    stop
	end if
      end if
      if(.not. allocated( multipole_radius_sq_trans))then
	allocate( multipole_radius_sq_trans(n_atoms),stat=info)
	if(info/=0)then
	    write(use_unit,*)'Error in allocation: multipole_radius_sq_trans'
	    stop
	end if
      end if
      if(.not. allocated( l_hartree_max_far_distance_trans))then
	allocate( l_hartree_max_far_distance_trans(n_atoms),stat=info)
	if(info/=0)then
	    write(use_unit,*)'Error in allocation: l_hartree_max_far_distance_trans'
	    stop
	end if
      end if
      if(.not. allocated( outer_potential_radius_trans))then
	allocate( outer_potential_radius_trans(0:l_pot_max, n_atoms),stat=info)
	if(info/=0)then
	    write(use_unit,*)'Error in allocation: outer_potential_radius_trans'
	    stop
	end if
      end if

      call multipole_expantion_rho &
	( hartree_partition_tab, rho, &
	delta_v_hartree_part_at_zero_trans, &
	delta_v_hartree_deriv_l0_at_zero_trans, &
	multipole_moments_trans, multipole_radius_sq_trans, &
	l_hartree_max_far_distance_trans, &
	outer_potential_radius_trans, free_rho_superpos, .true., .false., 1)


      if(.not. allocated( free_hartree_superpos_trans))then
	allocate( free_hartree_superpos_trans(n_full_points_esp),stat=info)
	if(info/=0)then
	   write(use_unit,*)'Error in allocation: free_hartree_superpos_trans'
	   stop
	end if
      end if


      call esp_grid_storage ( partition_tab_esp, &
                                free_hartree_superpos_trans,&
                                 .true. )
  !   if (flag_out_locpot_atom) then 
  !          call sum_up_potential &
!	( delta_v_hartree_part_at_zero_trans, &
!	  multipole_moments_trans, &
!	  partition_tab_esp, potential_trans, &
!	  multipole_radius_sq_trans, &
!	  l_hartree_max_far_distance_trans,  &
!	  outer_potential_radius_trans, &
!	  free_hartree_superpos_trans, &
!	  .false., .false.,use_dipole_correction)
 !   else 
     call sum_up_potential &
	( delta_v_hartree_part_at_zero_trans, &
	  multipole_moments_trans, &
	  partition_tab_esp, potential_trans, &
	  multipole_radius_sq_trans, &
	  l_hartree_max_far_distance_trans,  &
	  outer_potential_radius_trans, &
	  free_hartree_superpos_trans, &
	  .true., .false.,use_dipole_correction)
 !   endif 

    else
      if (allocated(potential_trans)) then
	deallocate(potential_trans)
      end if
      allocate( potential_trans( n_full_points_esp, n_spin),stat=info)
      if(info/=0)then
	write(use_unit,*)'Error in allocation: potential_trans'
	stop
      end if
      if (spin_treatment.eq.1) then
      	if (.not.allocated(hamiltonian_two)) then
	  allocate( hamiltonian_two(n_hamiltonian_matrix_size, n_spin),&
			stat=info)
	  call check_allocation(info, 'hamiltonian_two                   ')
	end if
      else
	if (.not.allocated(hamiltonian_two)) then
	  allocate(hamiltonian_two(1,1))
	end if
      endif
      if (myid.eq.0) then
	write (use_unit,'(2X,A)') &
	  "Calculating XC-Potential."
      end if
      call esp_get_n_compute_dens( partition_tab_esp )
      call get_vxc_esp &      
        ( KS_eigenvector, KS_eigenvector_complex, occ_numbers, &
	partition_tab_esp,  l_shell_max, &
	potential_trans,  hamiltonian(1,1),hamiltonian_two(1,1),1,&
	1,1,.true., 2)
    endif

    call esp_collect_pot(potential_trans(:,cube_state(1,i_cube)),partition_tab_esp,cube_output,&
                         cube_edge_steps(1,i_cube), cube_edge_steps(2,i_cube),&
                         cube_edge_steps(3,i_cube))
    if (allocated(delta_v_hartree_part_at_zero_trans)) then
	deallocate(delta_v_hartree_part_at_zero_trans)
    end if
    if (allocated(delta_v_hartree_deriv_l0_at_zero_trans)) then
	deallocate(delta_v_hartree_deriv_l0_at_zero_trans)
    end if
    if (allocated(multipole_moments_trans)) then
	deallocate(multipole_moments_trans)
    end if
    if (allocated(multipole_radius_sq_trans)) then
	deallocate(multipole_radius_sq_trans)
    end if
    if (allocated(l_hartree_max_far_distance_trans)) then
	deallocate(l_hartree_max_far_distance_trans)
    end if
    if (allocated(outer_potential_radius_trans)) then
	deallocate(outer_potential_radius_trans)
    end if
    if (allocated(radius_esp_min)) then
	deallocate(radius_esp_min)
    end if
    if (allocated(radius_esp_max)) then
	deallocate(radius_esp_max)
    end if
    if (allocated(partition_tab_esp)) then
	deallocate(partition_tab_esp)
    end if
    if (allocated(potential_trans)) then
	deallocate(potential_trans)
    end if
    if (allocated(free_hartree_superpos_trans)) then
	deallocate(free_hartree_superpos_trans)
    end if
    if (allocated(hamiltonian_two)) then
	deallocate(hamiltonian_two)
    end if
endif
   

     ! Output of the results - on central processor only 
     if (myid.eq.0) then

           !Open file
           open(unit= 10+i_cube, file=cube_filename(i_cube), ACTION='WRITE')
   
           !Write the header
          if (cube_format(i_cube).eq.'gOpenMol') flag_gOpenMol_format=.true.  !Kept for compatibility, just in case. 
          select case (cube_format(i_cube))
            case ('xsf')
              call write_cube_header_xsf(10+i_cube,i_cube)
              write(file_format,'(A)') '(1E13.5, $)'
            case('gOpenMol')
              call write_cube_header_gOpenMol(10+i_cube, cube_units,i_cube)
              write(file_format,'(A)') '(1E13.5)'
            case default
               call write_cube_header(10+i_cube,offset_coord - cube_origin(1:3,i_cube),cube_edge_steps(1:3,i_cube), i_cube)
              write(file_format,'(A)') '(1E16.8, $)'
           end select

         blocksize = 1  !This is incredibly inelegant, but I can't think of a better way, and omitting the 
                                  !carriage return at either point causes some visualizers to display nonsense results
                                  !if the output gets larger then 9000 points o.O

                  do i_point = 1,n_points,1
                      write (10+i_cube,fmt=file_format) &
                           cube_output(i_point)
                       if(mod(blocksize,6).eq.0) then
                          write(10+i_cube,*) ""
                          blocksize=0
                      elseif(mod(i_point,cube_edge_steps(3,i_cube)).eq.0) then
                           write(10+i_cube,*) ""
                           blocksize=0
                      endif
                      blocksize= blocksize+1
                  enddo
     !Add a carriage return at the end
      if(cube_format(i_cube).ne.'xsf')  write(10+i_cube,*) ""

        if (cube_format(i_cube).eq.'xsf') then 
            write(10+i_cube,*) 'END_DATAGRID_3D'
            write(10+i_cube,*) 'END_BLOCK_DATAGRID_3D'
        endif
        !close the file
        close(10+i_cube)
  endif !myid = 0
  !Secondary output of files, e.g. for stm data
  if (myid.eq.0) then
    if (cube_type(i_cube)=='stm') then 
       write(stm_filename,'(A,I5.5,A)') 'stm_z_map_', i_cube, '.cube'
       open(unit=9,file=stm_filename,ACTION='WRITE')
       select case (cube_format(i_cube))
            case ('xsf')
              call write_cube_header_xsf(9,i_cube)
              write(file_format,'(A)') '(1E13.5, $)'
            case('gOpenMol')
              call write_cube_header_gOpenMol(9, cube_units,i_cube)
              write(file_format,'(A)') '(1E13.5)'
            case default
               call write_cube_header(9,offset_coord - cube_origin(1:3,i_cube),cube_edge_steps(1:3,i_cube), i_cube)
              write(file_format,'(A)') '(1E13.5, $)'
       end select
       call write_wsxm_data(i_cube, cube_output, n_points)
       call output_cube_z_map(i_cube,n_points,file_format)
       if (cube_format(i_cube).eq.'xsf') then
            write(10+i_cube,*) 'END_DATAGRID_3D'
            write(10+i_cube,*) 'END_BLOCK_DATAGRID_3D'
      endif
      close(9)
    endif
  endif

  !Secondary output for average
  if(cube_average(i_cube).eqv..true.) then
    allocate(average(cube_edge_steps(3,i_cube)))
    average=0.0
    points_in_plane=cube_edge_steps(1,i_cube)*cube_edge_steps(2,i_cube)
    do i_plane=1,cube_edge_steps(3,i_cube),1
       do i_y=1,cube_edge_steps(2,i_cube)
        do i_x=1,cube_edge_steps(1,i_cube)
            i_point=(i_x-1)*(cube_edge_steps(3,i_cube)*cube_edge_steps(2,i_cube)) &
                           +(i_y-1)*cube_edge_steps(3,i_cube)+i_plane
            average(i_plane)=average(i_plane)+cube_output(i_point)
        enddo
       enddo
    enddo
    average=average/points_in_plane

    write(average_filename,'(A,A)') 'plane_average_', trim(cube_filename(i_cube))! 'potential'
    if (myid.eq.0) then
    !write(use_unit,*) ' Debug info: ', average_filename
    open(unit=9,file=average_filename) !TODO: Make Filename dynamic
    if (flag_out_elec_real) then 
       write (9,fmt='(A,F15.8)') "#Numerical average free-atom electrostatic potential in [eV]:", average_free_es_pot*hartree
       write (9,fmt='(A,F15.8)') "#Average real-space part of the electrostatic  potential in [eV]:", average_potential*hartree
       write(9,fmt='(2X,A,2X,A)') '#z-coordinate (Ang)', 'plane-average (hartree)'
    else
      write(9,fmt='(2X,A,2X,A)') 'z-coordinate (Ang)', 'plane-average (unit as in parent cube file)'
    endif
    do i_plane=1,cube_edge_steps(3,i_cube)
         current_coord=cube_origin(3,i_cube)+(i_plane-((cube_edge_steps(3,i_cube)+1)/2.0))*cube_units(3)
         current_coord=current_coord*bohr
         write(9,*) current_coord, average(i_plane) !TODO: Add z-coordinate
    enddo
    close(9)
    endif
    


    deallocate(average)
  endif



     !Finally, deallocate
     if (allocated(densmat))                      deallocate(densmat)
     if (allocated(densmat_sparse))               deallocate(densmat_sparse)
     if (allocated(KS_vec))                       deallocate(KS_vec)
     if (allocated(KS_vec_complex))               deallocate(KS_vec_complex)
     if (allocated(KS_vec_complex_im))            deallocate(KS_vec_complex_im)
     if (allocated(KS_vec_complex_re))            deallocate(KS_vec_complex_re)
     if (allocated(KS_vec_soc_up))                deallocate(KS_vec_soc_up)
     if (allocated(KS_vec_soc_dn))                deallocate(KS_vec_soc_dn)

     !if (allocated(points))                      deallocate(points)
     if (allocated(minicube_points))              call aims_deallocate(minicube_points, "minicube_points")
     if (allocated(cube_output))                  deallocate(cube_output)
     if (allocated(minicube_output))              call aims_deallocate(minicube_output, "minicube_output")
     if (allocated(minicube_kindens))             deallocate(minicube_kindens)
     if (allocated(minicube_wavefunc_complex))    call aims_deallocate(minicube_wavefunc_complex, "minicube_wavefunc_complex")
     if (allocated(minicube_wavefunc_complex_temp))    &
           call aims_deallocate(minicube_wavefunc_complex_temp, "minicube_wavefunc_complex_temp")
     if (allocated(minicube_density))             deallocate (minicube_density)
     if (allocated(minicube_density_gradient))    deallocate (minicube_density_gradient)
     if (allocated(minicube_initial_density))     deallocate(minicube_initial_density)
     if (allocated(minicube_wavefunc))            call aims_deallocate(minicube_wavefunc, "minicube_wavefunc")
     if (allocated(minicube_D))                   deallocate (minicube_D)
     if (allocated(minicube_D0))                  deallocate (minicube_D0)
     if (allocated(minicube_embedding_potential)) deallocate (minicube_embedding_potential)
     if (allocated(elf_dens))                     deallocate (elf_dens)
     if (allocated(elf_dens_grad))                deallocate (elf_dens_grad)
     if (allocated(elf_kindens))                  deallocate (elf_kindens)

     !SR: MPB solvation effects
     if (allocated(dielec_func)) deallocate(dielec_func)
     if (allocated(ion_dens)) deallocate(ion_dens)
     if (allocated(kappa_func)) deallocate(kappa_func)
     if (allocated(delta_v_cube)) deallocate(delta_v_cube)
     if (allocated(delta_v_gradient_abs)) deallocate(delta_v_gradient_abs)
     if (use_separate_alpha) then
        if (allocated(kappa_func_anion)) deallocate(kappa_func_anion)
     endif
     if (allocated(my_free_hartree_superpos)) deallocate(my_free_hartree_superpos)
     if (allocated(my_free_rho_superpos)) deallocate(my_free_rho_superpos)
     if (allocated(my_rho)) deallocate(my_rho)
     if (allocated(my_rho_gradient)) deallocate(my_rho_gradient)
     if (allocated(my_free_rho_gradient_superpos)) deallocate(my_free_rho_gradient_superpos)
     if (allocated(dielec_func_grad)) deallocate(dielec_func_grad)
     if (allocated(kappa_func_gradient)) deallocate(kappa_func_gradient)   
     if (use_separate_alpha) then
        if (allocated(kappa_func_gradient_anion)) deallocate(kappa_func_gradient_anion)
     end if
     if (allocated(minicube_density_gradient_tp)) deallocate(minicube_density_gradient_tp) 
    !Done with this cube - tell the user
    if (myid.eq.0) then
     write(info_str,fmt='(A,I0,A)') '   Cube file ', i_cube, ' finished.'
     call localorb_info(info_str)
     call localorb_info('')
    endif


   enddo !i_cube to n_cubes

end subroutine output_cube_files_p2_cpscf


