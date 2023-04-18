!****s* FHI-aims/output_cube_files_p1
!  NAME
!  output_cube_files_p1
!  SYNOPSIS
subroutine output_cube_files_p1 ()
  !  PURPOSE
  !  Plot charge density for periodic systems. This routine can still
  !  be very improved, see FIXME's below.
  !
  !  USES
  !
  use timing
  use mpi_tasks
  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use mpi_utilities
  use localorb_io
  use basis
  use cartesian_ylm
  use constants
  use plot
  use physics
  use synchronize_mpi
  use density_matrix_evaluation
  use gt
  use pbc_lists

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

  !  local variables
  real*8 inv_bohr_3
  real*8 sqrt_inv_bohr_3
  character*160 :: info_str

  real*8 coord_current(3), coord_current_temp(3)
  real*8 dist_tab(n_centers_basis_integrals)
  real*8 dist_tab_sq(n_centers_basis_integrals, cube_edge_steps(3,1))
  real*8 i_r(n_centers_basis_integrals)
  real*8 dir_tab(3, n_centers_basis_integrals, cube_edge_steps(3,1))
  real*8 dir_tab_norm(3, n_centers_basis_integrals)
  real*8 trigonom_tab(4, n_centers_basis_integrals)
  real*8 en_upper_limit,en_lower_limit

  real*8,dimension(:),allocatable:: radial_wave
  real*8,dimension(:,:),allocatable:: wave


  integer :: n_compute, n_compute_a
  integer,dimension(:),allocatable :: i_basis
  integer :: i_cube

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)



  integer :: n_compute_atoms, i_k_point, i_k

  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)


  integer :: spline_array_start(n_centers_integrals)
  integer :: spline_array_end(n_centers_integrals)

  integer :: i_debug

  !     other local variables
  integer, dimension(n_spin,n_k_points) :: max_occ_number
  real*8, dimension(n_states,n_spin) :: occ_numbers_sqrt

  integer :: l_ylm_max
  integer :: n_points

  logical, dimension(:),allocatable :: if_cube
  logical :: if_stm
  real*8,     dimension(:,:,:),allocatable :: KS_vec
  integer :: cube_edge_steps_max
  complex*16, dimension(:,:,:),allocatable :: KS_vec_complex
  integer :: eig_dens_count, n_compute_basis_local, n_compute_basis
  integer :: n_points_per_task
  real*8,     dimension(:,:,:),allocatable :: densmat
  real*8,     dimension(:,:),allocatable :: densmat_tmp
  real*8,     dimension(:,:),allocatable :: densmat_sparse
  real*8,     dimension(:),allocatable :: densmat_sparse_tmp
  real*8,     dimension(:,:),allocatable :: density_matrix_con
  real*8,     dimension(:,:),allocatable :: work
  real*8,     dimension(:,:,:),allocatable :: occ_num

  real*8,     dimension(:,:,:),allocatable :: KS_ev_compute
  complex*16, dimension(:,:,:),allocatable :: KS_ev_compute_complex


  real*8, dimension(:,:,:), allocatable ::  KS_orbital
  complex*16, dimension(:,:,:), allocatable ::  KS_orbital_complex
  real*8, dimension(:,:,:), allocatable :: local_rho
  real*8, dimension(:,:), allocatable :: local_free_rho
  real*8, dimension(:),   allocatable :: local_rho_temp
  real*8, dimension(:,:,:), allocatable ::  local_orb
  logical, dimension(:), allocatable :: point_in_cell

  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:), allocatable :: ylm_tab

  real*8, dimension(:,:,:), allocatable :: dir_tab_global
  real*8, dimension(:,:), allocatable :: dist_tab_sq_global
  character(LEN=15):: file_format

  real*8 offset_coord(3)
  real*8 cube_units(3)

  !character*100 cube_filename(i_cube)
  

  !     counters
  integer :: coord_x, coord_y, coord_z, i_z
  integer :: i_l
  integer :: i_m
  integer :: i_state
  integer :: i_point
  integer :: off_z, last_z
  integer :: i_full_points
  integer :: i_coord

  integer :: i_spin = 1
  integer :: i_index, i_bas
  integer:: info, n_batch, n_batch_size, i_batch_z

  character(*), parameter :: deffmt = '2X'
  real*8 :: cpu, clock
  real*8 :: cpu_tot, clock_tot
  real*8 :: cpu_calc, clock_calc, cpu_sync, clock_sync, cpu_write, clock_write

  integer :: n_syncs
  logical :: cube_parallelized

  cpu_calc = 0.d0
  clock_calc = 0.d0
  cpu_sync = 0.d0
  clock_sync = 0.d0
  cpu_write = 0.d0
  clock_write = 0.d0


  ! TODO: figure out automatically which region to plot
  ! FIXME: ATM all plots use the same grid  (grid of the first plot)
  ! FIXME: ATM only canonical cartesian basis is supported 

  !     begin work


 !TODO: Change here to allow non-standard grids
     cube_units(1)= cube_edge_unit(1,1,1)
     cube_units(2)= cube_edge_unit(2,2,1)
     cube_units(3)= cube_edge_unit(3,3,1)
     offset_coord = real(cube_edge_steps(1:3,1)-1)/2.d0* cube_units 


!OTH: Added to retain compitability between old p1-Version and output_cube_p2
  do i_cube = 0,n_cube,1
      if (cube_format(i_cube).eq.'gOpenMol') flag_gOpenMol_format=.true. 
      if (cube_type(i_cube).eq.'stm')  cube_type(i_cube)='eigenstate_density'
  enddo

  call get_timestamps(cpu_tot, clock_tot)

  if(.not. allocated(if_cube))then
     allocate( if_cube(n_cube),stat=info )
     if(info/=0)then
        write(use_unit,*)'Error in allocation: if_cube(n_cube)'
        stop
     end if
  end if

  if_cube=.true.


  cube_edge_steps_max=0
  do i_z=1,3,1
     cube_edge_steps_max=MAX(cube_edge_steps(i_z,1),cube_edge_steps_max)
  enddo

  if (myid.eq.0) then
     write(use_unit,'(A)')  "------------------------------------------------------------"
     write(use_unit,'(2X,A)') "Writing density files:"
  endif
  do i_cube = 1, n_cube, 1
     if (cube_type(i_cube).eq.'spin_density') then
        cube_type(i_cube)='total_density'
        cube_spin(1,i_cube)=1
        cube_spin(2,i_cube)=-1
     endif

     if (cube_type(i_cube).eq.'total_density') then
        if(myid.eq.0) then
           if(cube_spin(1,i_cube).eq.1.and.cube_spin(2,i_cube).eq.1) then
              cube_filename(i_cube) ="total_density.cube"
           elseif (cube_spin(1,i_cube).eq.1.and.cube_spin(2,i_cube).eq.-1) then
              write (unit=cube_filename(i_cube),fmt='(A)') 'spin_density.cube'
           elseif (cube_spin(1,i_cube).eq.1.and.cube_spin(2,i_cube).eq.0) then
              write (unit=cube_filename(i_cube),fmt='(A)') 'density_spin_up.cube'
           elseif (cube_spin(1,i_cube).eq.0.and.cube_spin(2,i_cube).eq.1) then
              write (unit=cube_filename(i_cube),fmt='(A)') 'density_spin_down.cube'
           else
              write (unit=cube_filename(i_cube),fmt='(A,I2.2,A)') 'density_user_defined_spin_mask', i_cube, '.cube'
           endif
        endif
     else if (cube_type(i_cube).eq.'delta_density') then
        if(myid.eq.0) cube_filename(i_cube) ="delta_density.cube"
     else if (cube_type(i_cube).eq.'eigenstate'.or.cube_type(i_cube).eq.'eigenstate_density') then
        if(myid.eq.0) then
           if (cube_type(i_cube).eq.'eigenstate') then
              if(n_periodic .gt. 0) then
                 write (unit=cube_filename(i_cube),fmt='(A,I5.5,A,I1.1,A,I4.4,A)') 'eigenstate_',cube_index(i_cube),"_spin_", &
                 cube_state(1,i_cube), "_k_point_", cube_state(2,i_cube), ".cube" 
                 if(.not.real_eigenvectors)then
                    write(use_unit,*) "Eigenstate is requested, but the states are complex.", &
                    "The real part will be printed."
                 endif
              else
                 cube_state(2,i_cube)=1
                 write (unit=cube_filename(i_cube),fmt='(A,I5.5,A,I1.1,A)') 'eigenstate_',cube_index(i_cube),"_spin_",&
                 cube_state(1,i_cube), ".cube"
              endif
           else 
              if(cube_stm(3,i_cube).lt.0) then
                 if(n_periodic .gt. 0) then
                    write (unit=cube_filename(i_cube),fmt='(A,I5.5,A,I1.1,A,I4.4,A)') 'eigenstate_density_', &
                    cube_index(i_cube),"_spin_", &
                    cube_state(1,i_cube), "_k_point_", cube_state(2,i_cube),".cube"
                 else
                    cube_state(2,i_cube)=1
                    write (unit=cube_filename(i_cube),fmt='(A,I5.5,A,I1.1,A)') 'eigenstate_density_', &
                    cube_index(i_cube),"_spin_", cube_state(1,i_cube), ".cube"
                 endif
              endif
           endif
        endif

        if(cube_stm(3,i_cube).gt.0) then 
           cube_stm(1,i_cube)=chemical_potential
           if(cube_stm(2,i_cube).le.0)then
              en_upper_limit=cube_stm(1,i_cube)
              en_lower_limit=cube_stm(1,i_cube)+cube_stm(2,i_cube)
           endif
           if(cube_stm(2,i_cube).ge.0)then
              en_upper_limit=cube_stm(1,i_cube)+cube_stm(2,i_cube)
              en_lower_limit=cube_stm(1,i_cube)
           endif

           if(myid.eq.0) then
              eig_dens_count=0
              do i_k_point=1,n_k_points,1
                 do i_state=1,n_states,1
                    do i_spin=1,n_spin,1
                       if(KS_eigenvalue(i_state,i_spin,i_k_point).le.en_upper_limit.and. &
                       KS_eigenvalue(i_state,i_spin,i_k_point).ge.en_lower_limit) then
                          eig_dens_count=eig_dens_count+1
                       endif
                    enddo
                 enddo
              enddo

              if(eig_dens_count.eq.0) then
                 write(use_unit,*) "Cube files for eigenstate densities are requested, but no states fall", &
                 " within the specified energy window. The output contains all zeros"
                 if_cube(i_cube)=.false.
              endif
              write(unit=cube_filename(i_cube),fmt='(A,I2.2,A)') &
              'stm_',cube_index(i_cube),'.cube'
           endif
        endif
     endif

     if(myid.eq.0) then
        !        if(if_cube(i_cube)) then

        write(use_unit,'(2X,A,1X,A)') "| Creating cube file   :", cube_filename(i_cube)

        open(unit= 10+i_cube, file=cube_filename(i_cube), ACTION='WRITE')

        if(flag_gOpenMol_format)then
           call write_cube_header_gOpenMol(10+i_cube, cube_units,i_cube)
           write(file_format,'(A)') '(1E13.5)'
        else

           call write_cube_header(10+i_cube,offset_coord - cube_origin(1:3,1),cube_edge_steps(1:3,1), i_cube)
           write(file_format,'(A)') '(1E13.5, $)'
        end if
        !        endif
     endif

  enddo


  if(myid.eq.0) then
     if_stm=.false.
     do i_cube=1,n_cube,1
        if(cube_stm(3,i_cube).gt.0) then
           if_stm=.true.
           open(unit=9,file='stm_z_map.cube',ACTION='WRITE')
           if(flag_gOpenMol_format)then
              call write_cube_header_gOpenMol(9, cube_units,i_cube)
           else
              call write_cube_header(9,offset_coord - cube_origin(1:3,1),cube_edge_steps(1:3,1), i_cube)              
           endif
           exit
        endif
     enddo
  endif

  inv_bohr_3 = 1.0d0/(bohr**3)
  sqrt_inv_bohr_3 = sqrt(inv_bohr_3)

  l_ylm_max = l_wave_max

  ! First allocations

  if(.not. allocated( ylm_tab))then
     allocate( ylm_tab( (l_ylm_max+1)**2,n_centers_integrals),stat=info )
     if(info/=0)then
        write(use_unit,*)'Error in allocation: ylm_tab'
        stop
     end if
  end if


  if(.not. allocated( index_lm))then
     allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),stat=info) 
     if(info/=0)then
        write(use_unit,*)'Error in allocation: index_lm'
        stop
     end if
  end if

  !     initialize index_lm
  i_index = 0
  do i_l = 0, l_ylm_max, 1
     do i_m = -i_l, i_l
        i_index = i_index + 1
        index_lm(i_m, i_l) = i_index
     enddo
  enddo

  do i_k = 1, n_k_points
     !     find the maximal occupation number
     do i_spin = 1, n_spin, 1
        ! initialize
        max_occ_number(i_spin,i_k) = 0
        do i_state = n_states, 1, -1
           if (dabs(occ_numbers(i_state,i_spin,i_k)).gt.0.d0) then
              max_occ_number(i_spin,i_k) = i_state
              exit
           endif
        enddo
     enddo
  end do

  if(.not. allocated(i_basis))then
     allocate(i_basis(n_centers_basis_T),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: i_basis'
        stop
     end if
  end if

  if(.not. allocated(occ_num))then
     allocate(occ_num(n_states,n_spin,n_k_points),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: occ_num'
        stop
     end if
  end if


  if(.not. allocated(densmat))then
     allocate(densmat(n_centers_basis_T,n_centers_basis_T,n_cube),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: densmat'
        stop
     end if
  end if
  if(.not. allocated(densmat_tmp))then
     allocate(densmat_tmp(n_centers_basis_T,n_centers_basis_T),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: densmat_tmp'
        stop
     end if
  end if
  if(.not. allocated(densmat_sparse))then
     allocate(densmat_sparse(n_hamiltonian_matrix_size,n_cube),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: densmat_sparse'
        stop
     end if
  end if
  if(.not. allocated(densmat_sparse_tmp))then
     allocate(densmat_sparse_tmp(n_hamiltonian_matrix_size),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: densmat_sparse'
        stop
     end if
  end if


  allocate( local_rho(1:cube_edge_steps(3,1),1:2,n_cube) ) 
  allocate( local_free_rho(1:cube_edge_steps(3,1),1:2) ) 
  allocate( local_rho_temp(cube_edge_steps(3,1)))
  allocate( point_in_cell(1:cube_edge_steps(3,1))) 
  ! allocate(local_orb(maxval(max_occ_number(1,:)), 1:cube_edge_steps(3,1),1:2) )
  ! local_orb = 0.d0

  ! Calculate density matrices if needed

  call kweight_occs('output_cube_files_p1', occ_numbers)
  do i_cube=1,n_cube,1
     if(cube_type(i_cube).ne.'eigenstate') then

        occ_num=0.0
        densmat(1:n_centers_basis_T,1:n_centers_basis_T,i_cube)=0.0
        densmat_sparse(1:n_hamiltonian_matrix_size,i_cube)=0.0

        if(cube_type(i_cube).eq.'eigenstate_density') then

           if(cube_stm(3,i_cube).gt.0) then
              cube_stm(1,i_cube)=chemical_potential
              if(cube_stm(2,i_cube).le.0)then
                 en_upper_limit=cube_stm(1,i_cube)
                 en_lower_limit=cube_stm(1,i_cube)+cube_stm(2,i_cube)
              endif
              if(cube_stm(2,i_cube).ge.0)then
                 en_upper_limit=cube_stm(1,i_cube)+cube_stm(2,i_cube)
                 en_lower_limit=cube_stm(1,i_cube)
              endif
              do i_k_point=1,n_k_points,1
                 do i_state=1,n_states,1
                    do i_spin=1,n_spin,1
                       if(KS_eigenvalue(i_state,i_spin,i_k_point).le.en_upper_limit.and. &
                       KS_eigenvalue(i_state,i_spin,i_k_point).ge.en_lower_limit) then
                          occ_num(i_state,i_spin,i_k_point)=1.0
                       endif
                    enddo
                 enddo
              enddo

           endif ! if(cube_stm)
        elseif (cube_type(i_cube).eq.'total_density'.or.cube_type(i_cube).eq.'delta_density') then
           occ_num=occ_numbers

           en_lower_limit=-1.d10
           en_upper_limit=1.d10
        endif ! if(cube_type(i_cube).eq.'eigenstate_density')

        if(cube_type(i_cube).eq.'total_density'.or.cube_type(i_cube).eq.'delta_density'.or.& 
           (cube_type(i_cube).eq.'eigenstate_density'.and.cube_stm(3,i_cube).gt.0)) then
           do i_spin=1,n_spin,1
              call evaluate_densmat_part(KS_eigenvector, KS_eigenvector_complex, occ_num, &
              densmat_tmp, densmat_sparse_tmp, i_spin, &
              KS_eigenvalue, en_lower_limit, en_upper_limit)
              if(packed_matrix_format == PM_none)then
                 densmat(1:n_centers_basis_T,1:n_centers_basis_T,i_cube)=&
                 densmat(1:n_centers_basis_T,1:n_centers_basis_T,i_cube)+&
                 cube_spin(i_spin,i_cube)*densmat_tmp(1:n_centers_basis_T,1:n_centers_basis_T)
              else
                 densmat_sparse(1:n_hamiltonian_matrix_size,i_cube)=&
                 densmat_sparse(1:n_hamiltonian_matrix_size,i_cube)+&
                 cube_spin(i_spin,i_cube)*densmat_sparse_tmp(1:n_hamiltonian_matrix_size)
              endif
           enddo
        endif
     endif

!---------------------------------------------------
!here uses scGW density if required
  if(use_scgw .or. use_scgw0)then
    if(cube_type(i_cube).eq.'total_density') then
      if(n_spin .ne. 1)then
        if(myid.eq.0)then
          write(use_unit,*) "No density plot for spin-polarized scGW! "
          stop
        endif
      endif
      if(myid.eq.0 .and. i_cube .eq.1)then
        write(use_unit,*) " ...preparing scGW density for cube file... "
      endif
      densmat (:,:,i_cube) = green_fn(:,:)*(2.d0/n_spin)
    endif
  endif
!----------------------------------------------------
  enddo
  call de_kweight_occs('output_cube_files_p1', occ_numbers)


  ! First, determine the size of the arrays that have to be allocated

  n_compute_basis_local=0
  n_compute_basis=0

  ! JW, 2010-08-19:
  !
  ! In principle, the code of this subroutine is parallelized.  Unfortunately,
  ! the parallelization happens in an inner loop and the corresponding
  ! synchronization can be much more expensive than what is gained by doing
  ! the calculations in parallel.  The proper fix would be to change the
  ! parallelization.  But on the other hand, calculation times might depend
  ! mainly on the cube grid size and rather softly on system size, so that a
  ! serial implementation might do.
  !
  ! For now, use a wild guess about what might be better.
  !
  ! On the long run, it might be better to first calculate everything and then
  ! to synchronize.  Memory should not be an issue.

  n_syncs = product(cube_edge_steps(1:2, 1))
  cube_parallelized = (20 * n_centers_integrals * (n_basis/n_atoms) > n_syncs)

  if (cube_parallelized) then
     call localorb_info('  Dividing the cube on processors')
     n_points_per_task = int(ceiling(dble(cube_edge_steps(3,1))/dble(n_tasks)))
     off_z = myid*n_points_per_task
     last_z = min((myid+1) * n_points_per_task, cube_edge_steps(3,1))
  else
     call localorb_info('  Calculating whole cube on root processor')
     n_points_per_task = cube_edge_steps(3,1)
     if (myid == 0) then
        ! All points
        off_z = 0
        last_z = cube_edge_steps(3,1)
     else
        ! No points
        off_z = 0
        last_z = 0
     end if
  end if


  do coord_x =   1,cube_edge_steps(1,1),1
     do coord_y = 1,cube_edge_steps(2,1),1
        i_k = 0
        !      local_rho = 0.d0
        !      local_free_rho = 0.d0
        point_in_cell = .true.


        i_full_points = 0
        n_compute = 0
        i_basis = 0
        i_point = 0

        if (last_z > off_z) then
           do coord_z = off_z+1, last_z
              i_point = i_point+1
              !    generate output grid 
              coord_current(1)=cube_units(1) *(coord_x-1)-offset_coord(1)+cube_origin(1,1)
              coord_current(2)=cube_units(2) *(coord_y-1)-offset_coord(2)+cube_origin(2,1)
              coord_current(3)=cube_units(3) *(coord_z-1)-offset_coord(3)+cube_origin(3,1)



              if(n_periodic > 0)then
                 coord_current_temp = coord_current
                 call map_to_center_cell(coord_current_temp(1:3) )   
                 if( abs(coord_current(1) -coord_current_temp(1))> 1e-3 .or. abs(coord_current(2)-coord_current_temp(2))>1e-3 .or. &
                 abs(coord_current(3) - coord_current_temp(3))>1e-3 )then
                    !                             point_in_cell(i_point) = .false.
                 end if
                 coord_current = coord_current_temp
              end if

              if(point_in_cell(i_point))then

                 !     compute atom-centered coordinates of current integration point as viewed from all atoms
                 call tab_atom_centered_coords_p0 &
                 ( coord_current,  &
                 dist_tab_sq(1,i_point),  &
                 dir_tab(1,1,i_point), &
                 n_centers_basis_integrals, centers_basis_integrals )


                 call prune_basis_p0 &
                 (dist_tab_sq(1,i_point), &
                 n_compute_a, n_compute, i_basis,  &
                 n_centers_basis_T, n_centers_integrals, inv_centers_basis_integrals )


              end if ! point_in_cell 
           enddo     ! coord_z

           n_max_compute_dens = MAX(n_compute, n_max_compute_dens)
           n_compute_basis_local=MAX(n_compute,n_compute_basis_local)

           n_points = i_point
           if (n_compute.gt.0) then

              ! Determine all radial functions, ylm functions and their derivatives that
              ! are best evaluated strictly locally at each individual grid point.
              i_point = 0

              do coord_z = off_z+1, last_z

                 i_point = i_point+1
                 n_compute_atoms = 0
                 n_compute_fns = 0
                 i_basis_fns_inv = 0


                 if(point_in_cell(i_point))then

                    ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                    ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                    ! without any copying and without doing any unnecessary operations. 
                    ! The price is that the interface is no longer explicit in terms of physical 
                    ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.
                    call prune_radial_basis_p0 &
                    ( dist_tab_sq(1,i_point),  &
                    dist_tab(1), &
                    dir_tab(1,1,i_point), &
                    n_compute_atoms, atom_index, atom_index_inv, &
                    n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                    i_atom_fns, spline_array_start, spline_array_end, &
                    n_centers_integrals, centers_basis_integrals)


                    n_max_compute_fns_dens = MAX(n_compute_fns, n_max_compute_fns_dens)

                 end if ! point_in_cell

              enddo ! coord_z


           end if !n_compute
        end if

     enddo ! coord_y
  enddo ! coord_x 


  call sync_find_max(n_compute_basis_local,n_compute_basis)

  ! now allocate necessary quantities for this subroutine

  if(.not.allocated(radial_wave))then
     allocate(radial_wave(n_max_compute_fns_dens),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: radial_wave'
        stop
     end if
  else
     deallocate(radial_wave)
     allocate(radial_wave(n_max_compute_fns_dens),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: radial_wave'
        stop
     end if
  end if

  if(.not. allocated(wave))then
     allocate(wave(n_compute_basis,  cube_edge_steps(3,1)),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: wave'
        stop
     end if
  else
     deallocate(wave)
     allocate(wave(n_compute_basis,  cube_edge_steps(3,1)),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: wave'
        stop
     end if
  end if

  if(.not. allocated(density_matrix_con))then
     allocate(density_matrix_con(n_compute_basis,n_compute_basis),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: density_matrix_con'
        stop
     end if
  else
     deallocate(density_matrix_con)
     allocate(density_matrix_con(n_compute_basis,n_compute_basis),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: density_matrix_con'
        stop
     end if
  end if

  if(.not. allocated(work))then
     allocate(work(n_compute_basis,cube_edge_steps(3,1)),stat=info)
     call check_allocation(info, 'work            ')
  end if


  if(real_eigenvectors)then

     if(.not. allocated( KS_ev_compute))then
        allocate( KS_ev_compute(n_states, n_centers_basis_T, n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_ev_compute'
           stop
        end if
     end if

     if(.not. allocated( KS_vec))then
        allocate( KS_vec(n_centers_basis_T, n_states, n_cube),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_vec'
           stop
        end if
     end if
     KS_vec = 0.d0

  else

     if(.not. allocated( KS_vec_complex))then
        allocate( KS_vec_complex(n_centers_basis_T,n_states,n_cube),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation:  KS_vec_complex'
           stop
        end if
     end if
     KS_vec_complex = (0.d0, 0.d0)


     if(.not. allocated( KS_ev_compute_complex))then
        allocate( KS_ev_compute_complex(n_states, n_centers_basis_T, n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_ev_compute_complex'
           stop
        end if
     end if


  end if

  if(real_eigenvectors)then
     allocate( KS_orbital(1,cube_edge_steps(3,1),n_cube) ) 
  else
     allocate( KS_orbital_complex(1,cube_edge_steps(3,1),n_cube) ) 
  end if

  ! ************************************************
  ! ************************************************

  do i_cube=1,n_cube,1
     if (cube_type(i_cube).eq.'eigenstate'.or.(cube_type(i_cube).eq.'eigenstate_density'.and.cube_stm(3,i_cube).lt.0)) then
        i_k=int(ceiling(dble(cube_state(2,i_cube))/dble(n_tasks)))
        if(real_eigenvectors)then
           if (myid.eq.  MOD(cube_state(2,i_cube), n_tasks)) then

              do i_bas = 1,  n_centers_basis_T, 1

                 KS_vec(i_bas,cube_index(i_cube),i_cube) =  KS_eigenvector(Cbasis_to_basis(i_bas), &
                 cube_index(i_cube),cube_state(1,i_cube),i_k) !* &
                 !dble(k_phase(center_to_cell(Cbasis_to_center(  i_bas  )),cube_state(2,i_cube)))
              enddo
           endif

           call sync_vector(KS_vec(1,cube_index(i_cube),i_cube),n_centers_basis_T)

        else

           if (myid.eq.  MOD(cube_state(2,i_cube), n_tasks)) then

              do i_bas = 1,  n_centers_basis_T, 1

                 KS_vec_complex(i_bas,cube_index(i_cube),i_cube) =  KS_eigenvector_complex(Cbasis_to_basis(i_bas), &
                 cube_index(i_cube),cube_state(1,i_cube),i_k) !* &
                 !dconjg(k_phase(center_to_cell(Cbasis_to_center(  i_bas  )),cube_state(2,i_cube)))
              enddo

           endif

           call sync_vector_complex(KS_vec_complex(1,cube_index(i_cube),i_cube),n_centers_basis_T)

        endif
     endif
  end do

  do coord_x =   1,cube_edge_steps(1,1),1
     do coord_y = 1,cube_edge_steps(2,1),1

        local_rho = 0.d0
        local_rho_temp = 0.d0
        local_free_rho = 0.d0
        work = 0.0
        if(real_eigenvectors) then
           KS_orbital=0.0
        else
           KS_orbital_complex=(0.d0,0.d0)
        endif

        point_in_cell = .true.


        i_full_points = 0
        n_compute = 0
        i_basis = 0
        i_point = 0


        call get_timestamps(cpu, clock)

        if (myid.lt.cube_edge_steps(3,1)) then
           do coord_z = off_z+1, last_z

              i_point = i_point+1
              !    generate output grid 
              coord_current(1)=cube_units(1) *(coord_x-1)-offset_coord(1)+cube_origin(1,1)
              coord_current(2)=cube_units(2) *(coord_y-1)-offset_coord(2)+cube_origin(2,1)
              coord_current(3)=cube_units(3) *(coord_z-1)-offset_coord(3)+cube_origin(3,1)


              if(n_periodic > 0)then
                 coord_current_temp = coord_current
                 call map_to_center_cell(coord_current_temp(1:3) )   
                 if( abs(coord_current(1) -coord_current_temp(1))> 1e-3 .or. abs(coord_current(2)-coord_current_temp(2))>1e-3 .or. &
                 abs(coord_current(3) - coord_current_temp(3))>1e-3 )then
                    !                             point_in_cell(i_point) = .false.
                 end if
                 coord_current = coord_current_temp
              end if


              if(point_in_cell(i_point))then

                 !     compute atom-centered coordinates of current integration point as viewed from all atoms
                 call tab_atom_centered_coords_p0 &
                 ( coord_current,  &
                 dist_tab_sq(1,i_point),  &
                 dir_tab(1,1,i_point), &
                 n_centers_basis_integrals, centers_basis_integrals )


                 cube: do i_cube = 1, n_cube
                    if (cube_type(i_cube).eq.'delta_density') then

                       call tab_global_geometry_p0 &
                       ( dist_tab_sq(1,i_point), &
                       dir_tab(1,1,i_point), &
                       dist_tab(1), &
                       i_r(1), &
                       dir_tab_norm(1,1), &
                       n_centers_basis_integrals,  centers_basis_integrals)

                       call evaluate_free_rho_sums_p0 (  dist_tab(1), i_r(1), &
                       local_free_rho(off_z+i_point,1),  &
                       n_centers_basis_integrals, centers_basis_integrals)

                       local_free_rho(off_z+i_point,1) =  &
                       local_free_rho(off_z+i_point,1)* pi4_inv
                       exit cube 
                    end if !delta_density

                 end do cube
                 !     determine which basis functions are relevant at current grid point,
                 !     and tabulate their indices

                 !     next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed

                 call prune_basis_p0 &
                 (dist_tab_sq(1,i_point), &
                 n_compute_a, n_compute, i_basis,  &
                 n_centers_basis_T, n_centers_integrals, inv_centers_basis_integrals )

              end if ! point_in_cell 
           enddo     !          end loop over the z component

           n_points = i_point
           if (n_compute.gt.0) then

              ! Determine all radial functions, ylm functions and their derivatives that
              ! are best evaluated strictly locally at each individual grid point.
              i_point = 0

              do coord_z = off_z+1, last_z

                 i_point = i_point+1
                 n_compute_atoms = 0
                 n_compute_fns = 0
                 i_basis_fns_inv = 0


                 if(.not. point_in_cell(i_point))then

                    wave(:,i_point) = 0.d0
                 else

                    ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                    ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                    ! without any copying and without doing any unnecessary operations. 
                    ! The price is that the interface is no longer explicit in terms of physical 
                    ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.
                    call prune_radial_basis_p0 &
                    ( dist_tab_sq(1,i_point),  &
                    dist_tab(1), &
                    dir_tab(1,1,i_point), &
                    n_compute_atoms, atom_index, atom_index_inv, &
                    n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                    i_atom_fns, spline_array_start, spline_array_end, &
                    n_centers_integrals, centers_basis_integrals)


                    ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                    ! for all atoms which are actually relevant

                    call tab_local_geometry_p0 &
                    ( dist_tab_sq(1, i_point), n_compute_atoms, atom_index, &
                    dir_tab(1,1,i_point), dist_tab(1),  &
                    i_r(1) &
                    )


                    ! Now evaluate radial functions u(r) from the previously stored compressed spline arrays  

                    call evaluate_radial_functions_p0 &
                    (   spline_array_start, spline_array_end, &
                    n_compute_atoms, n_compute_fns,  &
                    dist_tab(1), i_r(1), &
                    atom_index, i_basis_fns_inv, &
                    basis_wave_ordered, radial_wave(1), &
                    .false., n_compute, n_max_compute_fns_dens   &
                    )

                    call tab_trigonom_p0 &
                    ( n_compute_atoms, dir_tab(1,1,i_point),  &
                    trigonom_tab(1,1) &
                    )

                    ! tabulate distance and Ylm's w.r.t. other atoms            
                    call tab_wave_ylm_p0 &
                    ( n_compute_atoms, atom_index,  &
                    trigonom_tab(1,1), l_shell_max,  &
                    l_ylm_max, &
                    ylm_tab(1,1) )


                    ! tabulate total wave function value for each basis function in all cases -
                    ! but only now are we sure that we have ylm_tab ...

                    call evaluate_waves_p0 &
                    (l_ylm_max, ylm_tab(1,1),  &
                    dist_tab(1),  &
                    index_lm, n_compute,  &
                    i_basis, radial_wave(1),  &
                    wave(1,i_point), n_compute_atoms,   &
                    atom_index_inv, n_compute_fns,  &
                    i_basis_fns_inv, n_max_compute_fns_dens  &
                    )
                 end if ! .not. point_in_cell

              enddo !        end loop over z component

           end if ! if (n_compute.gt.0)
           ! write(use_unit,*) 'free',local_free_rho(:,1)

           do i_cube = 1, n_cube, 1

              if (n_compute.gt.0) then

                 if (cube_type(i_cube).eq.'total_density' .or.cube_type(i_cube).eq.'delta_density'.or.&
                 (cube_type(i_cube).eq.'eigenstate_density'.and.cube_stm(3,i_cube).gt.0)) then

                    if(packed_matrix_format /= PM_none )then
                       call  prune_density_matrix_sparse(densmat_sparse(1,i_cube), density_matrix_con, &
                       n_compute, i_basis)
                    else
                       call  prune_density_matrix(densmat(1,1,i_cube), density_matrix_con, &
                       n_compute, i_basis)
                    end if

                    call evaluate_KS_density_densmat(n_points, wave(1,1), n_compute,   &
                    local_rho(off_z+1:last_z, 1, i_cube), &
                    n_compute_basis, n_centers_basis_T, &
                    density_matrix_con, work)

                 else 
                    if(cube_type(i_cube).eq.'eigenstate'.or.cube_type(i_cube).eq.'eigenstate_density')then
                       if(real_eigenvectors)then
                          call evaluate_KS_orbital_density_p0 &
                          (cube_index(i_cube),n_points, wave(1,1), n_compute, &
                          i_basis, KS_vec(1,1,i_cube),  &
                          1, &
                          KS_orbital(1, off_z+1:last_z, i_cube),  &
                          local_rho(off_z+1:last_z, 1, i_cube), &
                          n_compute_basis, n_centers_basis_T)
                       else
                          call evaluate_KS_orbital_density_complex_p0 &
                          (cube_index(i_cube),n_points, wave(1,1), n_compute, &
                          i_basis, KS_vec_complex(1,1,i_cube),  &
                          1, &
                          KS_orbital_complex(1, off_z+1:last_z, i_cube),  &
                          local_rho(off_z+1:last_z, 1, i_cube), &
                          n_compute_basis, n_centers_basis_T)

                       endif
                    endif
                 end if

              else

                 if (cube_type(i_cube).eq.'total_density' .or.cube_type(i_cube).eq.'delta_density'.or.&
                 cube_type(i_cube).eq.'eigenstate_density') then
                    local_rho(off_z+1:last_z,1,i_cube) = 0.d0
                 elseif (cube_type(i_cube).eq.'eigenstate')then
                    if(real_eigenvectors) then
                       KS_orbital(1,off_z+1:last_z, i_cube) = 0.0
                    else
                       KS_orbital_complex(1,off_z+1:last_z, i_cube)=(0.0,0.0)
                    endif
                 endif
              end if ! n_compute
           end do ! i_cube
        endif ! myid

        call get_times(cpu, clock, cpu_calc, clock_calc)

        if (cube_parallelized) then
           call get_timestamps(cpu, clock)
           call sync_vector(local_rho(1,1,1),cube_edge_steps(3,1)*n_cube*2)
           if (real_eigenvectors) then
              call sync_vector(KS_orbital(1,1,1),cube_edge_steps(3,1)*n_cube)
           else
              call sync_vector_complex(KS_orbital_complex(1,1,1),cube_edge_steps(3,1)*n_cube)
           endif
           call sync_vector(local_free_rho(1,1),cube_edge_steps(3,1)*2)
           call get_times(cpu, clock, cpu_sync, clock_sync)
        end if

        call get_timestamps(cpu, clock)

        if(myid ==0)then

           ! plot data
           do i_z = 1,cube_edge_steps(3,1),1

              do i_cube = 1, n_cube, 1


                 if (cube_type(i_cube) == 'total_density' .or. cube_type(i_cube) == 'eigenstate_density') then
                    if(cube_stm(3,i_cube).lt.0)then
                       write (10+i_cube,fmt=file_format) &
                       local_rho(i_z,1,i_cube) * inv_bohr_3
                    else
                       write (10+i_cube,fmt=file_format) &
                       hartree*cube_stm(2,i_cube)*local_rho(i_z,1,i_cube) * inv_bohr_3

                    endif
                 elseif (cube_type(i_cube)=='delta_density') then

                    write (10+i_cube,fmt=file_format) &
                    !                               local_rho_free(i_z,1,i_cube) * inv_bohr_3
                    (local_rho(i_z,1,i_cube)-local_free_rho(i_z,1)) * inv_bohr_3

                 else if(cube_type(i_cube) == 'eigenstate') then

                    if(real_eigenvectors) then
                       write (10+i_cube,fmt=file_format) &
                       KS_orbital(1,i_z,i_cube) * sqrt_inv_bohr_3
                    else
                       write (10+i_cube,fmt=file_format) &
                       Real(KS_orbital_complex(1,i_z,i_cube)) * sqrt_inv_bohr_3
                    endif

                 endif

                 if(.not.flag_gOpenMol_format .and. mod(i_z-1,6).eq.5) &
                 write (10+i_cube,*) ""

              enddo !i_cube
              if(if_stm) then
                 write(9,fmt=file_format) dble(i_z)
                 if(.not.flag_gOpenMol_format .and. mod(i_z-1,6).eq.5) &
                 write (9,*) ""
              endif
           enddo ! i_z
        endif
        if(myid.eq.0) then
          do i_cube = 1, n_cube, 1
              if(.not.flag_gOpenMol_format) write (10+i_cube,*) ""
           enddo
           if(if_stm) then
              if(.not.flag_gOpenMol_format) write (9,*) ""
           endif
        endif
        call get_times(cpu, clock, cpu_write, clock_write)

     end do            !     end loop over y component integration loop
  end do        !     end loop over x component
  !end do ! i_spin

  !     finally, deallocate stuff.
  if (allocated(ylm_tab)) then
     deallocate(ylm_tab)
  end if
  if (allocated(index_lm)) then
     deallocate(index_lm)
  end if
  if (allocated(dir_tab_global)) then
     deallocate(dir_tab_global)
  end if
  if (allocated(dist_tab_sq_global)) then
     deallocate(dist_tab_sq_global)
  end if
  if (allocated(local_rho)) then
     deallocate(local_rho)
  end if
  if (allocated(KS_orbital)) then
     deallocate(KS_orbital)
  end if
  if (allocated(KS_orbital_complex)) then
     deallocate(KS_orbital_complex)
  end if
  if (allocated(local_orb)) then
     deallocate(local_orb)
  end if


  do i_cube = 1, n_cube, 1
     close(10+i_cube)
  enddo

  if(myid.eq.0) then
     if(if_stm) close(9)
  endif

  if (allocated(point_in_cell)) deallocate(point_in_cell)

  if(allocated(radial_wave))deallocate(radial_wave)

  if(allocated(wave)) deallocate(wave)

  if(allocated(i_basis)) deallocate(i_basis)

  if(allocated( KS_ev_compute)) deallocate( KS_ev_compute )

  if(allocated( KS_ev_compute_complex)) deallocate( KS_ev_compute_complex)

  if(allocated(local_free_rho)) deallocate(local_free_rho) 

  if(allocated( KS_vec))deallocate( KS_vec)

  if(allocated( KS_vec_complex))deallocate( KS_vec_complex)

  if(allocated( densmat))deallocate( densmat)

  if(allocated( densmat_tmp))deallocate( densmat_tmp)

  if(allocated( densmat_sparse))deallocate( densmat_sparse)

  if(allocated( densmat_sparse_tmp))deallocate( densmat_sparse_tmp)

  if(allocated( work))deallocate( work)

  if(allocated( density_matrix_con))deallocate( density_matrix_con)

  call get_times(cpu_tot, clock_tot)
  call localorb_info('')
  call output_timeheader(deffmt, 'Cube output timings')
  call output_times(deffmt, 'Calculation for cube output', cpu_calc, clock_calc)
  if (cube_parallelized) then
     call output_times(deffmt, 'Synchronization for cube output', cpu_sync, clock_sync)
  end if
  call output_times(deffmt, 'I/O for cube output', cpu_write, clock_write)
  call output_times(deffmt, 'Total cube calculations and output', cpu_tot, clock_tot)

end subroutine output_cube_files_p1
!****** 
!!------------------------------------------------------------------
!!****s* FHI-aims/write_cube_header_gOpenMol
!!  NAME
!!    write_cube_header_gOpenMol
!!  SYNOPSIS
!
!subroutine write_cube_header_gOpenMol(descriptor, cube_units)
!
!!  PURPOSE
!!  Writes cube-files header for gOpenMol format
!!
!!  USES
!!
!  use dimensions
!  use runtime_choices
!  use grids
!  use geometry
!  use species_data
!  use mpi_utilities
!  use constants
!  use plot
!  implicit none
!
!!  ARGUMENTS
!
!  integer  descriptor
!  real*8 cube_units(3)
!
!!  INPUTS
!!   o descriptor -- file number where results are written
!!   o cube_units -- scale of the units
!!   
!!  OUTPUT
!!    none
!!  AUTHOR
!!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!!  SEE ALSO
!!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!!    Computer Physics Communications (2008), submitted.
!!  COPYRIGHT
!!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!!   the terms and conditions of the respective license agreement."
!!  HISTORY
!!    Release version, FHI-aims (2008).
!!  SOURCE
!!
!
!
!  !     local
!  integer i_atom
!
!
!
!
!  write (descriptor,*) '3  2'
!
!  write (descriptor,'(3I4)')   cube_edge_steps(1,1), cube_edge_steps(2,1), cube_edge_steps(3,1)
!
!  write (descriptor,fmt='(6F12.6)') &
!       (-cube_units(1)* cube_edge_steps(1,1)/2 + cube_origin(1,1))*bohr, &
!       (cube_units(1)* cube_edge_steps(1,1)/2 + cube_origin(1,1))*bohr, &
!       (-cube_units(2)* cube_edge_steps(2,1)/2 + cube_origin(2,1))*bohr, &
!       (cube_units(2)* cube_edge_steps(2,1)/2 + cube_origin(2,1))*bohr, &
!       (-cube_units(3)* cube_edge_steps(3,1)/2 + cube_origin(3,1))*bohr, &
!       (cube_units(3)* cube_edge_steps(3,1)/2 + cube_origin(3,1))*bohr
!
!  open(88,file='coords.xyz')
!
!  write(88,*) n_atoms
!  write(88,*) 'cell'
!
!  do i_atom = 1, n_atoms
!
!     write(88,'(A,3F12.6)') trim(species_name(species(i_atom))), coords(3,i_atom),coords(2,i_atom),coords(1,i_atom)
!
!  end do
!
!
!  close(88)
!
!end subroutine write_cube_header_gOpenMol
!!******
!!------------------------------------------------------------------
!
!
!
!
