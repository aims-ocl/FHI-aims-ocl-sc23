subroutine get_esp_charges
!  PURPOSE
!   Wrapper function to calculate ESP charges.
!   Calculates transition denstiy, multipole moments, potential finally fits
!   ESP charges to the potential
!  USES

  use pbc_lists,only: kweight_occs
  use runtime_choices,only:out_esp_full,esp_k_point,esp_min_radius,&
                           esp_max_radius,esp_n_max_radius,esp_grid_type,&
                           real_eigenvectors
  use dimensions,only: n_esp,n_atoms,n_states,n_spin,n_k_points,n_periodic,&
                       n_species,n_hamiltonian_matrix_size,n_full_points,&
                       l_pot_max,n_basis
  use mpi_tasks,only:myid,check_allocation,n_tasks
  use runtime_choices,only:use_density_matrix
  use geometry,only: cell_volume,lattice_vector
  use localorb_io,only: use_unit
  use physics,only:KS_eigenvalue, occ_numbers, KS_eigenvector, &
                   KS_eigenvector_complex, partition_tab,&
                   hartree_partition_tab,hamiltonian,overlap_matrix,&
                   hartree_potential,rho_change,nlcc_forces_on,gga_forces_on,&
                   pulay_forces_on,pulay_forces,nlcc_forces,gga_forces,&
                   hellman_feynman_forces,rho,rho_gradient,free_rho_superpos
  use species_data,only:l_shell_max
  use esp_charges
  use mixing,only:  delta_rho_KS,delta_rho_gradient
  use esp_grids,only:  grid_partitioned_esp,n_radial_esp,n_full_points_esp,&
                       cleanup_grids_esp,allocate_grids_esp,get_grids_esp,&
                       esp_get_n_compute_dens
  use constants, only: bohr
  use synchronize_mpi_basic, only: sync_vector,sync_vector_complex
  implicit none

  !  ARGUMENTS

  !  INPUTS
  !    KS-coefficents from SCF-calculation are used
  !  OUTPUT
  !    ESP-charges are written to output
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
  !



  !  local variables
  integer :: state_one, state_two, i_k_point_one, i_esp, &
             i_esp_n_radius,i_esp_n_radius_prev, i_spin_esp
  integer :: n_esp_internal
  logical :: density_full
  real*8  :: esp_min, esp_max
  real*8  :: esp_min_prev, esp_max_prev
  integer :: info,i_state,j_state
  real*8, dimension(:,:), allocatable :: hamiltonian_two
  real*8, dimension(:, :),allocatable :: rho_trans
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
  real*8, target, dimension(:) ,allocatable        :: potential_trans
  real*8, target, dimension(:) ,allocatable        :: partition_tab_esp
  real*8, dimension(:) ,allocatable        :: radius_esp_min
  real*8, dimension(:) ,allocatable        :: radius_esp_max
  real*8, dimension(:),allocatable :: free_hartree_superpos_trans
  real*8,     dimension(:,:,:) ,allocatable     :: occ_numbers_esp
  real*8,     dimension(:,:) ,allocatable     :: esp_charges_fit
  real*8,     dimension(:,:) ,allocatable     :: esp_dipole_fit
  integer,     dimension(:,:) ,allocatable     :: esp_state_full
  real*8,     dimension(:,:,:,:) ,allocatable     ::   KS_eigenvector_esp
  complex*16,     dimension(:,:,:,:) ,allocatable     ::   KS_eigenvector_complex_esp
  real*8,     dimension(3)     :: off
  integer :: esp_n_state_min, esp_n_state_max, rm, km
  real*8  :: R_c
  logical :: output_esp
  integer :: pbc_method
  integer :: i_esp_equal_grid_x, i_esp_equal_grid_y, i_esp_equal_grid_z,&
             i_esp_equal_grid_x_prev, i_esp_equal_grid_y_prev, &
             i_esp_equal_grid_z_prev, i_esp_grid, esp_grid_prev
  real*8, dimension(:),allocatable                     :: esp_cube_output
  real*8, dimension(:,:),allocatable                     :: esp_cube_coord
  integer                                    :: esp_output_cube_n = 0
  logical                                    :: esp_equal_grid = .false.
  logical                                    :: esp_log_grid = .false.
  logical                                    :: cube_rad_grid = .false.
  logical                                    :: out_pot = .false.
  logical                                    :: use_dip_for_cube = .false.
  logical                                    :: pot_with_dip =.false.
  if(out_esp_full)then
    i_k_point_one = esp_k_point
    esp_min = esp_min_radius
    esp_max =  esp_max_radius
    esp_min_prev = esp_min_radius
    esp_max_prev =  esp_max_radius
    density_full = .false.
    i_esp_n_radius= esp_n_max_radius
    i_esp_n_radius_prev= esp_n_max_radius
    output_esp = .false.
    esp_output_cube_n = 0
    out_pot = .false.
    rm = 7
    km = 7
    R_c = 10.0/bohr
    pbc_method = 2
    i_esp_grid = esp_grid_type
    esp_grid_prev = i_esp_grid
    if (i_esp_grid.eq.1)then
      esp_log_grid = .false.
      esp_equal_grid = .false.
    elseif (i_esp_grid.eq.2)then
      esp_log_grid = .true.
      esp_equal_grid = .false.
    else
      esp_log_grid = .false.
      esp_equal_grid = .true.
    endif
    i_esp_equal_grid_x = esp_n_max_radius
    i_esp_equal_grid_y = esp_n_max_radius
    i_esp_equal_grid_z = esp_n_max_radius
    i_esp_equal_grid_x_prev = i_esp_equal_grid_x
    i_esp_equal_grid_y_prev = i_esp_equal_grid_y
    i_esp_equal_grid_z_prev = i_esp_equal_grid_z
    call esp_state_minmax(KS_eigenvalue, esp_n_state_min, esp_n_state_max)
    n_esp_internal = (((esp_n_state_max-esp_n_state_min+1)+1)*&
         (esp_n_state_max-esp_n_state_min+1)/2)
    if(.not. allocated( esp_state_full))then
      allocate( esp_state_full(2,n_esp_internal),stat=info)
      if(info/=0)then
        write(use_unit,*)'Error in allocation: esp_state_full'
        stop
      endif
    end if
    i_esp = 0
    do i_state =  esp_n_state_min, esp_n_state_max, 1
      do j_state = i_state, esp_n_state_max, 1
        i_esp = i_esp + 1
        esp_state_full(1,i_esp) = i_state
        esp_state_full(2,i_esp) = j_state
      enddo
    enddo
  else
    n_esp_internal = n_esp
    output_esp = .true.
  endif

  if(.not. allocated( esp_charges_fit))then
    allocate( esp_charges_fit(n_atoms,n_esp_internal),stat=info)
    if(info/=0)then
      write(use_unit,*)'Error in allocation: esp_charges_fit'
      stop
    end if
  end if
  if(.not. allocated( esp_dipole_fit))then
    allocate( esp_dipole_fit(3,n_esp_internal),stat=info)
    if(info/=0)then
      write(use_unit,*)'Error in allocation: esp_dipole_fit'
      stop
    end if
  end if

  if(.not. allocated( occ_numbers_esp))then
    allocate( occ_numbers_esp(n_states, n_spin, n_k_points),stat=info)
    if(info/=0)then
      write(use_unit,*)'Error in allocation: occ_numbers_esp'
      stop
    end if
  end if
   ! call kweight_occs('esp_charges', occ_numbers_esp)
  do i_esp = 1, n_esp_internal, 1
    if(out_esp_full)then
      state_one = esp_state_full(1,i_esp)
      state_two = esp_state_full(2,i_esp)
    else
      state_one = esp_state(1,i_esp)
      state_two = esp_state(2,i_esp)
      i_k_point_one = esp_kpoint(i_esp)
      i_spin_esp = esp_spin(i_esp)
      esp_min = esp_vdw_radius(1,i_esp)
      esp_max = esp_vdw_radius(2,i_esp)
      density_full = flag_density_full(i_esp)
      i_esp_n_radius= esp_n_radius(i_esp)
      rm = esp_rm(i_esp)
      km = esp_km(i_esp)
      R_c = esp_R_c(i_esp)/bohr
      pbc_method = esp_pbc_method(i_esp)
      i_esp_equal_grid_x = esp_equal_grid_n(1,i_esp)
      i_esp_equal_grid_y = esp_equal_grid_n(2,i_esp)
      i_esp_equal_grid_z = esp_equal_grid_n(3,i_esp)
      esp_output_cube_n = esp_output_cube(i_esp)
      out_pot = esp_output_pot(i_esp)
      use_dip_for_cube = esp_use_dip_for_cube(i_esp)
      i_esp_grid = esp_grid(i_esp)
      if (i_esp.eq.1) then
        esp_min_prev = esp_vdw_radius(1,i_esp)
        esp_max_prev = esp_vdw_radius(2,i_esp)
        i_esp_n_radius_prev= esp_n_radius(i_esp)
        i_esp_equal_grid_x_prev = esp_equal_grid_n(1,i_esp)
        i_esp_equal_grid_y_prev = esp_equal_grid_n(2,i_esp)
        i_esp_equal_grid_z_prev = esp_equal_grid_n(3,i_esp)
        esp_grid_prev = esp_grid(i_esp)
      else
        esp_min_prev = esp_vdw_radius(1,i_esp-1)
        esp_max_prev = esp_vdw_radius(2,i_esp-1)
        i_esp_n_radius_prev= esp_n_radius(i_esp-1)
        i_esp_equal_grid_x_prev = esp_equal_grid_n(1,i_esp-1)
        i_esp_equal_grid_y_prev = esp_equal_grid_n(2,i_esp-1)
        i_esp_equal_grid_z_prev = esp_equal_grid_n(3,i_esp-1)
        esp_grid_prev = esp_grid(i_esp-1)
      endif
      if (i_esp_grid.eq.1)then
        esp_log_grid = .false.
        esp_equal_grid = .false.
        cube_rad_grid = .false.
      elseif (i_esp_grid.eq.2)then
        esp_log_grid = .true.
        esp_equal_grid = .false.
        cube_rad_grid = .false.
      elseif (i_esp_grid.eq.3)then
        esp_log_grid = .false.
        esp_equal_grid = .true.
        cube_rad_grid = .false.
      elseif(i_esp_grid.eq.4)then
        esp_log_grid = .false.
        esp_equal_grid = .true.
        cube_rad_grid = .true.
      endif
    endif
    if (myid.eq.0) then
      write (use_unit,'(2X,A)') ' '
      write (use_unit,'(2X,A)') "ESP charge fitting starts."
      write(use_unit,'(2X,A,1X,I8.8,1X,A,1X,I8.8)') "| Calculation : ", &
           i_esp, " of ", n_esp_internal
      if (density_full)then
        write (use_unit,'(2X,A)') "| Full density"
      else
        write(use_unit,'(2X,A,1X,I8.8,1X,I8.8)') &
             "| Transistion Eigenstate i -> j: ", state_one, state_two
        write(use_unit,'(2X,A,1X,I8.8)') "| Spin: ", i_spin_esp
        write(use_unit,'(2X,A,1X,I8.8)') "| k-point: ", i_k_point_one
      endif
    endif
  ! Allocations
    if(.not.density_full)then
      ! Big One: Actually only needed if state_one.ne.state_two
      if(.not.density_full.and.state_one.ne.state_two)then
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

    if(.not. allocated( rho_trans))then
      allocate( rho_trans(n_spin, n_full_points),stat=info)
        if(info/=0)then
          write(use_unit,*)'Error in allocation: rho_trans'
          stop
        end if
      end if
    endif

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



  ! Calculte requested density
    if (density_full) then
      occ_numbers_esp = occ_numbers
      call kweight_occs('esp_charges', occ_numbers_esp)
    else
      occ_numbers_esp = 0.d0
      occ_numbers_esp(state_one, 1, i_k_point_one) = 1.d0
      occ_numbers_esp(state_two, 1, i_k_point_one) = 1.d0
    endif

    if ((esp_output_cube_n.le.1).or.(esp_output_cube_n.eq.6))then
      if(.not.density_full)then
        if(use_density_matrix)then
          if (myid.eq.0) then
            write (use_unit,'(2X,A)') &
                 "Computing density with Density matrix."
          end if
          call get_transition_density_densmat &
               ( KS_eigenvector, KS_eigenvector_complex, occ_numbers_esp, &
                 partition_tab, hartree_partition_tab,  l_shell_max, &
                 rho_trans,  hamiltonian(1,1),hamiltonian_two(1,1),state_one,&
                 state_two,i_k_point_one,density_full)
        else
          if (myid.eq.0) then
            write (use_unit,'(2X,A)') "Computing density for Cluster case."
          end if
          if(n_tasks.gt.1.and.state_one.ne.state_two)then
            if(real_eigenvectors)then
              if (.not.allocated(KS_eigenvector_esp)) then
                allocate( KS_eigenvector_esp(n_basis, n_states, n_spin,n_k_points),&
                     stat=info)
                call check_allocation(info, 'KS_eigenvector_esp                  ')
              end if
              KS_eigenvector_esp=0.d0
              if(myid.eq.0)then
                KS_eigenvector_esp(:,:,:,:)=KS_eigenvector(:,:,:,:)
              endif
              KS_eigenvector(:,:,:,:)=0.d0
              if(myid.eq.0)then
                KS_eigenvector(:,:,:,:)=KS_eigenvector_esp(:,:,:,:)
              endif
              call sync_vector(KS_eigenvector,n_basis*n_states*n_spin*n_k_points)
              deallocate( KS_eigenvector_esp)
            else
              if (.not.allocated(KS_eigenvector_complex_esp)) then
                allocate( KS_eigenvector_complex_esp(n_basis, n_states, n_spin,n_k_points),&
                     stat=info)
                call check_allocation(info, 'KS_eigenvector_complex_esp                  ')
              end if
              KS_eigenvector_complex_esp=0.d0
              if(myid.eq.0)then
                KS_eigenvector_complex_esp(:,:,:,:)=KS_eigenvector_complex(:,:,:,:)
              endif
              KS_eigenvector_complex(:,:,:,:)=0.d0
              if(myid.eq.0)then
                KS_eigenvector_complex(:,:,:,:)=KS_eigenvector_complex_esp(:,:,:,:)
              endif
              call sync_vector_complex(KS_eigenvector_complex,n_basis*n_states*n_spin*n_k_points)
              deallocate( KS_eigenvector_complex_esp)
            end if
          endif

          call get_transition_density &
               ( KS_eigenvector, KS_eigenvector_complex, &
                 occ_numbers, &
                 partition_tab, hartree_partition_tab, &
                 l_shell_max, &
                 rho_trans, &
                 state_one, state_two, &
                 i_k_point_one, i_k_point_one,density_full)
        endif
      endif

    ! Expand density into multipole moments
      if(density_full)then
        if (myid.eq.0) then
          write (use_unit,'(2X,A)') &
               "Computing multipole expansion for total density."
        endif
        call multipole_expantion_rho &
             ( hartree_partition_tab, rho, &
               delta_v_hartree_part_at_zero_trans, &
               delta_v_hartree_deriv_l0_at_zero_trans, &
               multipole_moments_trans, multipole_radius_sq_trans, &
               l_hartree_max_far_distance_trans, &
               outer_potential_radius_trans, free_rho_superpos, density_full,&
               .false., i_spin_esp)
      else
        if (myid.eq.0) then
          write (use_unit,'(2X,A)') &
               "Computing multipole expansion for transition density."
        endif
        call multipole_expantion_rho &
             ( hartree_partition_tab, rho_trans, &
               delta_v_hartree_part_at_zero_trans, &
               delta_v_hartree_deriv_l0_at_zero_trans, &
               multipole_moments_trans, multipole_radius_sq_trans, &
               l_hartree_max_far_distance_trans, &
               outer_potential_radius_trans, free_rho_superpos, density_full,&
              .false., i_spin_esp)
      endif
    endif

    if((i_esp.eq.1).or.(esp_min.ne.esp_min_prev).or.(esp_max.ne.esp_max_prev)&
       .or.(i_esp_n_radius.ne.i_esp_n_radius_prev).or.&
       (i_esp_equal_grid_x.ne.i_esp_equal_grid_x_prev).or.&
       (i_esp_equal_grid_y.ne.i_esp_equal_grid_y_prev).or.&
       (i_esp_equal_grid_z.ne.i_esp_equal_grid_z_prev).or.&
       (esp_grid_prev.ne.i_esp_grid))then
      grid_partitioned_esp = .false.
      !call cleanup_grids_esp( )
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
      if (myid.eq.0) then
        write (use_unit,'(2X,A)') &
             "Generating points for esp charges."
      end if
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
    ! Create points for the calculation of the potential
      call esp_radius(radius_esp_min, radius_esp_max, esp_min, esp_max)
      if(.not.esp_equal_grid)then
        call allocate_grids_esp( )
        n_radial_esp(:) = i_esp_n_radius
        call get_grids_esp(radius_esp_min, radius_esp_max,.False.,esp_log_grid)
      endif
      off = lattice_vector(1,1:3)
      call esp_partition_grid( radius_esp_min, radius_esp_max, .False.,&
                               1,esp_equal_grid,cube_rad_grid,&
                               i_esp_equal_grid_x, i_esp_equal_grid_y, i_esp_equal_grid_z, &
                               .false., lattice_vector,off)
      if(.not. allocated( partition_tab_esp))then
        allocate( partition_tab_esp(n_full_points_esp),stat=info)
        if(info/=0)then
          write(use_unit,*)'Error in allocation: partition_tab_esp'
          stop
        end if
      end if
      if(.not. allocated( potential_trans))then
        allocate( potential_trans(n_full_points_esp),stat=info)
        if(info/=0)then
          write(use_unit,*)'Error in allocation: potential_trans'
          stop
        end if
      end if
      if(.not. allocated( free_hartree_superpos_trans))then
        allocate( free_hartree_superpos_trans(n_full_points_esp),stat=info)
        if(info/=0)then
          write(use_unit,*)'Error in allocation: free_hartree_superpos_trans'
          stop
        end if
      end if
      if(density_full)then
        call esp_grid_storage ( partition_tab_esp, &
                                free_hartree_superpos_trans,&
                                esp_equal_grid )
      endif
      if(esp_equal_grid)then
        if(n_periodic.eq.0)then
          partition_tab_esp = 1d0
        else
          partition_tab_esp(:) = (1d0/dble(i_esp_equal_grid_x*i_esp_equal_grid_y*&
                                  i_esp_equal_grid_z))*cell_volume
        endif
      endif
    endif
  ! Sum up the potential from multipole moments at created points
    if (((esp_output_cube_n.eq.1).or.(esp_output_cube_n.eq.6)).and.use_dip_for_cube) then
      pot_with_dip =.true.
    else
      pot_with_dip =.false.
    endif
    if ((esp_output_cube_n.le.1).or.(esp_output_cube_n.eq.6))then
      if (myid.eq.0) then
        write (use_unit,'(2X,A)') &
             "Summing up potential from multipole expanded  density."
      end if
      call sum_up_potential &
           ( delta_v_hartree_part_at_zero_trans, &
             multipole_moments_trans, &
             partition_tab_esp, potential_trans, &
             multipole_radius_sq_trans, &
             l_hartree_max_far_distance_trans,  &
             outer_potential_radius_trans, &
             free_hartree_superpos_trans, &
             density_full, output_esp,pot_with_dip)
    else
      if (myid.eq.0) then
        if (esp_output_cube_n.eq.2)then
          write (use_unit,'(2X,A)') "Calculating XC-Potential."
        elseif (esp_output_cube_n.eq.3)then
          write (use_unit,'(2X,A)') "Calculating X-Potential."
        elseif (esp_output_cube_n.eq.4)then
           write (use_unit,'(2X,A)') "Calculating C-Potential."
        elseif (esp_output_cube_n.eq.5)then
           write (use_unit,'(2X,A)') "Calculating density."
        endif
      end if
      call esp_get_n_compute_dens( partition_tab_esp )
      call get_vxc_esp &
           ( KS_eigenvector, KS_eigenvector_complex, occ_numbers_esp, &
             partition_tab_esp,  l_shell_max, &
             potential_trans,  hamiltonian(1,1),hamiltonian_two(1,1),state_one,&
             state_two,i_k_point_one,density_full, esp_output_cube_n)
    endif
    if((esp_output_cube_n.gt.0).and.esp_equal_grid.and..not.out_pot)then
      if(.not. allocated( esp_cube_output))then
        if (myid.eq.0) then
          write (use_unit,'(2X,A)') "Writting potential to cube file."
        end if
        allocate( esp_cube_output(i_esp_equal_grid_x*i_esp_equal_grid_y*&
             i_esp_equal_grid_z),stat=info)
        if(info/=0)then
          write(use_unit,*)'Error in allocation: esp_cube_output'
          stop
        end if
      end if
      if(.not. allocated( esp_cube_coord))then
        allocate( esp_cube_coord(i_esp_equal_grid_x*i_esp_equal_grid_y*&
             i_esp_equal_grid_z,3),stat=info)
        if(info/=0)then
          write(use_unit,*)'Error in allocation: esp_cube_coord'
          stop
        end if
      end if
      if (esp_output_cube_n.lt.6) then
        call esp_collect_pot(potential_trans,partition_tab_esp,esp_cube_output,&
             i_esp_equal_grid_x,i_esp_equal_grid_y, i_esp_equal_grid_z )
      elseif (esp_output_cube_n.eq.6) then
        call esp_collect_pot_Hu(potential_trans,partition_tab_esp,esp_cube_output,&
             i_esp_equal_grid_x,i_esp_equal_grid_y, i_esp_equal_grid_z,esp_cube_coord )
      end if
      call esp_output_cubefile(esp_cube_output,esp_cube_coord,i_esp_equal_grid_x,&
           i_esp_equal_grid_y,i_esp_equal_grid_z, radius_esp_max,i_esp,&
           cube_rad_grid,esp_output_cube_n)
    endif

    if((esp_output_cube_n.eq.0).or.out_pot)then
    ! and fit ESP charges
      if (myid.eq.0.and.output_esp) then
        write (use_unit,'(2X,A)') "Fitting esp charges to potential."
      end if
      if (n_periodic.eq.0) then
        call esp_fit(potential_trans,partition_tab_esp, &
             esp_charges_fit(1:n_atoms,i_esp), &
             esp_dipole_fit(1:3,i_esp), &
             density_full,output_esp,out_pot)
      else
        select case(pbc_method)
          case(1)
            call esp_fit_pbc(potential_trans,partition_tab_esp, &
                 esp_charges_fit(1:n_atoms,i_esp), &
                 esp_dipole_fit(1:3,i_esp), &
                 density_full,output_esp,rm,km,R_c,out_pot,&
                 use_dip_for_cube)
          case(2)
            call esp_fit_pbc_two(potential_trans,partition_tab_esp, &
                 esp_charges_fit(1:n_atoms,i_esp), &
                 esp_dipole_fit(1:3,i_esp), &
                 density_full,output_esp,rm,km,R_c,out_pot,&
                 use_dip_for_cube)
          case(3)
            call esp_fit_pbc_three(potential_trans,partition_tab_esp, &
                 esp_charges_fit(1:n_atoms,i_esp), &
                 esp_dipole_fit(1:3,i_esp), &
                 density_full,output_esp,R_c,out_pot,&
                 use_dip_for_cube)
          case default
            call esp_fit_pbc(potential_trans,partition_tab_esp, &
                 esp_charges_fit(1:n_atoms,i_esp), &
                 esp_dipole_fit(1:3,i_esp), &
                 density_full,output_esp,rm,km,R_c,out_pot,&
                 use_dip_for_cube)
        endselect
      endif
      if(out_pot.and.esp_equal_grid)then
        if (myid.eq.0) then
          write (use_unit,'(2X,A)') "Writting fitted potential to cube file."
        end if
        if(.not. allocated( esp_cube_output))then
          allocate( esp_cube_output(i_esp_equal_grid_x*i_esp_equal_grid_y*&
               i_esp_equal_grid_z),stat=info)
          if(info/=0)then
            write(use_unit,*)'Error in allocation: esp_cube_output'
            stop
          end if
        end if
        if(.not. allocated( esp_cube_coord))then
          allocate( esp_cube_coord(i_esp_equal_grid_x*i_esp_equal_grid_y*&
                                  i_esp_equal_grid_z,3),stat=info)
          if(info/=0)then
            write(use_unit,*)'Error in allocation: esp_cube_coord'
            stop
          end if
        end if
        call esp_collect_pot(potential_trans,partition_tab_esp,esp_cube_output,&
             i_esp_equal_grid_x,i_esp_equal_grid_y,&
             i_esp_equal_grid_z )
        call esp_output_cubefile(esp_cube_output,esp_cube_coord,i_esp_equal_grid_x,&
             i_esp_equal_grid_y,i_esp_equal_grid_z,&
             radius_esp_max,i_esp,cube_rad_grid,esp_output_cube_n)
      endif
    endif
  ! Deallocation
    if (allocated(rho_trans)) then
      deallocate(rho_trans)
    end if
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
    if (allocated(esp_cube_output)) then
      deallocate(esp_cube_output)
    end if
    if (allocated(esp_cube_coord)) then
      deallocate(esp_cube_coord)
    end if
    if (allocated(hamiltonian_two)) then
      deallocate( hamiltonian_two)
    end if
  enddo
   if(out_esp_full.and.(myid.eq.0))then
     call output_esp_charges(esp_charges_fit,esp_dipole_fit, &
                 KS_Eigenvalue, esp_k_point, esp_n_state_min, esp_n_state_max)
   endif
   if(.not.esp_equal_grid)then
     call cleanup_grids_esp( )
   endif
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
   if (allocated(occ_numbers_esp)) then
     deallocate(occ_numbers_esp)
   end if
   if (allocated(esp_charges_fit)) then
     deallocate(esp_charges_fit)
   end if
   if (allocated(esp_dipole_fit)) then
     deallocate(esp_dipole_fit)
   end if
   if (allocated(esp_state_full)) then
     deallocate(esp_state_full)
   end if
   if (myid.eq.0.and.output_esp) then
     write (use_unit,'(2X,A)')  "Done."
   end if
endsubroutine get_esp_charges
