!  FIXME: For optimum integrations, the present subroutine should also interpolate all
!         wave functions as splines onto the radial integration grid, and _then_ perform
!         all normalisation integrals on this grid. THEN, we would have a numerically
!         well-defined Hamiltonian matrix later on.
!         The present case will produce exact orthonormality on the log grid, but only
!         approximate orthonormality on the later radial integration grid.
!****s* FHI-aims/shrink_fixed_basis_phi_thresh
!  NAME
!   shrink_fixed_basis_phi_thresh
!  SYNOPSIS
subroutine shrink_fixed_large_basis ( )
!  PURPOSE
!   Subroutine shrink_fixed_large_basis:
!   * orthogonalizes extra basis functions to atomic ones, for each species.
!   * throws out all basis functions which are substantially identical to another one
!   * gathers basis functions in a compressed array for numerics
!   * implicitly sort all basis functions into blocks of identical l per species
!   * sorts those blocks according to basis-dependent cutoff potential width 
!   * labels all basis functions to remember their origin
!  USES
  use constants,       only : light_speed_sq
  use dimensions,      only : n_max_basis_fns, n_species, l_ext_max,  n_basis_fns, &
                              n_ext, n_ext_fns, n_max_grid, n_max_spline, n_states, &
                              use_basis_gradients, use_ext_basis, use_plus_u, &
                              calculate_all_eigenstates
  use runtime_choices, only : flag_rel, force_smooth_cutoff, n_empty, out_basis, &
                              REL_atomic_zora, REL_own, smooth_cutoff_limit, &
                              wave_threshold
  use grids,           only : n_grid, r_grid, r_radial, r_grid_inc
  use basis            ! TODO: Add
  use species_data     ! TODO: Add
  use spline,          only : val_spline_deriv, cubic_spline
  use free_atoms,      only : free_potential_spl
  use mpi_tasks,       only : myid, aims_stop_coll
  use localorb_io,     only : use_unit, localorb_info
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
  !  imported variables
  !     input
  !     most "input" is found in module species_data
  !     output
  !     all output is found in module basis
  !  local variables
  !     n_max : aux variables to determine loop counter maxima
  !     i_first_fn : First basis function for current atom and current angular momentum shell
  !     i_reject_wave : If the innermost maximum of a given wave function lies inside of
  !                  r_radial(i_reject_wave), then that wave function will not be used.  
  integer   :: n_max
  integer   :: i_first_fn, i_last_fn
  character :: l_shell_str
  integer   :: i_reject_wave
  real*8    :: r_inner_max(n_max_basis_fns)
  integer   :: function_index( n_species, n_max_basis_fns)
  !     Arrays for intermediate storage/handling of radial functions before spline
  !     basis_wave is compacted array of basis functions
  !     basis_deriv is compacted array of basis function derivatives (if needed)
  !     basis_kinetic : directly calculated kinetic energy term for all basis functions
  !                which are not free-atom like basis functions. [Treat free-atom like basis functions
  !                separately so we can subtract out the free-atom basis-defining potential later.]
  real*8, dimension(:,:), allocatable :: basis_wave
  real*8, dimension(:,:), allocatable :: basis_deriv
  real*8, dimension(:,:), allocatable :: basis_kinetic
  real*8, dimension(:,:), allocatable :: basis_kinetic_scaled_zora
  !     Arrays for a possible reindexing of all spline arrays
  integer :: i_species_index(n_max_basis_fns)
  integer :: perm_ext_fns(n_max_basis_fns)
  integer :: perm_ext_fns_inv(n_max_basis_fns)
  real*8  :: temp_radius(n_max_basis_fns)
  integer :: n_function, i_offset_spl, i_basis_2
  integer :: i_offset_spl_density
  integer :: i_offset
  integer :: species_offset
  integer :: outer_point
  !  counters
  integer :: i_species, i_basis, i_function, i_atom, i_shell
  integer :: i_m, i_l, i_grid, i_type, i_function2
  integer :: i_ionic, i_conf, i_hydro, i_gaussian
  integer :: i_atomic
  integer :: i_spline
  !  functions
  real*8    :: int_log_mesh
  character :: l_to_str
  real*8    :: get_inner_max
  integer,     dimension(n_species)       :: species_first_function, species_last_function   
  integer,     dimension(n_max_basis_fns) :: fn_species
  real*8,      dimension(n_max_basis_fns) :: fn_eigenval
  real*8,      dimension(n_max_basis_fns) :: fn_cutoff_radius
  character*8, dimension(n_max_basis_fns) :: fn_type
  integer,     dimension(n_max_basis_fns) :: fn_l
  integer,     dimension(n_max_basis_fns) :: fn_n
  character*16                            :: output_name
  real*8,      dimension(n_max_basis_fns) :: scalar_prod
  integer                                 :: i_prev
  real*8 :: write_radius, write_function, write_i_r
  ! Paula: these are for ordering the basis functions, so that the memory use in
  ! packing would be smaller.
  integer,    allocatable, dimension(:)   :: perm_ext, perm_ext_inv
  integer,    allocatable, dimension(:)   :: work
  character*8,allocatable, dimension(:)   :: work_c
  real*8,     allocatable, dimension(:,:) :: work_r2D
  real*8,     allocatable, dimension(:)   :: work_r
  real*8,     allocatable, dimension(:)   :: radius
  real*8 :: dummy(1)
  real*8 :: V_radial_deriv, basis_radial_deriv
  character*100 :: info_str
  
  ! for atom bsse: may need cleanup
  integer :: j_basis, basis_counter, empty_basis  

! DIRAC SOLVER:  ADD CODE TO HANDLE DIRAC ORBITALS

!  Sanity check
  if(.not. use_ext_basis) then
    write(use_unit,*) "ERROR: Something is wrong! We do not use large basis, but"
    write(use_unit,*) "       shrink_fixed_large_basis was called!!"
    write(use_unit,*) "       will stop now!"
    stop
  end if

  call localorb_info("Assembling extended basis for auxiliary basis.",use_unit,'(2X,A)')
  !     allocations - local variables for basis function storage
  !     basis_deriv is strictly _not_ needed unless use_basis_gradients is true.
  !                 Allocate it here anyway (perhaps unnecessarily) so that
  !                 orthonormalize_basis_fn can be called properly
  allocate ( basis_wave    (n_max_grid, n_max_basis_fns) )
  allocate ( basis_deriv   (n_max_grid, n_max_basis_fns) )
  allocate ( basis_kinetic (n_max_grid, n_max_basis_fns) )
  if (flag_rel .eq. REL_atomic_zora) then
     allocate( basis_kinetic_scaled_zora(n_max_grid,n_max_basis_fns))
  end if

  !     determine basis _functions_ first
  !     i_function counts distinct radial functions (at, type, n, l)
  i_function             = 0
  i_basis                = 0
  n_ext                  = 0
  fn_eigenval            = 0.d0
  function_index         = 0
  species_first_function = 0
  species_last_function  = 0
  r_inner_max(:)         = 0d0
  do i_species = 1, n_species, 1
     ! species dependent counters & such
     species_first_function(i_species) = i_function + 1           ! remember where in the basis function list a certain species starts.
     i_reject_wave                     = innermost_max(i_species) ! index for wave function rejection
     i_basis                           = 0                        ! radial function counter for species
     fn_cutoff_radius                  = r_cutoff(i_species)      ! initialize cutoff radius for potential sorting later: this is the biggest value
     ! reat each l shell as a block
     do i_l = 0, ext_l_shell_max(i_species), 1
        i_type = 0
        i_first_fn = i_function+1    ! first function in a given l-shell of a given species, need to remember that throughout this loop!
        if (include_min_basis(i_species)) then
           ! treat atomic-like wave functions first
           i_type = i_type + 1
           do i_atomic = 1, n_atomic(i_species), 1
              if (atomic_l(i_species,i_atomic).eq.i_l) then
                 i_function                   = i_function+1
                 fn_species(i_function)       = i_species
                 fn_type(i_function)          = "atomic"
                 fn_l(i_function)             = i_l
                 fn_n(i_function)             = atomic_n(i_species,i_atomic)
                 fn_cutoff_radius(i_function) = atomic_outer_radius(i_species,i_atomic)
                 if (core_n_max(i_l,i_species) .ge. &
                      atomic_n(i_species,i_atomic)) then
                    ! only for our own relativity, store the correct eigenvalue -
                    ! for all but core functions, this has to be zero.
                    fn_eigenval(i_function) = &
                         atomic_eigenval(i_species,i_atomic)
                 end if
                 !               store current basis function in basis_wave
                 !               also store kinetic energy term of basis function
                 do i_grid = 1, n_grid(i_species),1
                    basis_wave(i_grid, i_function) = &
                         atomic_wave(i_grid, i_species, i_atomic)
                    basis_kinetic(i_grid, i_function) = &
                         atomic_kinetic(i_grid, i_species, i_atomic)
                 enddo
                 if (use_basis_gradients) then
                    do i_grid = 1, n_grid(i_species),1
                       basis_deriv(i_grid, i_function) = &
                            atomic_wave_deriv(i_grid, i_species, i_atomic)
                    enddo
                 end if
                 r_inner_max(i_function) = get_inner_max ( basis_wave(1,i_function), &
                      r_grid(1,i_species), n_grid(i_species) )

                 !               if a minimal basis function does not meet the criterion, stop the calculation
                 !               entirely - the grids must be chosen differently.
                 if (r_inner_max(i_function).lt.r_radial(i_reject_wave,i_species)) &
                      then
                    l_shell_str = l_to_str(i_l)
                    if (myid.eq.0) then
                       write(use_unit,'(1X,A,A)') &
                            "* Warning - the innermost maximum ", &
                            "of the minimal basis function"
                       write(use_unit,'(1X,A,I2,3A)') &
                            "* ", atomic_n(i_species,i_atomic), &
                            l_shell_str, &
                            " of species ",trim(species_name(i_species))
                       write(use_unit,'(1X,A,I1,A)') &
                            "* has its innermost maximum inside the ", &
                            i_reject_wave, "th radial integration shell."
                       write(use_unit,'(1X,A,A)') &
                            "* Adjust your integration grid to be ", &
                            "accurate enough before continuing."
                    end if
                    stop
                 end if
              end if
           enddo
        end if
        
        !     treat confined basis functions next
        if (n_conf(i_species).gt.0) then
           i_type = i_type+1
           do i_conf = 1, n_conf(i_species), 1
              if (conf_l(i_species,i_conf).eq.i_l) then
                 i_function                   = i_function+1
                 i_shell                      = conf_n(i_species, i_conf)
                 fn_species(i_function)       = i_species
                 fn_type(i_function)          = "confined"
                 fn_l(i_function)             = i_l
                 fn_n(i_function)             = i_shell
                 fn_cutoff_radius(i_function) = confined_outer_radius(i_species, i_conf)
                 if (core_n_max(i_l,i_species) .ge. i_shell) then
                    ! only for our own relativity, store the correct eigenvalue -
                    ! for all but core functions, this has to be zero.
                    fn_eigenval(i_function) = &
                         confined_eigenval(i_species,i_conf)
                 end if
                 ! copy current basis function to basis_wave
                 do i_grid = 1, n_grid(i_species), 1
                    basis_wave(i_grid, i_function)    = confined_wave(i_grid, i_species, i_conf)
                    basis_kinetic(i_grid, i_function) = confined_kinetic(i_grid, i_species, i_conf)
                 enddo
                 if (use_basis_gradients) then
                    do i_grid = 1, n_grid(i_species), 1
                       basis_deriv(i_grid, i_function) = &
                            confined_wave_deriv(i_grid, i_species, i_conf)
                    enddo
                 end if
                 r_inner_max(i_function) = get_inner_max ( basis_wave(1,i_function), &
                      r_grid(1,i_species), n_grid(i_species) )
              end if
           enddo
        end if           ! end confined function
        
        ! treat ionic basis functions next
        if (n_ionic(i_species).gt.0) then
           i_type = i_type+1
           do i_ionic = 1, n_ionic(i_species), 1
              if (ionic_l(i_species, i_ionic).eq.i_l) then
                 i_function                   = i_function+1
                 i_shell                      = ionic_n(i_species, i_ionic)
                 fn_species(i_function)       = i_species
                 fn_type(i_function)          = "ionic"
                 fn_l(i_function)             = i_l
                 fn_n(i_function)             = i_shell
                 fn_cutoff_radius(i_function) = ionic_outer_radius(i_species, i_ionic)
                 if (core_n_max(i_l,i_species) .ge. &
                      ionic_n(i_species,i_ionic)) then
                    ! only for our own relativity, store the correct eigenvalue -
                    ! for all but core functions, this has to be zero.
                    fn_eigenval(i_function) = &
                         ionic_eigenval(i_species,i_ionic)
                 end if
                 !               copy current basis function to basis_wave
                 do i_grid = 1, n_grid(i_species), 1
                    basis_wave(i_grid, i_function) = &
                         ionic_wave(i_grid, i_species, i_ionic)
                    basis_kinetic(i_grid, i_function) = &
                         ionic_kinetic(i_grid, i_species, i_ionic)
                 enddo
                 if (use_basis_gradients) then
                    do i_grid = 1, n_grid(i_species), 1
                       basis_deriv(i_grid, i_function) = &
                            ionic_wave_deriv(i_grid, i_species, i_ionic)
                    enddo
                 end if
                 r_inner_max(i_function) = get_inner_max ( basis_wave(1,i_function), &
                      r_grid(1,i_species), n_grid(i_species) )
              end if
           enddo
           !     end ionic function
        end if
        
        !     treat hydrogenic basis functions next
        !     use all hyrogenic functions both normal and auxliary
        if (n_hydro(i_species).gt.0 ) then
           i_type = i_type+1
           do i_hydro = 1, n_hydro(i_species), 1
              if (hydro_l(i_species, i_hydro).eq.i_l) then
                 i_function                   = i_function+1
                 i_shell                      = hydro_n(i_species, i_hydro)
                 fn_species(i_function)       = i_species
                 fn_type(i_function)          = "hydro"
                 fn_l(i_function)             = i_l
                 fn_n(i_function)             = i_shell
                 fn_cutoff_radius(i_function) = hydro_outer_radius(i_species,i_hydro)
                 !     copy current basis function to basis_wave
                 do i_grid = 1, n_grid(i_species), 1
                    basis_wave(i_grid, i_function) = &
                         hydro_wave(i_grid, i_species, i_hydro)
                    basis_kinetic(i_grid, i_function) = &
                         hydro_kinetic(i_grid, i_species, i_hydro)
                 enddo
                 if (use_basis_gradients) then
                    do i_grid = 1, n_grid(i_species), 1
                       basis_deriv(i_grid, i_function) = &
                            hydro_wave_deriv(i_grid, i_species, i_hydro)
                    enddo
                 end if
                 r_inner_max(i_function) = get_inner_max ( basis_wave(1,i_function), &
                      r_grid(1,i_species), n_grid(i_species) )
              end if
           enddo
           !         end hydrogenic function
        end if
        
        !         treat Gaussian basis functions next
        if (n_gaussian(i_species).gt.0) then
           i_type = i_type+1
           do i_gaussian = 1, n_gaussian(i_species), 1
              if (gaussian_l(i_species,i_gaussian).eq.i_l) then
                 i_function             = i_function + 1
                 i_shell                = gaussian_n(i_species, i_gaussian)
                 fn_species(i_function) = i_species
                 fn_type(i_function)    = "gaussian"
                 fn_l(i_function)       = i_l
                 fn_n(i_function)       = i_shell
                 fn_cutoff_radius(i_function) = gaussian_outer_radius(i_species,i_gaussian)
                 !               copy current basis function to basis_wave
                 do i_grid = 1, n_grid(i_species), 1
                    basis_wave(i_grid, i_function) = &
                         gaussian_wave(i_grid, i_species, i_gaussian)
                    basis_kinetic(i_grid, i_function) = &
                         gaussian_kinetic(i_grid, i_species, i_gaussian)
                 enddo
                 if (use_basis_gradients) then
                    do i_grid = 1, n_grid(i_species), 1
                       basis_deriv(i_grid, i_function) = &
                            gaussian_wave_deriv(i_grid, i_species, i_gaussian)
                    enddo
                 end if
                 r_inner_max(i_function) = get_inner_max ( basis_wave(1,i_function), &
                      r_grid(1,i_species), n_grid(i_species) )
              end if
           enddo
           !         end Gaussian function
        end if
        
        ! now have the basis functions and their associated information sorted into different arrays. 
        ! sort the cutoff potentials for each l-shell and then normalize according to cutoff!
        ! i_first_fn indexes the first function in a certain l-shell
        ! i_function indexes the last function in a given l-shell 
        i_last_fn  = i_function
        i_function = i_first_fn

        if (use_plus_u) then
           if (basis_dep_cutoff_thresh(i_species) .ne. 0d0) then
              basis_dep_cutoff_thresh(i_species) = 0d0
              write(info_str,'(1X,A)') "* Warning - Since '+U' treatment was requested,"
              call localorb_info(info_str)
              write(info_str,'(1X,A)') "* the basis-dependent cutoff is turned off automatically."
              call localorb_info(info_str)
           end if
        end if

!test
!        if (myid.eq.0) then
!          write(use_unit,*) 
!          write(use_unit,'(2X,A,I3)')  "l channel: ", i_l
!          write(use_unit,*) 
!          do i_function2 = i_function,i_last_fn,1
!            write(use_unit,'(2X,I5,F15.8,A)') i_function2, fn_cutoff_radius(i_function2)*bohr, " AA"
!          end do
!        end if 
!end test

        if (basis_dep_cutoff_thresh(i_species).gt.0d0) then
           ! allocate temporary storage arrays and switches if not already done so
           if (.not.allocated(perm_ext))     allocate(perm_ext    (n_max_basis_fns))
           if (.not.allocated(perm_ext_inv)) allocate(perm_ext_inv(n_max_basis_fns))
           if (.not.allocated(work))           allocate(work          (n_max_basis_fns))             ! integer
           if (.not.allocated(work_r))         allocate(work_r        (n_max_basis_fns))             ! real
           if (.not.allocated(work_c))         allocate(work_c        (n_max_basis_fns))             ! character
           if (.not.allocated(work_r2D))       allocate(work_r2D      (n_max_grid, n_max_basis_fns)) ! basis function
           ! sort single l-shell according to array fn_cutoff_radius: from i_function to i_last_fn         
           call insertionsort(fn_cutoff_radius(i_function:i_last_fn),i_last_fn-i_function+1, &
                perm_ext(i_function:i_last_fn), perm_ext_inv(i_function:i_last_fn) )
           ! all entries in perm_ext start counting @ 1 - but I would really like them to start at i_function. 
           perm_ext_inv(i_function:i_last_fn) = perm_ext_inv(i_function:i_last_fn) + i_function - 1 
           ! permute all other necessary data according to sorting:
           do i_function2 = i_function, i_last_fn
              work      (i_function2) = fn_n        (perm_ext_inv(i_function2))
              work_r    (i_function2) = r_inner_max (perm_ext_inv(i_function2))
              work_c    (i_function2) = fn_type     (perm_ext_inv(i_function2))
              work_r2D(:,i_function2) = basis_wave(:,perm_ext_inv(i_function2)) 
           end do
           fn_n        (i_function:i_last_fn) = work      (i_function:i_last_fn)
           fn_type     (i_function:i_last_fn) = work_c    (i_function:i_last_fn)
           r_inner_max (i_function:i_last_fn) = work_r    (i_function:i_last_fn)
           basis_wave(:,i_function:i_last_fn) = work_r2D(:,i_function:i_last_fn)
           do i_function2 = i_function, i_last_fn
              work_r2D(:,i_function2) = basis_deriv(:,perm_ext_inv(i_function2))
           end do
           basis_deriv(:,i_function:i_last_fn) = work_r2D(:,i_function:i_last_fn)
           do i_function2 = i_function, i_last_fn
              work_r2D(:,i_function2) = basis_kinetic(:,perm_ext_inv(i_function2))
           end do
           basis_kinetic(:,i_function:i_last_fn) = work_r2D(:,i_function:i_last_fn)
        end if

        ! keep for output in this particular l-shell
        l_shell_str = l_to_str(i_l)
        
        ! loop through all basis functions in the current l-shell, to orthonormalize the bunch according to current sorting 
        ! take special attention to the fact that their number might be reduced as we reject one here or there. 
        i_function = i_first_fn

        do while (i_function.le.i_last_fn)
           ! check inner radius
           if (r_inner_max(i_function).ge.r_radial(i_reject_wave,i_species)) then
              ! if alright, orthonormalize with the rest
              call orthonormalize_basis_fn &
                   ( i_function, i_first_fn, basis_wave, &
                   n_grid(i_species), &
                   r_grid(1,i_species), basis_acc(i_species), &
                   basis_kinetic, basis_deriv, &
                   function_index(i_species, i_function), .false., dummy)
           else
              function_index(i_species,i_function) = -1   ! if not, reject function
           end if  
           ! evaluate function according to the outcome of the orthornomalization procedure
           if ( function_index(i_species,i_function).gt.0) then
              ! accepted, that's good I suppose
              write(info_str,'(2X,3A,A8,A,I3,3A)') &
                   "| Species ", trim(species_name(i_species)), " : ", trim(fn_type(i_function)), &
                   " orbital ", fn_n(i_function)," ",l_shell_str," accepted."
              call localorb_info(info_str,use_unit,'(A)')                   
              i_basis    = i_basis + (2*i_l+1)
              i_function = i_function + 1
           else 
              if ( function_index(i_species,i_function).eq.0 ) then
                 write(info_str,'(2X,3A,A8,A,I3,1X,A,A)') &
                      "| Species ", trim(species_name(i_species)), " : ", trim(fn_type(i_function)), &
                      " orbital ", fn_n(i_function), l_shell_str, &
                      " rejected: Linear dependence."
                 call localorb_info(info_str,use_unit,'(A)')
              else
                 write(info_str,'(2X,3A,A8,A,I3,1X,A,A)') &
                      "| Species ", trim(species_name(i_species)), " : ",trim(fn_type(i_function)),&
                      " orbital ", fn_n(i_function), l_shell_str, &
                      " rejected: Insufficient integration grid."
                 call localorb_info(info_str,use_unit,'(A)')
              end if
              ! Function got rejected outright. 
              ! if this was an atomic function, there is something seriously wrong. Notify user and quit.
              if (fn_type(i_function).eq.'atomic') then
                 call localorb_info("* WARNING: Attempting to reject atomic basis function.",use_unit,'(A)')
                 call localorb_info("*          You may be running with a very large basis set, ",use_unit,'(A)')
                 call localorb_info("*          but non-zero basis_dep_cutoff in control.in.",use_unit,'(A)')
                 call localorb_info("*          If so, set basis_dep_cutoff to zero in control.in and retry.",use_unit,'(A)')
                 call localorb_info("*          Please also consult the manual on basis_dep_cutoff .",use_unit,'(A)')
                 call localorb_info("*          Check settings and restart.",use_unit,'(A)')
                 stop
              end if
              ! Before continuing, sort every other function in this
              ! l-shell back one step. In particular, this requires moving ...
              !    basis_wave(n_max_grid, n_max_basis_fns)
              !    basis_kinetic(n_max_grid, n_max_basis_fns)
              !    basis_deriv(n_max_grid, n_max_basis_fns)
              do i_function2 = i_function + 1, i_last_fn
                 basis_wave   (:,i_function2-1) = basis_wave   (:,i_function2)
                 basis_kinetic(:,i_function2-1) = basis_kinetic(:,i_function2)
                 basis_deriv  (:,i_function2-1) = basis_deriv  (:,i_function2)
                 r_inner_max  (  i_function2-1) = r_inner_max  (  i_function2)
                 fn_n         (  i_function2-1) = fn_n         (  i_function2)
                 fn_type      (  i_function2-1) = fn_type      (  i_function2)
                 fn_l         (  i_function2-1) = fn_l         (  i_function2)
                 fn_species   (  i_function2-1) = fn_species   (  i_function2)
              end do
              ! finally, the last function must be decreased by 1 as we took one out in the middle. 
              i_last_fn = i_last_fn - 1
           end if
        end do     ! loop over basis functions within a given l-shell
        i_function = i_function - 1  ! This points to the last known accepted function. 
     enddo      ! end loop over l-shells
     species_last_function(i_species) = i_function    ! remember where a species stops, required for labeling later.
     n_ext = n_ext + i_basis * atoms_in_structure(i_species)
  enddo    ! first loop over species 

  ! deallocate work array from function sorting within a single l-shell
  if (allocated(perm_ext))     deallocate(perm_ext)
  if (allocated(perm_ext_inv)) deallocate(perm_ext_inv)
  if (allocated(work))           deallocate(work)
  if (allocated(work_r))         deallocate(work_r)
  if (allocated(work_c))         deallocate(work_c)
  if (allocated(work_r2D))       deallocate(work_r2D)

  !     now we know the actual basis size, and hence all necessary array dimensions
  n_ext_fns = i_function

  !     verify n_states against n_basis
  if (n_states.gt.n_ext) then
     if (myid.eq.0) then
        if (n_empty /= -1 .or. calculate_all_eigenstates) then    ! explicitly requested
           write(use_unit,"()")
           write(use_unit,"(2X,A,I8,A)") "* You requested ", n_states, &
           & " Kohn-Sham states in the calculation,"
           write(use_unit,"(2X,A,I8,A)") "* but there are only ", n_ext, &
           & " available basis states."
           write(use_unit,"(2X,A,A,I8,A)") "* Reducing total number of ",&
           " Kohn-Sham states to ", n_ext, "."
           write(use_unit,*)
        else
           write(use_unit,"(2X,A,A,I8,A)") "Reducing total number of ",&
           " Kohn-Sham states to ", n_ext, "."
        end if
     end if
!     n_states = n_ext
  end if

  ! VB: Deal with specialities for kinetic energy here. If atomic_zora needed, must amend
  !     basis_kinetic by various extra terms.
  !
  if (flag_rel .eq. REL_atomic_zora)then
     do i_function = 1, n_ext_fns, 1
        do i_grid = 1, n_grid(fn_species(i_function)),1
           ! first, obtain radial derivative of atomic potential
           V_radial_deriv =  val_spline_deriv( dble(i_grid), &
                free_potential_spl(1,1, fn_species(i_function)), &
                n_grid(     fn_species(i_function)          ) ) &
                / (log(r_grid_inc( fn_species(i_function))) &
                * r_grid(i_grid, fn_species(i_function)))
           ! ... and radial derivative of this radial function
           basis_radial_deriv = basis_deriv(i_grid, i_function)
           ! first, use unmodified kinetic energy part to create renormalization
           ! expression for "scaled ZORA"
           basis_kinetic_scaled_zora(i_grid,i_function) = &
                basis_kinetic(i_grid,i_function) &
                * 2 * light_speed_sq &
                / (  2 * light_speed_sq &
                - free_potential_spl(1,i_grid,fn_species(i_function)))**2 &
                !
                -  2*light_speed_sq &
                !
                * V_radial_deriv &
                * (   basis_radial_deriv &
                -    basis_wave(i_grid,i_function) &
                / r_grid(i_grid, fn_species(i_function))) &
                !
                / (  2*light_speed_sq &
                - free_potential_spl(1,i_grid,fn_species(i_function)))**3
           
           ! now, modify the kinetic energy part itself for on-site ("atomic") ZORA corrections
           basis_kinetic(i_grid,i_function) = &
                basis_kinetic(i_grid,i_function) &
                * 2 * light_speed_sq &
                / (  2 * light_speed_sq &
                - free_potential_spl(1,i_grid,fn_species(i_function))) &
                !
                - light_speed_sq &
                / (  2*light_speed_sq &
                - free_potential_spl(1,i_grid,fn_species(i_function)))**2 &
                !
                * V_radial_deriv &
                * (   basis_radial_deriv &
                -    basis_wave(i_grid,i_function) &
                / r_grid(i_grid, fn_species(i_function)))
           
        end do
     end do
     
  end if ! REL_atomic_zora
  
  if (flag_rel .eq. REL_own) then
     do i_function = 1, n_ext_fns, 1
        ! amend kinetic energy by core eigenvalue (all other eigenvales are set to zero above)
        
        ! Notice that this treatment is still outrightly wrong. Really, we should be
        ! setting the correct kinetic energy term already when the basis functions are
        ! generated. Then, we could simply use the orthonormalizaton above to project out all core
        ! contributions correctly. Now, confined and ionic core radial functions will
        ! get the wrong kinetic energy ... as will any accidentally overlapping core contributions from
        ! the original non-orthonormalized basis functions.
        
        do i_grid = 1, n_grid(fn_species(i_function)),1
           ! first, obtain radial derivative of atomic potential
           V_radial_deriv =  val_spline_deriv( dble(i_grid), &
                free_potential_spl(1,1, fn_species(i_function)), &
                n_grid(     fn_species(i_function)          ) ) &
                / (log(r_grid_inc( fn_species(i_function))) &
                * r_grid(i_grid, fn_species(i_function)))
           ! ... and radial derivative of this radial function
           basis_radial_deriv = basis_deriv(i_grid, i_function)
           ! now, modify the kinetic energy part itself for on-site ("atomic") ZORA corrections
           basis_kinetic(i_grid,i_function) = &
                basis_kinetic(i_grid,i_function) &
                * 2 * light_speed_sq &
                / (  2 * light_speed_sq + fn_eigenval(i_function) &
                - free_potential_spl(1,i_grid,fn_species(i_function))) &
                !
                - light_speed_sq &
                / (  2*light_speed_sq + fn_eigenval(i_function) &
                - free_potential_spl(1,i_grid,fn_species(i_function)))**2 &
                !
                * V_radial_deriv &
                * (   basis_radial_deriv &
                -    basis_wave(i_grid,i_function) &
                / r_grid(i_grid, fn_species(i_function)))
        end do
     end do
  end if ! REL_own
  
  if (out_basis) then
     !       Output accepted basis functions for consistency
     do i_function = 1, n_ext_fns, 1
        l_shell_str = l_to_str(fn_l(i_function))
        if (i_function.lt.10) then
           if (fn_n(i_function).lt.10) then
              write(output_name,'(A2,A1,I1,A1,I1,A1,I1,A1,A1,A4)') &
                   fn_type(i_function),"_",fn_species(i_function),"_", &
                   i_function, "_", &
                   fn_n(i_function),"_",l_shell_str,".dat"
           else if  (fn_n(i_function).lt.100) then
              write(output_name,'(A2,A1,I1,A1,I1,A1,I2,A1,A1,A4)') &
                   fn_type(i_function),"_",fn_species(i_function),"_", &
                   i_function, "_", &
                   fn_n(i_function),"_",l_shell_str,".dat"
           end if
        else if (i_function.lt.100) then
           if (fn_n(i_function).lt.10) then
              write(output_name,'(A2,A1,I1,A1,I2,A1,I1,A1,A1,A4)') &
                   fn_type(i_function),"_",fn_species(i_function),"_", &
                   i_function, "_", &
                   fn_n(i_function),"_",l_shell_str,".dat"
           else if  (fn_n(i_function).lt.100) then
              write(output_name,'(A2,A1,I1,A1,I2,A1,I2,A1,A1,A4)') &
                   fn_type(i_function),"_",fn_species(i_function),"_", &
                   i_function, "_", &
                   fn_n(i_function),"_",l_shell_str,".dat"
           end if
        end if
        
        if (myid.eq.0) then
           open(50, file=output_name)
           write(50,*) "# ",n_grid(fn_species(i_function))
           do i_grid = 1, n_grid(fn_species(i_function)), 1
              write(50,*) &
                   r_grid(i_grid, fn_species(i_function)), &
                   basis_wave(i_grid, i_function)
           enddo
           close(50)
        end if
        !         output kinetic energy term also
        write(output_name,'(A4,A2,A1,I1,A1,I1,A1,A1,A4)') &
             "kin_", &
             fn_type(i_function),"_",fn_species(i_function),"_", &
             fn_n(i_function),"_",l_shell_str,".dat"
        open(50, file=output_name)
        write(50,*) "# ",n_grid(fn_species(i_function))
        do i_grid = 1, n_grid(fn_species(i_function)), 1
           write(50,*) &
                r_grid(i_grid, fn_species(i_function)), &
                basis_kinetic(i_grid, i_function)
        enddo
        close(50)
     enddo
  end if

  ! Well, n_basis is already set.  But resetting doesn't hurt.
  call get_bas_dimensions(n_ext_fns, &
  &                       fn_l(1:n_ext_fns), fn_species(1:n_ext_fns), &
  &                       max_ext_L, max_n_ext_fnLsp, n_ext)


  ! VB: BEFORE this point, anything to do with basis function manipulation is handled.
  !     AFTER this point, only splining, and manipulation and reorganizing of splines.

  !     can now allocate the actual basis storage arrays from module basis .
  call allocate_ext()
  !     and can store the splined versions of each basis function.
  !
  !     Note. We also "doctor" the splines, so that they are strictly zero outside
  !     a given radius, and so that they provide a smooth non-oscillating transition to zero.
  !     This is not identical with outer_radius further below, which is detemined by a threshold
  !     parameter, rather than a "strictly zero" criterion.
  !
  do i_function = 1, n_ext_fns, 1
     ! determine first "strictly zero" point on logarithmic grid
     outer_point = n_grid(fn_species(i_function))
     do while ( ( outer_point.ge.1 ) &
          .and.(basis_wave(outer_point,i_function).eq.0.d0) )
        outer_point = outer_point - 1
     enddo
     if (outer_point.lt.n_grid(fn_species(i_function))) then
        outer_point = outer_point+1
     end if
     ! create radial function spline
     ext_wave_spl(1:n_max_spline,1:n_max_grid,i_function) = 0.d0
     call cubic_spline &
          ( basis_wave(1,i_function), outer_point, &
          ext_wave_spl(1,1,i_function) )
     ! fudge outermost segment of spline so that it goes to zero smoothly as
     ! a parabola (but with a discontinuous forst derivative at the first non-zero
     ! logarithmic grid point)
     ext_wave_spl(2,outer_point-1,i_function) = &
          -2.d0 * ext_wave_spl(1,outer_point-1,i_function)
     ext_wave_spl(3,outer_point-1,i_function) = &
          ext_wave_spl(1,outer_point-1,i_function)
     ext_wave_spl(4,outer_point-1,i_function) = 0.d0
     ! create radial kinetic energy part spline
     ext_kinetic_spl(1:n_max_spline,1:n_max_grid,i_function) = 0.d0
     call cubic_spline &
          ( basis_kinetic(1,i_function), outer_point, &
          ext_kinetic_spl(1,1,i_function) )
     ! fudge outermost segment of spline so that it goes to zero smoothly as
     ! a parabola (but with a discontinuous forst derivative at the first non-zero
     ! logarithmic grid point)
     ext_kinetic_spl(2,outer_point-1,i_function) = &
          -2.d0 * ext_kinetic_spl(1,outer_point-1,i_function)
     ext_kinetic_spl(3,outer_point-1,i_function) = &
          ext_kinetic_spl(1,outer_point-1,i_function)
     ext_kinetic_spl(4,outer_point-1,i_function) = 0.d0
     if (use_basis_gradients) then
        ! create radial derivative spline
        ext_deriv_spl(1:n_max_spline,1:n_max_grid,i_function) = 0.d0
        call cubic_spline &
             ( basis_deriv(1,i_function), outer_point, &
             ext_deriv_spl(1,1,i_function) )
        ! fudge outermost segment of spline so that it goes to zero smoothly as
        ! a parabola (but with a discontinuous forst derivative at the first non-zero
        ! logarithmic grid point)
        ext_deriv_spl(2,outer_point-1,i_function) = &
             -2.d0 * ext_deriv_spl(1,outer_point-1,i_function)
        ext_deriv_spl(3,outer_point-1,i_function) = &
             ext_deriv_spl(1,outer_point-1,i_function)
        ext_deriv_spl(4,outer_point-1,i_function) = 0.d0
     end if
     
     if (flag_rel .eq. REL_atomic_zora) then
        ! create spline of kinetic energy part needed in scaled zora
        ext_kinetic_scaled_zora_spl &
             (1:n_max_spline,1:n_max_grid,i_function) = 0.d0
        call cubic_spline &
             ( basis_kinetic_scaled_zora(1,i_function), &
             outer_point, &
             ext_kinetic_scaled_zora_spl(1,1,i_function) )
        ! fudge outermost segment of spline so that it goes to zero smoothly as
        ! a parabola (but with a discontinuous forst derivative at the first non-zero
        ! logarithmic grid point)
        ext_kinetic_scaled_zora_spl(2,outer_point-1,i_function) = &
             -2.d0 * ext_kinetic_scaled_zora_spl &
             (1,outer_point-1,i_function)
        ext_kinetic_scaled_zora_spl(3,outer_point-1,i_function) = &
             ext_kinetic_scaled_zora_spl(1,outer_point-1,i_function)
        ext_kinetic_scaled_zora_spl(4,outer_point-1,i_function)=0.d0
     end if
  enddo ! end loop over basis fn splines
  
  ! If requested, do a separate verification of the basis function behavior near the cutoff.
  if (force_smooth_cutoff) then
     do i_function = 1, n_ext_fns, 1
        ! determine first "strictly zero" point on logarithmic grid
        outer_point = n_grid(fn_species(i_function))
        do while ( ( outer_point.ge.1 ) &
             .and.(basis_wave(outer_point,i_function).eq.0.d0) )
           outer_point = outer_point - 1
        enddo
        if (outer_point.lt.n_grid(fn_species(i_function))) then
           outer_point = outer_point+1
        end if
        ! radial function first
        if ( ( basis_wave(outer_point-2,i_function) &
             .gt.smooth_cutoff_limit) .or. &
             ( basis_wave(outer_point-1,i_function) &
             .gt.smooth_cutoff_limit) ) then
           write(use_unit,'(1X,A)') &
                "* After radial function orthonormalization / packing:"
           write(use_unit,'(1X,A,A,I5,A)') &
                "* Warning: You requested strict checking of the radial ", &
                "function cutoff, and function # ", i_function, &
                " exceeds the threshold at its two outermost points."
           stop
        end if
        ! kinetic function next
        if ( ( basis_kinetic(outer_point-2,i_function) &
             .gt.smooth_cutoff_limit) .or. &
             ( basis_kinetic(outer_point-1,i_function) &
             .gt.smooth_cutoff_limit) ) then
           
           write(use_unit,'(1X,A)') &
                "* After radial function orthonormalization / packing:"
           write(use_unit,'(1X,A,A,I5,A)') &
                "* Warning: You requested strict checking of the radial ", &
                "function cutoff, and kinetic # ", i_function, &
                " exceeds the threshold at its two outermost points."
           stop
        end if
        ! radial derivative next
        if (use_basis_gradients) then
           if ( ( basis_deriv(outer_point-2,i_function) &
                .gt.smooth_cutoff_limit) .or. &
                ( basis_deriv(outer_point-1,i_function) &
                .gt.smooth_cutoff_limit) ) then
              write(use_unit,'(1X,A)') &
                   "* After radial function orthonormalization / packing:"
              write(use_unit,'(1X,A,A,I5,A)') &
                   "* Warning: You requested strict checking of the radial ", &
                   "function cutoff, and deriv # ", i_function, &
                   " exceeds the threshold at its two outermost points."
              stop
           end if
        end if
     enddo
  end if
  ! end cutoff verification.
  
  ! get the maximium l value for aulixiary basis set
  l_ext_max = maxval(fn_l(1:n_ext_fns))

  ! Generate basis function indexing
  extfn_l = fn_l(1:n_ext_fns)
  extfn_species = fn_species(1:n_ext_fns)
  extfn_type = fn_type(1:n_ext_fns)
  extfn_n = fn_n(1:n_ext_fns)
  call generate_full_bas(n_ext_fns, extfn_l, extfn_species, &
  &                      max_ext_L, max_n_ext_fnLsp, n_ext, &
  &                      ext_atom, ext_l, ext_m, ext_fn, &
  &                      Lsp2n_ext_fnLsp, Lsp2ext_fn, Lsp2ext_sp, &
  &                      atom2ext_off, sp2n_ext_sp)
  max_n_ext_sp = maxval(sp2n_ext_sp)

  
!! ATOM BSSE:
!!  create the basis_mapping here:
!  if (calculate_atom_bsse) then

!  if (current_atom_num_bsse==1) then
!      write(info_str,'(2X,A)') "Creating the basis mapping for fresh BSSE geometry"
!      call localorb_info(info_str,use_unit,'(A)',OL_norm)
!      if (myid.eq.0) then
!        write(use_unit,'(2X,A,I8)') &
!           "| Atom number for atom-based counterpoise correction: ", current_atom_num_bsse
!      endif
!      do i_basis = 1, n_basis, 1
!        basis_mapping(i_basis) = i_basis
!      enddo
!  endif !loop for current_atom_num_bsse<1
  
!  if (current_atom_num_bsse>1) then
!      write(info_str,'(2X,A)') "Creating the basis mapping for fresh BSSE geometry"
!      call localorb_info(info_str,use_unit,'(A)',OL_norm)
!      if (myid.eq.0) then
!        write(use_unit,'(2X,A,I8)') &
!           "| Atom number for atom-based counterpoise correction: ", current_atom_num_bsse
!      endif

!      j_basis=0
!      basis_counter=0
!      empty_basis=0
      
!  ! for the 1st atoms previous to the _atom_num_bsse : just update the basis counter
!   do i_atom = 2, current_atom_num_bsse, 1
!     do i_function = species_first_function(species(i_atom)), species_last_function(species(i_atom))

!        i_l = fn_l(i_function)
!        do i_m = -i_l, i_l
!           empty_basis             = empty_basis+1
!        end do
!     end do
!    enddo      
!      do i_function = species_first_function(species(1)), &
!                     species_last_function(species(1))

!        i_l = fn_l(i_function)
!        do i_m = -i_l, i_l
!           empty_basis             = empty_basis+1
!           basis_counter= basis_counter+1
!           basis_mapping(basis_counter) = empty_basis
!        end do
!      end do  
      
!   do i_atom = 2, current_atom_num_bsse, 1
!     do i_function = species_first_function(species(i_atom)), species_last_function(species(i_atom))

!        i_l = fn_l(i_function)
!        do i_m = -i_l, i_l
!           j_basis             = j_basis+1
!           basis_counter= basis_counter+1
!           basis_mapping(basis_counter) = j_basis
!         end do
!     end do
!    enddo
    
!   do i_atom = current_atom_num_bsse+1,n_atoms, 1
!     do i_function = species_first_function(species(i_atom)), species_last_function(species(i_atom))

!        i_l = fn_l(i_function)
!        do i_m = -i_l, i_l
!           empty_basis             = empty_basis+1
!           basis_counter= basis_counter+1
!           basis_mapping(basis_counter) = empty_basis
!        end do
!     end do
!    enddo
    
!  endif ! for current_atom_num_bsse>1    

!  endif !  end stuff for atom_bsse
  
  !  All basis functions are stored and indexed. Now verify consistency.
  
  !     determine the outer radius of each basis function u(r)
  !     [i.e. the radius outside of which all integrations may be skipped
  !      because u(r) is practically zero]
  !     MUST CHECK BOTH ACTUAL FN AND SECOND DERIVATIVE
  do i_function = 1, n_ext_fns, 1
     i_grid = n_grid(fn_species(i_function))
     do while ( (abs(basis_wave(i_grid,i_function)).le.wave_threshold) .and. &
          (abs(basis_kinetic(i_grid,i_function)).le.wave_threshold) &
          .and.(i_grid.gt.1) )
        i_grid = i_grid-1
     enddo
     if (i_grid.le.1) then
        
        if (myid.eq.0) then
           write(use_unit,'(1X,A,A)') &
                "* Warning - a basis function is ",&
                "lower that the requested ", &
                "threshold value for integrations everywhere."
           write(use_unit,'(1X,A,A)') &
                "Species : ", trim(species_name(fn_species(i_function)))
           write(use_unit,'(1X,A,A)') &
                "Type    : ", fn_type(i_function)
           write(use_unit,'(1X,A,I3,A,I3)') &
                "(n,l)   : ", fn_n(i_function), ",", fn_l(i_function)
        end if
        
        stop
     end if
     if (i_grid.eq.n_grid(fn_species(i_function))) then
        outer_ext_radius(i_function) = &
             r_grid(i_grid,fn_species(i_function))
     else
        outer_ext_radius(i_function) = &
             r_grid(i_grid,fn_species(i_function))
     end if
  enddo
  
  if (use_basis_gradients) then
     !       also check first derivative
     do i_function = 1, n_ext_fns, 1
        i_grid = n_grid(fn_species(i_function))
        do while ( (abs(basis_deriv(i_grid,i_function)).le.wave_threshold) &
             .and. (i_grid.gt.1) )
           i_grid = i_grid-1
        enddo
        if (i_grid.eq.n_grid(fn_species(i_function))) then
           if ( r_grid(i_grid,fn_species(i_function)) &
                .gt.outer_ext_radius(i_function) ) then
              outer_ext_radius(i_function) = &
                   r_grid(i_grid,fn_species(i_function))
           end if
        else
           if ( r_grid(i_grid+1,fn_species(i_function)) &
                .gt.outer_ext_radius(i_function) ) then
              outer_ext_radius(i_function) = &
                   r_grid(i_grid,fn_species(i_function))
           end if
        end if
     enddo
  end if
  
  !------------ end order of basis functions----------
  
  atom_ext_radius_sq = 0.d0  
  do i_function = 1, n_ext_fns, 1
     outer_ext_radius_sq(i_function) = outer_ext_radius(i_function)**2
     atom_ext_radius_sq(fn_species(i_function)) = &
          max (atom_ext_radius_sq(fn_species(i_function)), &
          outer_ext_radius_sq(i_function) )
  enddo

  ! Make sure that multipole_radius_free >= atom_radius := max(outer_radius)
  !multipole_radius_free_sq = max(multipole_radius_free_sq, atom_ext_radius_sq)
  !multipole_radius_free = sqrt(multipole_radius_free_sq)
  atom_ext_radius = sqrt(atom_ext_radius_sq)

  ! Construct arrays to reconstruct fn<->species relations without
  ! basisfn_species.
  ext_fn_atom = .false.
  do i_basis_2 = 1, n_ext, 1
     ext_fn_atom(ext_fn(i_basis_2),ext_atom(i_basis_2)) = &
          .true.
  enddo
  i_offset_spl = 0
  species_offset = 0
  do i_species = 1, n_species, 1
     ext_fn_start_spl(i_species) = i_offset_spl + 1
     i_function2 = 0
     i_species_index = 0
     do i_basis_2 = 1, n_ext_fns, 1
        if (fn_species(i_basis_2).eq.i_species) then
           i_function2 = i_function2 + 1
           i_species_index(i_function2) = i_basis_2
        end if
     enddo
     n_function = i_function2
     n_ext_fn_species(i_species) = n_function
     do i_basis_2 = 1, n_function, 1
        temp_radius(i_basis_2) = &
             outer_ext_radius_sq(i_species_index(i_basis_2))
     enddo
     
     call insertionsort( temp_radius, n_function, perm_ext_fns, &
          perm_ext_fns_inv )
     ! Index array linking the actual order of the spline array to
     ! the original order of radial functions
     ! for Hamiltonian evaluation
     do i_spline = 1, n_function, 1
        i_offset = i_spline + species_offset
        ! store index of spline function as a function of the radial fn index
        perm_ext_fns_spl(i_offset) = &
             species_offset + &
             perm_ext_fns_inv(i_spline)
             i_ext_radial_fn(i_offset) = 0
        ! store the inverse, i.e. the index of the radial function as a
        ! function of the index used in the array of splined radial functions
        i_ext_radial_fn(i_offset) = species_offset + &
             perm_ext_fns(i_spline)
     enddo
  
     ! store species offset and increment it for the next species
     ext_spline_offset(i_species) = species_offset
     species_offset = species_offset + n_function
  
     
     ! and store all basis function splines in arrays in order of
     ! increasing outer_radius per species
     ! in principle, this loop also runs over _spline_ functions, not over
     ! radial functions in their original order ...
     do i_basis_2 = 1, n_function, 1
        i_offset_spl = i_offset_spl + 1
        ext_wave_ordered(i_offset_spl,1:4,1:n_grid(i_species)) = &
             ext_wave_spl(1:4,1:n_grid(i_species), &
             i_species_index(perm_ext_fns_inv(i_basis_2)))
        ext_kinetic_ordered(i_offset_spl,1:4,1:n_grid(i_species)) &
             = ext_kinetic_spl(1:4,1:n_grid(i_species), &
             i_species_index(perm_ext_fns_inv(i_basis_2)))
             
        if(flag_rel==REL_atomic_zora)then
           ext_kinetic_scaled_zora_ordered &
                (i_offset_spl,1:4,1:n_grid(i_species)) &
                = ext_kinetic_scaled_zora_spl &
                (1:4,1:n_grid(i_species), &
                i_species_index(perm_ext_fns_inv(i_basis_2)))
        end if
        if (use_basis_gradients) then
           ext_deriv_ordered(i_offset_spl,1:4,1:n_grid(i_species)) &
                = ext_deriv_spl(1:4,1:n_grid(i_species), &
                i_species_index(perm_ext_fns_inv(i_basis_2)))
        end if
     enddo
  enddo
    
  if (myid.eq.0) then
     write(use_unit,*)
     write(use_unit,'(2X,A)') &
          "Extended basis size parameters after reduction:"
     write (use_unit,'(2X,A,I8)') &
          "| Total number of radial functions: ", &
          n_ext_fns
     write (use_unit,'(2X,A,I8)') &
          "| Total number of basis functions : ", &
          n_ext
     write(use_unit,*)
  end if
  if ( allocated(basis_wave) ) then
     deallocate ( basis_wave )
  end if
  if ( allocated(basis_deriv) ) then
     deallocate ( basis_deriv )
  end if
  if ( allocated(basis_kinetic) ) then
     deallocate ( basis_kinetic )
  end if
  if(allocated ( basis_kinetic_scaled_zora))then
     deallocate ( basis_kinetic_scaled_zora)
  end if
  return
end subroutine shrink_fixed_large_basis
!******
!
!****s* FHI-aims/copy_fixed_basis_to_ext
!  NAME
!   copy_fixed_basis_to_ext
!  SYNOPSIS
subroutine copy_fixed_basis_to_ext_arrays ( )
!  PURPOSE
!   If there is no auxialiary basisis function present 
!   we copy data from normal basis to directly to auxiliary basis
!    
!   This could be done also using pointers as then we would not need to
!   copy anything so we could save some memory.
!
!  USES
    use dimensions,      only : l_ext_max, l_wave_max, n_atoms, n_basis, &
                                n_basis_fns, n_ext, n_ext_fns, n_max_grid, &
                                n_max_spline, n_species, use_basis_gradients
    use runtime_choices, only : flag_rel, REL_atomic_zora
    use basis,           only : basisfn_l, basis_wave_spl, ext_deriv_spl, &
                                basis_kinetic_scaled_zora_spl, basis_wave_ordered, &
                                basis_deriv_ordered, basis_deriv_ordered, &
                                basis_kinetic_scaled_zora_ordered, outer_radius, &
                                basis_atom, basisfn_species,  basis_kinetic_ordered, &
                                basis_l, basis_m, basis_fn, basisfn_species, &
                                basisfn_species, basisfn_n, Lsp2n_basis_fnLsp, &
                                Lsp2basis_fn, Lsp2basis_fn, Lsp2basis_sp, &
                                atom2basis_off, atom2basis_off, outer_radius_sq, &
                                atom_radius, atom_radius, perm_basis_fns_spl, &
                                basisfn_type, sp2n_basis_sp, atom_radius_sq, &
                                i_radial_fn, basis_fn_start_spl, n_basis_fn_species, &
                                basis_fn_atom, n_basis_atom, spline_offset, &
                                basis_mapping, atom2basis_off, atom_radius, &
                                atom_radius_sq, atom2ext_off, atom_ext_radius, &
                                atom_ext_radius_sq, ext_atom, ext_deriv_ordered, & 
                                ext_fn, ext_fn_atom, ext_fn_start_spl, sp2n_ext_sp, &
                                ext_kinetic_ordered, ext_kinetic_scaled_zora_ordered, &
                                ext_kinetic_scaled_zora_spl, ext_kinetic_spl, &
                                ext_l, ext_m, ext_mapping, ext_spline_offset, &
                                ext_wave_ordered, ext_wave_spl, extfn_l, extfn_n, &
                                extfn_species, extfn_type, i_ext_radial_fn, Lsp2ext_fn, &
                                Lsp2ext_sp, Lsp2n_ext_fnLsp, max_ext_L, max_n_basis_sp, &
                                max_n_ext_fnlsp, max_n_ext_sp, n_ext_atom, n_ext_fn_species, &
                                outer_ext_radius, outer_ext_radius_sq, perm_ext_fns_spl, &
                                allocate_ext
    implicit none


    l_ext_max = l_wave_max
    n_ext = n_basis
    n_ext_fns = n_basis_fns
    max_n_ext_sp = max_n_basis_sp

    call get_bas_dimensions(n_ext_fns, &
    &                       basisfn_l(1:n_basis_fns), basisfn_species(1:n_basis_fns), &
    &                       max_ext_L, max_n_ext_fnLsp, n_ext)

    call allocate_ext()

    ext_wave_spl = basis_wave_spl(1:n_max_spline,1:n_max_grid, 1:n_ext_fns)
    if (use_basis_gradients) then
        ext_deriv_spl = ext_deriv_spl(1:n_max_spline,1:n_max_grid,1:n_ext_fns)
    end if
    ext_kinetic_spl = basis_wave_spl(1:n_max_spline,1:n_max_grid,1:n_ext_fns)
    if (flag_rel==REL_atomic_zora) then
        ext_kinetic_scaled_zora_spl =  & 
            basis_kinetic_scaled_zora_spl(1:n_max_spline,1:n_max_grid,1:n_ext_fns)
    end if      
    ext_wave_ordered = basis_wave_ordered(1:n_ext_fns,1:n_max_spline,1:n_max_grid)
    if (use_basis_gradients) then
        ext_deriv_ordered = & 
            basis_deriv_ordered(1:n_ext_fns,1:n_max_spline,1:n_max_grid )
    end if
    ext_kinetic_ordered = basis_kinetic_ordered(1:n_ext_fns,1:n_max_spline, 1:n_max_grid)
    if (flag_rel==REL_atomic_zora) then
        ext_kinetic_scaled_zora_ordered = &
            basis_kinetic_scaled_zora_ordered(1:n_ext_fns,1:n_max_spline,1:n_max_grid)
    end if
    outer_ext_radius = outer_radius(1:n_ext_fns)
    ext_atom = basis_atom(1:n_ext)
    ext_l = basis_l(1:n_ext)
    ext_m = basis_m (1:n_ext)
    ext_fn = basis_fn(1:n_ext)
    extfn_species = basisfn_species(1:n_ext_fns)
    extfn_l = basisfn_l(1:n_ext_fns)
    extfn_type = basisfn_type(1:n_ext_fns)
    extfn_n = basisfn_n(1:n_ext_fns)
    Lsp2n_ext_fnLsp = Lsp2n_basis_fnLsp(0:max_ext_L, 1:n_species)
    Lsp2ext_fn = Lsp2basis_fn(1:max_n_ext_fnLsp, 0:max_ext_L,1:n_species)
    Lsp2ext_sp = Lsp2basis_sp(1:max_n_ext_fnLsp, 0:max_ext_L, 1:n_species)
    atom2ext_off = atom2basis_off(1:n_atoms)
    sp2n_ext_sp = sp2n_basis_sp(1:n_species)
    outer_ext_radius_sq = outer_radius_sq(1:n_ext_fns)
    atom_ext_radius = atom_radius(1:n_species)
    atom_ext_radius_sq = atom_radius_sq(1:n_species)
    perm_ext_fns_spl = perm_basis_fns_spl(1:n_ext_fns)
    i_ext_radial_fn = i_radial_fn(1:n_ext_fns)
    ext_fn_start_spl = basis_fn_start_spl(1:n_species)
    n_ext_fn_species = n_basis_fn_species(1:n_species)
    ext_fn_atom = basis_fn_atom(1:n_ext_fns,1:n_atoms)
    n_ext_atom = n_basis_atom(1:n_atoms)
    ext_spline_offset = spline_offset(1:n_species)
    ext_mapping = basis_mapping(1:n_ext)
    
 
end subroutine copy_fixed_basis_to_ext_arrays    
