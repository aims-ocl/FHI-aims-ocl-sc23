!****h* FHI-aims/precondition
!  NAME
!    precondition 
!  SYNOPSIS

module precondition
!  PURPOSE
!    Charge density preconditioner according to a Kerker-type scheme
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  SOURCE

  logical :: kerker_preconditioner_on
  logical :: preconditioner_first_warnings = .true.
  integer, private :: n_spline_grid_dim, prec_l_dim

  real*8,  dimension(:,:,:)   , allocatable, private :: green_I, green_K  ! radial Green functions of modified Helmholtz Eqn, fct of r, l, and species
  real*8,  dimension(:,:,:)   , allocatable, private :: R_multipole       ! angular and radially resolved multipole components
  real*8,  dimension(:,:,:,:) , allocatable, private :: R_multipole_spl   ! the spline functions of above multipole components - FIXME: use same mem as hartree pot!!
  real*8,  dimension(:,:)     , allocatable, private :: R_prec_multipole  ! the multipole version of the output residual
  real*8,  dimension(:,:)     , allocatable, private :: aux_R_spline            ! temporarily contains the splines for one angular momentum shell on one atom: rad grid
  real*8,  dimension(:,:)     , allocatable, private :: aux_R_prec_spline       ! same for logarithmic grid
  real*8,  dimension(:,:)     , allocatable, private :: angular_integral_log    ! temporarily contains the values of R_spl(alpha,l,m,r)
  real*8,  dimension(:)       , allocatable, private :: integral_zero_r         ! running variables for Green function integration
  real*8,  dimension(:)       , allocatable, private :: integral_r_infty        
  real*8,  dimension(:,:,:)   , allocatable, private :: current_R_multipole_spl ! to store the multipole splines of one given atom for calculation of final answer
  real*8,  dimension(:)       , allocatable, private :: max_lm_spl_rad_sq       ! stores the maximum relevant radius of the logarithmically splined R_prec(l,m)
  real*8,  dimension(:)       , allocatable, private :: aux_R_prec_result       ! stores the spline result for R_prec in the final calculation
  real*8,  dimension(:,:)     , allocatable, private :: inner_radial_spl        ! stores the splines inside the two innermost radial points in real space
  integer, dimension(:,:)     , allocatable, private :: index_lm                ! (m,l)->1D translation array
  integer, dimension(:)       , allocatable, private :: prec_max_l_species     
  integer, dimension(:)       , allocatable, private :: prec_l_dim_species
 

!******

 contains
!------------------------------------------------------------------------------
!****s* preconditioner/prepare_preconditioner
!  NAME
!    prepare_preconditioner
!  SYNOPSIS
subroutine prepare_preconditioner
!  PURPOSE
!    This takes care of initialize all global variables for preconditioning.
!  USES
  use runtime_choices
  use dimensions
  use grids
  use localorb_io
  use species_data

!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SOURCE
  implicit none
  integer :: i_species, i_grid, i_m, i_l, i_index
  character*100 :: info_str


  ! ---------- Allocation and initialization ----------------------------------
  call deallocate_preconditioner    ! Just in case the spline distribution & other things have failed ... 

  ! go through the various species to check on the size of the angular momentum expansion
  if (.not.allocated(prec_max_l_species))      allocate(prec_max_l_species(n_species))
  if (.not.allocated(prec_l_dim_species))      allocate(prec_l_dim_species(n_species))

  prec_max_l_species(:) = precondition_max_l    ! initialize species-dependent l-values
  i_l = 0
  do i_species = 1, n_species
     if (l_hartree(i_species).lt.precondition_max_l) then
        write(info_str,'(2X,2A)') '* WARNING!!!! l_hartree < precondition_max_l for species ', species_name(i_species)
        call localorb_info(info_str)
        write(info_str,'(2X,A)')  '* The preconditioner would crash under this condition.'
        call localorb_info(info_str)
        write(info_str,'(2X,A)')  '* Setting precondition_max_l = l_hartree for this species to prevent difficulties.'
        call localorb_info(info_str)
        prec_max_l_species(i_species) = l_hartree(i_species)
     end if
     i_l = max(i_l,prec_max_l_species(i_species))
  end do
  if (i_l.lt.precondition_max_l) then
     write(info_str,'(2X,A)') '* The above warnings cause AIMS to reduce precondition_max_l to an overall'
     call localorb_info(info_str)
     write(info_str,'(2X,A,I4,A)') '* value of ',i_l,' for overall consisitency.'
     call localorb_info(info_str)
     precondition_max_l = i_l
  end if

  prec_l_dim_species(:) = (prec_max_l_species(:)+1)**2

  ! R_multipole is only used in this routine, but it does contain a runtime choice of the max number of angular momentum shells
  ! same goes for R_multipole_spl
  prec_l_dim = (precondition_max_l+1)**2                        ! total number of (l,m) pairs - this limit is easier to handle
  if (.not.allocated(aux_R_spline))       allocate(aux_R_spline(n_max_spline, n_max_radial+1))
  if (.not.allocated(aux_R_prec_spline))  allocate(aux_R_prec_spline(n_max_spline, n_max_grid))
  if (.not.allocated(R_multipole))        allocate(R_multipole(n_max_radial+1, prec_l_dim, n_atoms))   
  ! this array is used twice in a row with different meanings each time, but one meaning after the other
  ! hence, over-allocation is not exactly necessary
  n_spline_grid_dim = max(n_max_grid,n_max_radial+1)

  if (.not.allocated(R_multipole_spl)) allocate(R_multipole_spl(prec_l_dim, n_max_spline, n_spline_grid_dim, n_spline_atoms))
  if (.not.allocated(green_I))                 allocate(green_I(n_max_grid, precondition_max_l+1, n_species))
  if (.not.allocated(green_K))                 allocate(green_K(n_max_grid, precondition_max_l+1, n_species))
  if (.not.allocated(angular_integral_log))    allocate(angular_integral_log(prec_l_dim, n_max_grid ))
  if (.not.allocated(integral_zero_r))         allocate(integral_zero_r(prec_l_dim))
  if (.not.allocated(integral_r_infty))        allocate(integral_r_infty(prec_l_dim))
  if (.not.allocated(R_prec_multipole))        allocate(R_prec_multipole(n_max_grid, prec_l_dim))
  if (.not.allocated(current_R_multipole_spl)) allocate(current_R_multipole_spl(prec_l_dim, n_max_spline, n_spline_grid_dim))
  if (.not.allocated(max_lm_spl_rad_sq))       allocate(max_lm_spl_rad_sq(prec_l_dim))
  if (.not.allocated(aux_R_prec_result))       allocate(aux_R_prec_result(prec_l_dim))
  if (.not.allocated(index_lm))                allocate(index_lm(-precondition_max_l:precondition_max_l,0:precondition_max_l))
  if (.not.allocated(inner_radial_spl))        allocate(inner_radial_spl(3,prec_l_dim))

  ! map each (m,l) to an index and keep in storage for later use when splining
  i_index = 0
  index_lm = 0
  do i_l = 0, precondition_max_l, 1
     do i_m = -i_l, i_l
        i_index = i_index + 1
        index_lm(i_m, i_l) = i_index
     enddo
  enddo

  ! initialize the radial Green functions for this particular value of q0; might possibly change throughout calculation ??????
  do i_species = 1, n_species
     do i_grid = 1, n_grid(i_species)
        
        ! calculate all l-values at once - for the radial argument (q_0 r)
        ! FIXME: investigate whether or not is is possible to switch arguments 1 and 2 in green_I and green_K for speed reasons?
        call bessel_I(precondition_kerker_q0*r_grid(i_grid,i_species),precondition_max_l, green_I(i_grid, :, i_species))
        call bessel_K(precondition_kerker_q0*r_grid(i_grid,i_species),precondition_max_l, green_K(i_grid, :, i_species))
        green_I(i_grid,:,i_species) = green_I(i_grid,:,i_species)/sqrt(precondition_kerker_q0*r_grid(i_grid,i_species))
        green_K(i_grid,:,i_species) = green_K(i_grid,:,i_species)/sqrt(precondition_kerker_q0*r_grid(i_grid,i_species))
     end do         ! radial loop, initialization of bessel green function  
  end do            ! species loop, initialization of bessel fct Green function
end subroutine prepare_preconditioner
!******

!------------------------------------------------------------------------------
!****s* preconditioner/deallocate_preconditioner
!  NAME
!    deallocate_preconditioner
!  SYNOPSIS
subroutine deallocate_preconditioner 
!  PURPOSE
!    Clean up all of the variables allocated above
!  USES
  use runtime_choices
  use dimensions
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SOURCE
  implicit none 
  if (allocated(inner_radial_spl))        deallocate(inner_radial_spl)
  if (allocated(prec_max_l_species))      deallocate(prec_max_l_species)
  if (allocated(prec_l_dim_species))      deallocate(prec_l_dim_species)
  if (allocated(aux_R_prec_result))       deallocate(aux_R_prec_result)
  if (allocated(max_lm_spl_rad_sq))       deallocate(max_lm_spl_rad_sq)
  if (allocated(current_R_multipole_spl)) deallocate(current_R_multipole_spl)
  if (allocated(R_prec_multipole))        deallocate(R_prec_multipole)
  if (allocated(integral_r_infty))        deallocate(integral_r_infty)
  if (allocated(integral_zero_r))         deallocate(integral_zero_r)
  if (allocated(angular_integral_log))    deallocate(angular_integral_log)
  if (allocated(green_I))                 deallocate(green_I)
  if (allocated(green_K))                 deallocate(green_K)
  if (allocated(R_multipole_spl))         deallocate(R_multipole_spl)
  if (allocated(R_multipole))             deallocate(R_multipole)
  if (allocated(aux_R_spline))            deallocate(aux_R_spline)
  if (allocated(aux_R_prec_spline))       deallocate(aux_R_prec_spline)
  if (allocated(index_lm))                deallocate(index_lm)
!******
end subroutine deallocate_preconditioner


!------------------------------------------------------------------------------
!****s* preconditioner/precondition_kerker
!  NAME
!    precondition_kerker
!  SYNOPSIS
subroutine precondition_kerker(R_in, hartree_partition_tab)    
!  PURPOSE
!    actual preconditioning routine
!
! Routine that takes an input residual R_in and preconditions it to yield some output 
! R_prec, in fourier space this is described as 
! 
!       R_prec(q) = q^2 R_in(q)/(q^2 + q_zero^2) 
!                 = (1 + q_zero^2/(q^2+q_zero^2))R_in(q)
!                 = R_in(q) + q_zero^2 R_corr(q)
!
! The correction is the solution to the PDE
!
!      (\nabla^2-q_zero^2) R_corr = -R_in
!
! which is solved via a multipole expansion & green function integration
!
! The actual residuals are of course spin resolved, but that problem is handled in the calling sequence of this
! routine, this part of the code does not need to know anything about spin at all. 
!  USES
  use grids           ! information on integration grids and such 
  use dimensions      ! here each thread learns about its local memory structure
  use constants       ! 4 pi required in the very last line ... 
  use localorb_io     ! I/O only on thread 0, i.e. good for parallel playing
  use runtime_choices ! all control.in choices
  use mpi_tasks       ! contains subroutine get_my_task() and the actual ID if that is ever needed in here. 
  use mpi_utilities   ! required for some spline storage information
  use spline          ! contains the actual spline routines
  use pbc_lists       ! for atom-atom distance tab center_to_atom - which should work for PBC as well ...
  use species_data
  use synchronize_mpi
  use geometry, only: species

  implicit none

!  ARGUMENTS

  real*8, dimension(n_full_points)      :: R_in                ! input residual from previous scf cycle
  real*8, dimension(n_full_points)      :: hartree_partition_tab       

!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!   o R_in (input/output)   - on input: The charge density residual over the entire grid
!   o hartree_partition_tab - partition tab for the integrations
!  OUTPUT
!   o R_out (input/output)  - on output: preconditioned density residual
!  SOURCE



  real*8, dimension(n_full_points)      :: R_prec              ! preconditioned output residual: to be calculated

  ! various counter and integer tracking variables, as well as real number buffer variables
  integer :: i_my_batch, i_index, current_atom, current_radial, current_angular
  integer :: i_full_points, i_atom, current_spl_atom, i_m, i_l, i_radial, i_spline
  integer :: i_species, mpi_err, prev_atom, current_center, i_grid, i_center, n_grid_limit
  integer :: n_inner_grid
  real*8  :: radius, alpha, i_r_outer, delta, max_spl_rad_sq, dir_tab(3), dist_tab_sq
  real*8  :: coord_current(3), trigonom_tab(4), delta_2, delta_3, ddot, multipole_thresh
  real*8  :: integral_R_prec, integral_volume, integral_charge, r1, r2, f1, f2, denominator
  real*8  :: radius_sq
  character*80 :: info_str
  real*8, dimension(prec_l_dim) :: ylm_tab
 
  ! FIXME: make explicit parameter if necessary - threshold to determine extend of multipole expansion
  multipole_thresh = 1d-12

  R_multipole      = 0d0                                         ! better safe than sorry ... 
  R_multipole_spl  = 0d0
  R_prec_multipole = 0d0

  ! multipole expansion is done according to Delley (1990), modifying Eqn (11) according to this particular purpose:
  ! R_multipole(s) = Int d\Omega ylm(r-r_atom) hartree_partition_tab R_in(r-r_atom)
  ! i.e.  R_in = partition_function * R_total already contains the partitioned residual for each atom,
  !       this routine will only have to play with the atomic contributions and add up all
  !       points in a given atomic integration shell later. 

  i_full_points = 0                                          ! this counts the points on a thread - needed as index for the partition tab
  do i_my_batch = 1, n_my_batches, 1                          ! loop through all total batches

        do i_index = 1, batches(i_my_batch)%size, 1          ! loop through all points in batch

           ! get the technical details of current grid point
           current_atom    = batches(i_my_batch) % points(i_index) % index_atom
           current_radial  = batches(i_my_batch) % points(i_index) % index_radial
           current_angular = batches(i_my_batch) % points(i_index) % index_angular
           i_full_points   = i_full_points + 1
 
           ! treat grid point in question to obtain multipole expansion 
           if (hartree_partition_tab(i_full_points).gt.0.d0) then ! - well, only if there is a non-zero summation weight of course 

              ! look up the sperical harmonics for this point .... 
              ylm_tab (1:prec_l_dim_species(species(current_atom))) = &
                   local_ylm_tab(1:prec_l_dim_species(species(current_atom)),current_angular, &      ! local_ylm_tab contains information for all lm values serially
                   lebedev_grid_index(current_radial,species(current_atom)))

              ! ... to calculate the multipole expansion contribution by implied looping over all (l,m)
              R_multipole(current_radial, 1:prec_l_dim_species(species(current_atom)), current_atom) = &
                   R_multipole(current_radial, 1:prec_l_dim_species(species(current_atom)), current_atom) + &          ! other points are in this integration shell too don't forget them
                   ylm_tab(1:prec_l_dim_species(species(current_atom)))*hartree_partition_tab(i_full_points)&
                   * R_in(i_full_points)

           end if       !  do we need to play with this grid point, i.e. (hartree_partition_tab(i_full_points).gt.0.d0) ?
        end do          ! index = 1, batches(i_mybatch)%size
     ! end if             ! (myid.eq.batch_task_list(i_batch))
  end do                ! i_batch = 1, n_grid_batches, 1
  
  ! synchronize multipole expansion across all threads
  call sync_kerker_multipole(R_multipole, n_max_radial+1, prec_l_dim)


  ! spline the above multipole expansion for use with the logarithmic grid
  !       code is the same as in update_hartree_potential_p1; except that two array arguments in the splines are switched as suggested in that routine
  !   
  ! this part of the code runs over all atoms - we are no longer on the original grid because there now is a multipole 
  ! expansion that has a basis which is independent of the grid points. 
  ! This also means that each thread must know whether or not the atom under consideration is in its task list ...

  ! The same loop also calculates the atomic multipole expansion of the preconditioned residual akin to routine integrate_hartree_log_grid_p1
    
  ! loop over all atoms
  do i_atom = 1, n_atoms

     ! set end of multipole expansion to zero explicitly, - to ensure that there is no diverging charge density
     ! outside the known grid quantities - this should be done on all atoms; this is done again on the finally splined multipole density later
     R_multipole(n_radial(species(i_atom))+1,:,i_atom) = 0.d0
     
     ! atomic task list distribution - why work with atomic tasks
     if (myid.eq.task_list(i_atom)) then

        ! spline_atom_storage(n_atoms) is defined on each thread and contains the order of atoms which are treated by this 
        ! particular thread - in case that distributed spline_storage is required, the (atomic) index for a spline 
        ! in the spline storage array is given by spline_atom_storage(i_atom)
        current_spl_atom = spline_atom_storage(i_atom)
        
        ! do the splining itself: separately, for each atom, and each angular momentum shell. 
        ! the treatment of all radial shells is done in the splining routine
        do i_l = 0, prec_max_l_species(species(i_atom)), 1
           do i_m = - i_l, i_l, 1

              call cubic_spline &
                   (R_multipole(:,index_lm(i_m, i_l),i_atom), &
                   n_radial(species(i_atom))+1, &
                   aux_R_spline )
              
              ! copy the array of spline coefficients to its eventual location ...
              do i_radial = 1, n_radial(species(i_atom)), 1
                 do i_spline = 1, n_max_spline, 1
                    R_multipole_spl( index_lm(i_m, i_l), i_spline, i_radial+1, current_spl_atom) &
                         = aux_R_spline(i_spline,i_radial)
                 enddo
              enddo
           end do
        end do              ! end splining loop
        
        ! doctor the atomic splines such that they don't yield infinite residuals in the far field - just in case ... 
        ! find outermost radial grid point that is possibly non-zero
        i_radial = n_radial(species(i_atom))
        do while ( ( r_radial(i_radial,species(i_atom)) .ge. &
             multipole_radius_free(species(i_atom)) ) &
             .and.(i_radial.gt.1) )
           R_multipole_spl(:,:,i_radial,current_spl_atom) = 0.d0                 ! set whole spline to zero for non-valid points
           i_radial = i_radial - 1                                               ! next innermost radial point
        enddo

        ! Outermost atom radius in units of the radial integration grid  - this is for the free atom ... 
        ! FH: I did not go through the following spline-doctoring procedure in detail, but copied it from update_hartree_potential_p1 - 
        !     trusting that it was tested there in detail. 
        i_r_outer = invert_radial_grid &
             ( multipole_radius_free(species(i_atom)), &
             n_radial(species(i_atom)), &
             scale_radial(species(i_atom)) ) + 1

        ! difference from outermost density point to outermost free atom point
        delta = i_r_outer - i_radial
        delta_2 = delta*delta
        delta_3 = delta_2*delta

        ! doctor the spline coefficients at the outermost finite value
        ! i_radial to go smoothly to zero at multipole_radius_free
        do i_l = 0, prec_max_l_species(species(i_atom)), 1
           do i_m = - i_l, i_l, 1
              R_multipole_spl( index_lm(i_m, i_l), 3, i_radial, current_spl_atom ) = &
                   - 3.d0 / delta_2 * &
                   R_multipole_spl( index_lm(i_m, i_l), 1, i_radial, current_spl_atom ) &
                   - 2.d0 / delta * &
                   R_multipole_spl( index_lm(i_m, i_l), 2, i_radial, current_spl_atom )
              R_multipole_spl( index_lm(i_m, i_l), 4, i_radial, current_spl_atom ) = &
                   2.d0 / delta_3 * &
                   R_multipole_spl( index_lm(i_m, i_l), 1, i_radial, current_spl_atom ) &
                   + 1.d0 / delta_2 * &
                   R_multipole_spl( index_lm(i_m, i_l), 2, i_radial, current_spl_atom )
           enddo
        enddo

        ! figure out parabolic splines for the inside of r_radial(2,species(i_atom))
        ! These are REAL-SPACE splines that are INDEPENDENT of any mapping of an index/radius onto a 
        ! nonlinear grid and are there to take out a divergence of the radial grid mapping near the 
        ! origin.
        
        ! all the splines are written as parabolas to make evaluation a lot simpler 
        inner_radial_spl = 0d0
        ! calculate the inner spline for l = 0, m = 0: linear, going through f(r1) and f(r2)
        r1 = r_radial(1,species(i_atom))
        r2 = r_radial(2,species(i_atom))
        f1 = R_multipole(1,index_lm(0,0),i_atom)
        f2 = R_multipole(2,index_lm(0,0),i_atom)
        inner_radial_spl(2,index_lm(0,0)) = (f2-f1)/(r2-r1)     ! linear term 
        inner_radial_spl(1,index_lm(0,0)) = f1-r1*inner_radial_spl(2, index_lm(0,0)) ! constant term

        ! calculate the inner spline for l > 0, m - parabolic, going through f(r1), f(r2), and (0,0)
        denominator = (r1*r1*r2-r2*r2*r1)
        do i_l = 1, prec_max_l_species(species(i_atom)), 1
           do i_m = -i_l, i_l
              f1 = R_multipole(1, index_lm(i_m, i_l), i_atom)
              f2 = R_multipole(2, index_lm(i_m, i_l), i_atom)
              ! constant term = 0
              inner_radial_spl(2, index_lm(i_m, i_l)) = (r1*r1*f2-r2*r2*f1)/denominator  ! linear term
              inner_radial_spl(3, index_lm(i_m, i_l)) = (r2*f1-f2*r1)/denominator        ! parabolic term
           end do
        end do


        ! integrate the Green function integral (12b) in Delley, 1990; with the modified Green functions I_{l+1/2}(q_0 r) and K_{l+1/2}(q_0 r)
        ! for each atom; for now on the logarithmic grid and using the splines that were just computed. 
        ! what follows is also done in the routine integrate_hartree_log_grid

        ! FIXME: why do we even have to store all the splines? In case that this is a memory bottle-neck, there might be a large scale 
        !        workaround by working on each (l,m) value separately, store the temporary splines and then continue to evaluate the 
        !        completed R_prec(alpha,l,m,r) ... which would then have to be synchronized again, of course 
        
        ! FIXME: how are the speed, memory requirements, and accuracy of the overall result affected by the calculation of the 
        !        splines vs the use of the regular radial grids for the integrations? The latter should be much faster but should also
        !        yield a lot less accurate results. Does that matter for a preconditioning correction to the density????
        
        angular_integral_log = 0d0
        alpha                = log(r_grid_inc(species(i_atom)))   ! prefactor for the calculation of the grid radius & the integration weight
        
        
        ! calculate the indexes for integration limits
        ! n_inner_grid = the point where the second radial shell starts, i.e. where we go to the regular splines. 
        n_inner_grid = invert_log_grid(r_radial(2,species(i_atom)), &
             r_grid_min(species(i_atom)),r_grid_inc(species(i_atom)))

        ! n_grid_limit = the point after which there is no density anyway, i.e. the free atom multipole radius
        n_grid_limit = invert_log_grid(multipole_radius_free(species(i_atom))+extra_adding_to_hartree_potential_distance, &
             r_grid_min(species(i_atom)),r_grid_inc(species(i_atom)))

        ! first spline evaluation loop: 
        ! everything inside the innermost two radial shells.
        do i_grid = 1, n_inner_grid
           
           ! this is the actual radius known to the parabolic splines !
           radius = r_grid(i_grid,species(i_atom))
           radius_sq = radius*radius

           ! evaluate real-space parabolic spline for all (l,m)
           angular_integral_log(:, i_grid) = &
                inner_radial_spl(1,:)+inner_radial_spl(2,:)*radius+inner_radial_spl(3,:)*radius_sq
           
           ! transform log grid onto uniform index space for integration
           angular_integral_log(:,i_grid) = alpha*(radius**3)*angular_integral_log(:,i_grid)
           
        end do

        ! second spline evaluation loop on logarithmic grid:
        ! everything for which we know to have well-defined splines. 
        do i_grid = n_inner_grid, n_grid_limit + 1
           
              
           ! calculate index of grid point on radial grid, where the splines reside
           radius = invert_radial_grid( r_grid(i_grid,species(i_atom)), n_radial(species(i_atom)), &
                scale_radial(species(i_atom)))+1d0
           
           ! evaluate all (l,m) splines for this particular grid point    - arguments in the declaration of spline_vector ... 
           call spline_vector( radius,                                  & !  real*8  :: r_output
                R_multipole_spl(:,:,:,current_spl_atom), & !  real*8  :: spl_param(n_l_dim,4,n_grid_dim) - for this particular atom only
                n_max_radial+1,                          & !  integer :: n_grid_dim  
                prec_l_dim,                              & !  integer :: n_l_dim
                n_radial(species(i_atom))+1,             & !  integer :: n_points
                prec_l_dim_species(species(i_atom)),     & !  integer :: n_vector -- the number of points actually evaluated
                angular_integral_log(:,i_grid))            !  real*8  :: out_result(n_vector)                           
           
           ! transform log grid onto uniform index space for integration
           angular_integral_log(:,i_grid) = alpha*(r_grid(i_grid,species(i_atom))**3)*angular_integral_log(:,i_grid)
           
        end do   ! end evaluation of splines across the whole grid, for atom current_spl_atom
        
        ! we now know, for i_atom, the values of the splined multipole residual R_in on the logarithmic grid 
        !                                            - let's calculate R_prec_multipole on the same logarithmic grid. 
        ! This is done using the Adams-Moulton linear multistep integrator, using up to 4 terms. 
        !    see http://en.wikipedia.org/wiki/Linear_multistep_method or Abramowitz/Stegun p. 896 for details. 

        ! as in integrate_hartree_log_grid, do so in two parts
        ! part (1) contains the integral from 0 up to the grid radius r, with I_{l+1/2} in the integrand and 4 pi K_{l+1/2} as the prefactor
        ! part (2) is the other way around, i.e. from r to "infinity"    with K_{l+1/2} in the integrand and 4 pi I_{l+1/2} as the prefactor

        ! first integral: from zero to radius r, for all possible r

        ! first three terms: do one-by-one with increasing order of integration method. 
        integral_zero_r  = 0d0
        do i_l = 0, prec_max_l_species(species(i_atom)), 1
           do i_m = - i_l, i_l, 1

             ! TERM 1 : Integral_1 = h*f_1; but h = 1
             integral_zero_r(index_lm(i_m, i_l)) =  & 
                   angular_integral_log(index_lm(i_m, i_l),1) * green_I(1,i_l+1,species(i_atom))
             R_prec_multipole(1,index_lm(i_m,i_l)) =  &
                   integral_zero_r(index_lm(i_m, i_l) ) * green_K(1,i_l+1,species(i_atom))

             ! TERM 2 : Integral_2 = Integral_1 + h(f_2+f_1)/2 
             integral_zero_r(index_lm(i_m, i_l)) =  integral_zero_r(index_lm(i_m, i_l)) & 
                 + ( 1d0*angular_integral_log(index_lm(i_m, i_l),2) * green_I(2,i_l+1,species(i_atom)) &
                   + 1d0*angular_integral_log(index_lm(i_m, i_l),1) * green_I(1,i_l+1,species(i_atom)))/2d0
             R_prec_multipole(2,index_lm(i_m,i_l)) =  &
                   integral_zero_r(index_lm(i_m, i_l) ) * green_K(2,i_l+1,species(i_atom))

             ! TERM 3 : Integral_3 = Integral_2 + h(5f_3 + 8f_2 - f_1)/12
             integral_zero_r(index_lm(i_m, i_l)) =  integral_zero_r(index_lm(i_m, i_l)) & 
                 + ( 5d0*angular_integral_log(index_lm(i_m, i_l),3) * green_I(3,i_l+1,species(i_atom)) &
                   + 8d0*angular_integral_log(index_lm(i_m, i_l),2) * green_I(2,i_l+1,species(i_atom)) &
                   - 1d0*angular_integral_log(index_lm(i_m, i_l),1) * green_I(1,i_l+1,species(i_atom)))/12d0
             R_prec_multipole(3,index_lm(i_m,i_l)) =  &
                   integral_zero_r(index_lm(i_m, i_l) ) * green_K(3,i_l+1,species(i_atom))
           end do
        end do

        ! TERM i_grid > 4 : Integral_i = Integral_(i-1) + h[9 f_i + 19 f_(i-1) - 5 f_(i-2) + f_(i-3)]/24
        do i_grid = 4, n_grid(species(i_atom)), 1

           if ( r_grid(i_grid,species(i_atom)) .lt. &
                (multipole_radius_free(species(i_atom))+extra_adding_to_hartree_potential_distance)) then
              !       radial integration weight on the logarithmic grid alpha*r times usual radial 
              !       integration weight from integral r^2 dr

              ! loop over all angular momentum components
              do i_l = 0, prec_max_l_species(species(i_atom)), 1
                 do i_m = - i_l, i_l, 1
                    
                    integral_zero_r(index_lm(i_m, i_l)) =   integral_zero_r(index_lm(i_m, i_l)) &
                       + (  9d0*angular_integral_log(index_lm(i_m, i_l),i_grid  )*green_I(i_grid  ,i_l+1,species(i_atom)) &
                         + 19d0*angular_integral_log(index_lm(i_m, i_l),i_grid-1)*green_I(i_grid-1,i_l+1,species(i_atom)) &
                         -  5d0*angular_integral_log(index_lm(i_m, i_l),i_grid-2)*green_I(i_grid-2,i_l+1,species(i_atom)) &
                         +  1d0*angular_integral_log(index_lm(i_m, i_l),i_grid-3)*green_I(i_grid-3,i_l+1,species(i_atom)))/24d0
                    
                    R_prec_multipole(i_grid,index_lm(i_m,i_l)) =  &                        ! spice with outside constant Green function 
                         integral_zero_r(index_lm(i_m, i_l) ) * green_K(i_grid,i_l+1,species(i_atom))
                    
                 end do
              end do
           else
              do i_l = 0, prec_max_l_species(species(i_atom)), 1
                 do i_m = - i_l, i_l, 1
                    R_prec_multipole(i_grid,index_lm(i_m,i_l)) =  &   
                        integral_zero_r(index_lm(i_m, i_l) ) * green_K(i_grid,i_l+1,species(i_atom))
                 end do
              end do
           end if
        end do        ! end radial loop for calculation of first integral

        ! second integral: from radius r to infinity - well, do it backwards and start by figuring out where the 
        !    integrand starts being non-zero, i.e. where we have to start worrying about it. 
        n_grid_limit = n_grid(species(i_atom))
        do while ( r_grid(n_grid_limit,species(i_atom)) .gt. &
                   (multipole_radius_free(species(i_atom))+extra_adding_to_hartree_potential_distance)) 
           n_grid_limit = n_grid_limit - 1
        end do
        
        ! start integrating from the outside in, again using the Adams-Moulton linear multistep integrator.
        ! the first three terms warrant special treatment, similar to the above integral. 
        integral_r_infty = 0d0
        do i_l = 0, prec_max_l_species(species(i_atom)), 1
           do i_m = - i_l, i_l, 1

              ! TERM 1 : Integral_N = h*f_N; but h = 1            
              integral_r_infty(index_lm(i_m, i_l)) = & 
                   angular_integral_log(index_lm(i_m, i_l),n_grid_limit) * green_K(n_grid_limit,i_l+1,species(i_atom)) 
              R_prec_multipole(n_grid_limit,index_lm(i_m,i_l)) =   &        
                   R_prec_multipole(n_grid_limit,index_lm(i_m,i_l)) & 
                   + integral_r_infty(index_lm(i_m, i_l) ) * green_I(n_grid_limit,i_l+1,species(i_atom))

              ! TERM 2 : Integral_(N-1) = Integral_N + h(f_(N-1)+f_N)/2 
              integral_r_infty(index_lm(i_m, i_l)) = integral_r_infty(index_lm(i_m, i_l)) + & 
                   ( 1d0*angular_integral_log(index_lm(i_m, i_l),n_grid_limit-1)*green_K(n_grid_limit-1,i_l+1,species(i_atom)) &
                   + 1d0*angular_integral_log(index_lm(i_m, i_l),n_grid_limit  )*green_K(n_grid_limit  ,i_l+1,species(i_atom)))/2d0
              R_prec_multipole(n_grid_limit-1,index_lm(i_m,i_l)) =   &        
                   R_prec_multipole(n_grid_limit-1,index_lm(i_m,i_l)) & 
                   + integral_r_infty(index_lm(i_m, i_l) ) * green_I(n_grid_limit-1,i_l+1,species(i_atom))

              ! TERM 3 : Integral_(N-2) = Integral_(N-1) + h(5f_(N-2) + 8f_(N-1) - f_N)/12
              integral_r_infty(index_lm(i_m, i_l)) = integral_r_infty(index_lm(i_m, i_l)) + & 
                   ( 5d0*angular_integral_log(index_lm(i_m, i_l),n_grid_limit-2)*green_K(n_grid_limit-2,i_l+1,species(i_atom)) &
                   + 8d0*angular_integral_log(index_lm(i_m, i_l),n_grid_limit-1)*green_K(n_grid_limit-1,i_l+1,species(i_atom)) &
                   - 1d0*angular_integral_log(index_lm(i_m, i_l),n_grid_limit  )*green_K(n_grid_limit  ,i_l+1,species(i_atom)))/12d0
              R_prec_multipole(n_grid_limit-2,index_lm(i_m,i_l)) =   &        
                   R_prec_multipole(n_grid_limit-2,index_lm(i_m,i_l)) & 
                   + integral_r_infty(index_lm(i_m, i_l) ) * green_I(n_grid_limit-2,i_l+1,species(i_atom))
          end do
        end do
        
        ! all remaining terms
        ! Integral_i = Integral_(i+1) + h[9 f_i + 19 f_(i+1) - 5 f_(i+2) + f_(i+3)]/24
        do i_grid = n_grid_limit-3, 1, -1
           do i_l = 0, prec_max_l_species(species(i_atom)), 1
              do i_m = - i_l, i_l, 1
              integral_r_infty(index_lm(i_m, i_l)) = integral_r_infty(index_lm(i_m, i_l)) + & 
                   (  9d0*angular_integral_log(index_lm(i_m, i_l),i_grid  )*green_K(i_grid  ,i_l+1,species(i_atom)) &
                   + 19d0*angular_integral_log(index_lm(i_m, i_l),i_grid+1)*green_K(i_grid+1,i_l+1,species(i_atom)) &
                   -  5d0*angular_integral_log(index_lm(i_m, i_l),i_grid+2)*green_K(i_grid+2,i_l+1,species(i_atom)) &
                   +  1d0*angular_integral_log(index_lm(i_m, i_l),i_grid+3)*green_K(i_grid+3,i_l+1,species(i_atom)))/24d0
                 R_prec_multipole(i_grid,index_lm(i_m,i_l)) =   &        
                      R_prec_multipole(i_grid,index_lm(i_m,i_l)) & 
                      + integral_r_infty(index_lm(i_m, i_l) ) * green_I(i_grid,i_l+1,species(i_atom))
              end do
           end do
        end do ! end calculation of second integral 
 
        ! set end of multipole expansion to zero explicitly, - to ensure that there is no diverging charge density
        ! outside the known grid quantities - this should be done on all atoms; this is done again on the finally splined multipole density later 

        ! VB: The following is an incorrect solution and therefore commented out. 
        !     If the spline really went past the grid boundary of the logarithmic grid, we
        !     would need to add the same logic as for the long-range part of the Hartree potential.
        !     The only thing that possibly saves us is the fact that the preconditioner
        !     is screened. This still needs to be checked.
        !
        !     And anyway: n_radial is the completelt wrong quantity to use. This particular
        !     quantity happens on the logarithmic grid, not the 'radial' grid.

              ! R_prec_multipole(n_radial(species(i_atom))+1,:) = 0.d0

        ! spline_atom_storage(n_atoms) is defined on each thread and contains the order of atoms which are treated by this 
        ! particular thread - in case that distributed spline_storage is required, the (atomic) index for a spline 
        ! in the spline storage array is given by spline_atom_storage(i_atom)
        
        ! do the splining itself: separately, for each atom, and each angular momentum shell. 
        ! the treatment of all radial shells is done in the splining routine
        do i_l = 0, prec_max_l_species(species(i_atom)), 1
           do i_m = - i_l, i_l, 1
              
              aux_R_prec_spline = 0d0
              call cubic_spline &
                   (R_prec_multipole(:,index_lm(i_m, i_l)), &
                   n_grid(species(i_atom)), &
                   aux_R_prec_spline )
              
              ! copy the array of spline coefficients to its eventual location ...
              do i_radial = 1, n_grid(species(i_atom)), 1
                 do i_spline = 1, n_max_spline, 1
                    R_multipole_spl( index_lm(i_m, i_l), i_spline, i_radial, current_spl_atom) &
                         = aux_R_prec_spline(i_spline,i_radial)
                 enddo
              enddo
           end do
        end do              ! end splinig loop
     end if                 ! atomic task list distribution in spline loop

  end do                    ! spline and integration loop over all atoms

  ! if every cpu keeps its own version of the splines, distribute them to 
  ! all the CPU's before doing much more ... 
  if (.not.use_distributed_spline_storage) &
       call sync_kerker_multipole_splines(R_multipole_spl,prec_l_dim,n_max_spline,n_spline_grid_dim,n_spline_atoms)

  ! FINALLY: sum up the solved and splined preconditioned residual from all different atoms
  !   over all grid points.
  !  need to pay attention to the memory location of the various splines.... 
  
  ! intialize variables: R_prec has not been used so far - I'd rather not want to know what it currently contains ... 
  current_spl_atom = 0
  R_prec           = 0d0

  ! for any given atom, sum over all centers that relate to this atom - especially important for PBC's where centers.ne.atoms in general
  ! this IS the exact same case as for the Hartree potential calculation, so we might as well use the variables initialized there. 

  do i_center = 1, n_centers_hartree_potential, 1

     prev_atom        = current_spl_atom
     current_center   = centers_hartree_potential(i_center)
     current_spl_atom = center_to_atom(current_center)                     ! this is the actual atom number, .ne. current_center for PBC's
     
     ! obtain the atomic spline for current center from where it is actually stored ... 
     if (current_spl_atom.ne.prev_atom) then
        ! distribute the information about that particular atom to all the threads
        ! simultaneously, each one is going to use it in order to compute its part 
        ! of the new residual ... 

        current_R_multipole_spl = 0d0
        call get_R_prec_spline_mpi(current_spl_atom,prec_l_dim,current_R_multipole_spl,R_multipole_spl, n_spline_grid_dim)

        ! determine maximal relevant radius for the center in question: only needs to be done for a new center,
        ! but it should be a function of l and m, after whose calculation the maximal value is picked out as well.          
        max_spl_rad_sq = 0d0
        do i_l = 0, prec_max_l_species(species(current_spl_atom))
           do i_m = -i_l, i_l

              ! start with the furthest point on each grid and move inwards until you find non-zero splines.
              i_grid = n_grid(species(current_spl_atom))+1

              ! check through splines ... see if they are equal to zero, decrease i_grid if they are 
              ! the check is representative when we square the multipole coefficients. N
              ! FIXME: one might also consider a threshold for delta below which 'there is no far field'
              delta = 0d0
              do while ((delta.le.multipole_thresh).and.(i_grid.gt.1))   ! FIXME: this should of course be a threshold 1E-12
                 i_grid = i_grid - 1
                 do i_spline = 1, n_max_spline
                    ! FIXME: this would be good enough if we only used the i_spline = 1 (or 0?) value !!!
                    delta = delta + current_R_multipole_spl(index_lm(i_m,i_l),i_spline,i_grid)**2
                 end do
              end do

              ! set maximal radius for the atomic preconditioned residual multipole decomposition 
              max_lm_spl_rad_sq(index_lm(i_m,i_l)) = r_grid(i_grid, species(current_spl_atom))**2

              if (max_lm_spl_rad_sq(index_lm(i_m,i_l)).gt.max_spl_rad_sq) then 
                 max_spl_rad_sq = max_lm_spl_rad_sq(index_lm(i_m,i_l))
              end if
           end do           
        end do

        
        ! output warning, atom number etc if the maximal spline radius is exactly equal to the outermost 
        !             grid radius

        if ((max_spl_rad_sq.gt.(multipole_radius_free(species(current_spl_atom)) & 
             +extra_adding_to_hartree_potential_distance)**2).and. &
             preconditioner_first_warnings) then

           ! FIXME: VB suggests that you should stop right here, but I don't want to do this until the routine is known to be working.
           write(info_str,'(2X,A)'   ) '* WARNING: Multipole expansion of preconditioned residual '
           call localorb_info(info_str,use_unit,'(A)')
           write(info_str,'(2X,A,I4)') '* extends too far at atom ', current_spl_atom
           call localorb_info(info_str,use_unit,'(A)')
        end if
        if (max_spl_rad_sq.gt.(multipole_radius_free(species(current_spl_atom)) & 
             +extra_adding_to_hartree_potential_distance)**2) then
           
             ! stop integrating there anyway to be in line with what the Hartree potential does
             max_spl_rad_sq = (multipole_radius_free(species(current_spl_atom))+extra_adding_to_hartree_potential_distance)**2
        end if

     end if

     ! FIXME: This is done without the long-range contribution in periodic cases, does that matter?
     ! what are the physical implementations of the 'Far-field' ? Does the filtering of low-frequencies
     ! enhance the regions which a given atom influences? - see comments in header of this file 
     ! loop over all BATCHES 
     i_full_points = 0
     do i_my_batch = 1, n_my_batches, 1

           do i_index = 1, batches(i_my_batch)%size, 1

              ! grid point index in the final preconditioned residual 
              i_full_points   = i_full_points + 1

              ! obtain distance^2 between integration point and current_center as well as a direction vector
              coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
              call  tab_single_atom_centered_coords_p0( current_center, coord_current, &
                   dist_tab_sq, dir_tab )

              if (dist_tab_sq.le.max_spl_rad_sq) then            ! is current grid point within the realm of current_spl_atom?

                 ! the sum to be done here is the same as Delley (1990) Eqn (12c), give or take integration factors:
                 ! R_prec(r) = sum_{alpha,l,m} Y_lm (r-R_alpha) R_prec_{alpha,l,m}(|r-R_alpha|)

                 ! trigonom_tab contains trig functions of polar angles in dir_tab:
                 call tab_single_trigonom_p0( dir_tab, trigonom_tab )

                 ! these can be used in the ylm tab calculation routine to get all necessary spherical harmonics:
                 ! INVESTIGATE/FIXME: this might be faster if we use the proper (and already known!) (l,m) dependence 
                 call tab_single_wave_ylm_p0(current_center,trigonom_tab, &
                      prec_max_l_species,precondition_max_l, ylm_tab )
   
                 ! obtain spline-interpolated values of the multipole components
                 ! of the preconditioned residual, splined on the logarithmic
                 ! integration grid
                 radius = sqrt(dist_tab_sq)    ! this is the actual distance from point to atom: spline this ... 

                 ! DEBUG: keep the actual radius for a moment, need it for output in a few lines
!                 delta = radius

                 ! need to invert on the logarithmic grid to find the index of the point in question and where we should be 
                 ! looking between the two grid points 
                 radius = invert_log_grid(radius,                                &
                                          r_grid_min(species(current_spl_atom)), &
                                          r_grid_inc(species(current_spl_atom)))

                 ! ... for all (m,l) values at once 
                 ! INVESTIGATE/FIXME: this might be faster if we use the proper (and already known!) (l,m) dependence !!!!!!!
                 call spline_vector                             &
                      ( radius,                                 &
                        current_R_multipole_spl,                &  
                        n_spline_grid_dim,                      &
                        prec_l_dim,                             &   ! FIXME: this is only necessary for some maximum l, should be faster
                        n_grid(species_center(current_center)), &
                        prec_l_dim_species(species_center(current_center)), &
                        aux_R_prec_result)

                 ! add contribution from this particular source to R_prec
                 R_prec(i_full_points) = R_prec(i_full_points) + & 
                              ddot ( prec_l_dim_species(species_center(current_center)), &
                              aux_R_prec_result(:), 1, ylm_tab(:), 1 ) 

              else

                 ! Far field: do nothing for the moment 

                 ! FIXME:
                 ! far-field calculation for grid point and current center goes in here if necessary

              end if    ! does the current integration point fall within the relevant preconditioning multipole radius???

           end do    ! loop over points in a given batch 
    
        !end if       ! checking if this CPU is responsible for current batch?

     end do          ! loop over batches

  end do          ! loop over relevant Hartree centers


  ! warnings should only be output once
  preconditioner_first_warnings = .false.

  ! where is the factor 4 pi, in the partition tab or not ???
  if (use_kerker_factor) then
    R_in(:) = R_in(:) -  kerker_coeff*precondition_kerker_q0**2*R_prec(:)
  else
    R_in(:) = R_in(:) -  precondition_kerker_q0**2*R_prec(:)
  endif
end subroutine precondition_kerker
!******

!------------------------------------------------------------------------------
!****s* preconditioner/bessel_I
!  NAME
!    bessel_I
!  SYNOPSIS

subroutine bessel_I(r,Lmax,y)

!  PURPOSE
! compute the modified spherical bessel functions I_(1/2)(r), I_(3/2)(r), ... , I_(Lmax+1/2)(r) and 
! store in the array y(0:Lmax) in that order. Routine *should* be stable for the purposes of AIMS.
! method: two different cases, see below
!  USES
  use constants
  implicit none

!  ARGUMENTS

  integer Lmax
  double precision r

!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!   o r    - argument of the functions
!   o Lmax - maximal angular momentum shell for which to compute
!  OUTPUT
!   o y(0:Lmax) - values of the Bessel functions
!  SOURCE


  integer Ndiff, i, j
  double precision y(Lmax+1),r_thresh,y_next,y_this,y_prev,buf
  double precision term_thresh, thresh_used
  parameter (Ndiff = 100)
  parameter (r_thresh = 5d-1)
  parameter (term_thresh = 1d-100)
  ! two cases: 
  if (r.ge.r_thresh) then
     ! large r >= r_thresh
     ! start Ndiff functions above Lmax and iterate downwards with the recursion formula 
     !     I(n-1) = I(n+1) + 2 n I(n)/r (Abramowitz, Stegun, Eqn 10.2.18)
     ! this gets you within a constant factor of the actual result - the latter is obtained
     ! by comparing with I(0) = sqrt(2 pi/r) sinh(r) which evaluates quite nicely
     !
     ! in my empirical investigation, this fixed-point type evaluation converged for r >~ 0.1; 
     ! r_thresh is set well above that value
     y_next = 1d-10
     y_this = 0d0
     do i = Lmax+Ndiff, Lmax, -1
        y_prev = y_next + 2*(dble(i+1)+5d-1)*y_this/r
        y_next = y_this
        y_this = y_prev
     end do
     y(Lmax+1) = y_next + 2*(dble(Lmax+1)+5d-1)*y_this/r
     if (Lmax.gt.0) then
        y(Lmax  ) = y_this + 2*(dble(Lmax)+5d-1)*y(Lmax+1)/r
        if (Lmax.gt.1) then
           do i = Lmax-2, 0, -1
              y(i+1) = y(i+3) + 2*(dble(i+1)+5d-1)*y(i+2)/r
           end do
        end if
     end if
     buf = sinh(r)*sqrt(2d0/(pi*r))/y(1)
     y(:) = buf*y(:)         
  else
     ! small r <= r_thresh
     ! direct implementation of Abramowitz, Stegun; Eqn (10.2.5).
     ! the summation is cut off when the terms in the series are smaller than 
     ! either r^6 or term_thresh (defined above), whichever is smaller. 
     thresh_used = r**6
     if (thresh_used.gt.term_thresh) thresh_used = term_thresh
     do i = 0, Lmax
        buf = 1d0
        y(i+1) = buf
        j    = 1
        do while (buf.gt.term_thresh) 
           buf = buf*r*r/(2d0*dble(j)*dble(2*i+2*j+1))
           y(i+1) = y(i+1) + buf
           j = j + 1
        end do
        buf = 1d0
        do j = 0, i
           buf = buf*(2d0*dble(j)+1d0)
        end do
        y(i+1) = r**(dble(i)+5d-1)*y(i+1)*sqrt(2/pi)/buf
     end do
  end if
end subroutine bessel_I
!******

!****s* preconditioner/bessel_K
!  NAME
!    bessel_K
!  SYNOPSIS
subroutine bessel_K(r,Lmax,y) 
!  PURPOSE
! compute the bessel functions I_(1/2)(r), I_(3/2)(r), ... , I_(Lmax+1/2)(r) and 
! store in the array y(0:Lmax) in that order. Routine *should* be stable for the purposes of AIMS.
! calculation is based on I_(1/2)(r) and I_(3/2)(r) which is then paired with the usual upwards iteration
!  USES
  use constants
  implicit none
! ARGUMENTS


  integer Lmax
  double precision r,y(Lmax+1)

!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o r    - argument of the functions
!    o Lmax - maximal angular momentum shell for which to compute
!  OUTPUT
!    y(0:Lmax) - values of the Bessel functions
!  SOURCE

  integer i

  y(1) = sqrt(pi/(2d0*r))*exp(-r)
  if (Lmax.ge.1) then
     y(2) = (1d0+1d0/r)*y(1)
     if (Lmax.ge.2) then
        do i = 1, Lmax-1
           y(i+2) = y(i)+2d0*(dble(i)+5d-1)*y(i+1)/r
        end do
     end if
  end if

end subroutine bessel_K

end module precondition
