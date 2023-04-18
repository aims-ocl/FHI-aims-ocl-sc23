!****s* FHI-aims/integrate_v_times_radialwaves
!  NAME
!   integrate_v_times_radialwaves
!  SYNOPSIS

      subroutine integrate_v_times_radialwaves &
      ( v_times_radialwaves_spl)

!  PURPOSE
!  Subroutine integrate_v_times_radialwaves
!
!  computes the integration over the Coulomb potential and the radial wave
!  function:
!    v_times_radial_wave (r) = \int r'^2 dr' u_il (r')/r'|r-r'|
!     = 4*pi/(2*l+1) {\int_0^r r'2 dr' u_il(r')/r' r'^l/r^(l+1)  +
!                 \int_r^infty r'2 dr' u_il(r')/r' r^l/(r'^(l+1) }
!
!  This will be used in the evaluation of the coulomb interaction matrix
!    <i| v(r-r')|j>.
!
!  USES

      use dimensions
      use basis
      use runtime_choices
      use grids
      use geometry
      use species_data
      use spline
      use prodbas
      use constants
      use mpi_tasks
      use basbas_fn_coulomb
      use localorb_io, only : use_unit
      implicit none

!  ARGUMENTS
      real*8, dimension(n_max_spline, n_hartree_grid, n_loc_prodbas) &
               ::  v_times_radialwaves_spl

!  INPUTS 
!   o none  
!  OUTPUTS
!   o v_times_radialwaves_spl -- the integration over Coulomb interaction v
!         and the radial part of the auxiliary basis function  
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
!  variables for hse_correction_term & testvariables
      real*8, dimension(:,:,:,:), allocatable :: &
                hse_matrix
      real*8 :: hse_correction_term
      integer :: max_prodbas_l
      integer :: j_grid
      integer :: n_grid_inf


      real*8, dimension(:), allocatable :: v_times_radial_waves

      real*8, dimension(:,:), allocatable ::   radial_waves


      real*8  :: i_r

      real*8  :: r_l
      real*8  :: r_neg_l1

!   coefficients for the integration weight on logrithmic grid.
      real*8  :: alpha
!   integration weight  r^2 * alpha*r
      real*8  :: dr_coef
      real*8 :: multipole

!  counters

      integer i_loc_prodbas
      integer i_atom, i_species, i_basbas_fn
      integer i_grid
      integer i_l
      integer i_task
      integer i_basbas
      integer i_radial
      integer j_radial
      integer j_test
      integer i_order

      integer :: info
      character(*), parameter :: func = 'integrate_v_times_radialwaves'
!  begin work


      allocate(radial_waves(n_hartree_grid, n_loc_prodbas))
      allocate(v_times_radial_waves(n_hartree_grid))

!
! start allocate variables for hse_correction_term
!
      if ((use_hse.and.hse_omega_hf.ne.0.0d0.and..not. use_gw_and_hse)&
         .or. lrc_pt2_started) then
         if (.not. use_logsbt_for_radial_hse_integration) then
            max_prodbas_L = maxval(basbas_l)
            allocate(hse_matrix(n_hartree_grid, n_hartree_grid, &
            &                   max_prodbas_L+1, n_species), stat=info)
            call check_allocation(info, 'hse_matrix', func)
            do i_species = 1, n_species
               call integrate_errorfunction(max_prodbas_l, &
               & n_grid(i_species), n_hartree_grid, r_grid(:, i_species), &
               & hse_omega_hf, hse_matrix(:,:,:, i_species)) 
            end do
         end if
      endif

! end allocating errorfunction expansion


      do i_task = 1, n_tasks
        do i_loc_prodbas = 1, n_loc_prodbas

           i_basbas=map_prodbas(i_loc_prodbas, i_task)
           if(i_basbas.gt.0 .and. myid.eq.i_task-1) then

              i_atom = basbas_atom(i_basbas)
              i_l = basbas_l(i_basbas)
              i_species = species(i_atom)
              i_basbas_fn = basbas_fn(i_basbas)
              v_times_radial_waves = 0.d0

              if (force_hartree_log_grid) then

                 do i_grid = 1, n_grid(i_species)
                    radial_waves (i_grid, i_loc_prodbas) = &
                    & basbas_wave_spl( 1, i_grid, i_basbas_fn)
                 enddo

                 if ((.not.(use_hse.and.(hse_omega_hf.ne.0.0d0)) .or. &
                     use_gw_and_hse) .and. .not. lrc_pt2_started) then
                    if(Adams_Moulton_integrator)then

                       call adams_moulton_wave_integrator( &
                       & i_l, radial_waves(:, i_loc_prodbas), &
                       & n_grid(i_species), r_grid(:, i_species), &
                       & v_times_radial_waves, multipole)

                    else

                       call trapezoidal_wave_integrator( &
                       & i_l, radial_waves(:, i_loc_prodbas), &
                       & n_grid(i_species), r_grid(:, i_species), &
                       & v_times_radial_waves, multipole)

                    endif

                 else

                    if (use_logsbt_for_radial_hse_integration) then

                       call hse_logsbt_integrator_grid(i_l, &
                       & basbas_wave_spl(:,:,i_basbas_fn), n_grid(i_species), &
                       & r_grid_min(i_species), r_grid_inc(i_species), &
                       & hse_omega_hf, v_times_radial_waves)

                    else
                       
                       call hse_wave_integrator(i_l, &
                       & radial_waves(:, i_loc_prodbas), n_grid(i_species), &
                       & n_hartree_grid, r_grid(:, i_species), &
                       & hse_matrix(:,:, i_l+1, i_species), &
                       & v_times_radial_waves)

                    end if

                 endif

                 call cubic_spline &
                 ( v_times_radial_waves, n_grid(i_species), &
                 v_times_radialwaves_spl( :, :,  i_loc_prodbas))

                 !  obtain the outer radius for v_times_radialwaves
                 i_grid = n_grid(i_species)
                 do while ((abs(v_times_radial_waves(i_grid)) .le. &
                 &          wave_threshold) .and. &
                 &         (i_grid.gt.1))
                    i_grid = i_grid-1
                 enddo
                 if(i_grid.lt.1) then
                    write(use_unit,'(1X,A,A)') &
                    "* Warning - a basis function is lower that the ", &
                    "requested threshold value for integrations everywhere."
                    write(use_unit,'(1X,A,A)') &
                    "Species : ", species_name(i_species)
                    write(use_unit,'(1X,A,I3)') &
                    "(l)   : ", basbas_l(i_basbas)
                    stop
                 end if

              else  ! radial grid

                 call radial_wave_integrator(i_l, n_radial(i_species), &
                 & r_radial(:, i_species), w_radial(:, i_species), &
                 & n_grid(i_species), &
                 & r_grid_min(i_species), r_grid_inc(i_species), &
                 & basbas_wave_spl(:,:, basbas_fn(i_basbas)), &
                 & v_times_radial_waves)

                 call cubic_spline &
                 ( v_times_radial_waves, &
                 n_radial(species(i_atom))+1, &
                 v_times_radialwaves_spl( :, :,  i_loc_prodbas))

              end if


             endif

! end of i_loc_prodbas
        enddo
! end of loop i_task
      enddo
!      if(myid.eq.0) then
!      open(105,file='radial_integral.dat')
!      do i_grid=1, n_grid(1), 1
!       write(use_unit,'(2X,2f16.10)')r_grid(i_grid,1),v_times_radialwaves_spl(1,i_grid,10)
!       write(105,'(2X,51f16.10)') r_grid(i_grid,1), &
!                  (v_times_radialwaves_spl(1,i_grid,i_basbas), i_basbas=1, 50) 
!       enddo
!      endif

      deallocate (v_times_radial_waves)
      deallocate (radial_waves)
      if(allocated(hse_matrix))then
         deallocate(hse_matrix)
      endif
      return

      end subroutine integrate_v_times_radialwaves
!******
