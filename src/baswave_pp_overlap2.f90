!****s* FHI-aims/baswave_pp_overlap2
!  NAME
!    baswave_pp_overlap2
!  SYNOPSIS


!TODO: - when calling get_fnL_to_rowcol bas2rowcol is just the identity.

  subroutine baswave_pp_overlap2()

    !  PURPOSE
    !
    !    Calculates two-center integrals between all different basisfctns
    !    and pseudowavefctns.
    !
    !  USES
    use species_data
    use dimensions
    use sbt_overlap
    use basis
    use runtime_choices     ! for eg. N, lnr0 etc.
    use grids               ! n_grid etc.
    use mpi_tasks           ! check_allocation
    use logsbt
    use pseudodata
    use geometry
    use sbt_overlap_aims
    use prodbas

    implicit none

    !  ARGUMENTS


    !  INPUTS
    !    o pseudo_chi_spl, pseudo_grid_spl -- pseudo wave fctns on logarithmic grid
    !    o outer_radius, pp_outer_radius  -- outer radius at which the radial part has a
    !                                         magnitude of smaller than 1.0d-6 
    !  OUTPUTS
    !    o basiswav_pp_overlap(n_basis, n_pp_basis) -- overlap matrix
    !          --> (i_basis,:) : basis fctn
    !          --> (:,i_pp_basis) : pseudo valence fctn
    !              <basis(0)|pp_basis(Rvec)>
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE


!     local variables

    integer :: N, max_L
    real*8 :: lnr0, lnrange, lnk0
    real*8 :: power_bias_int
    real*8 :: Rabs_sq

    real*8,dimension(:,:), allocatable :: baswave      ! wave functions on finer logsbt grid
    real*8,dimension(:,:), allocatable :: pseudowave   ! pseudo wave functions on finer logsbt grid



!    integer, dimension(:,:,:), allocatable :: fnL_to_basfn
    integer, dimension(n_basis) :: basis2row
    integer, dimension(n_pp_basis) :: basis2col
    real*8, dimension(3) :: Rvec
!    integer, dimension(:,:), allocatable :: n_fnL
    integer :: pp_fnL_to_col(-max_pp_basis_l:max_pp_basis_l, max_n_pp_basis_fnLsp, 0:max_pp_basis_l )
    integer, dimension(:,:,:), allocatable :: fnL_to_row
    real*8, dimension(n_basis_fns) :: moment_1
    real*8, dimension(n_pp_basis_fns) :: moment_2
    real*8 :: times(4, SBT_N_TIMES)
    character(*), parameter :: func = 'basiswave_pp_overlap'
    character*130 :: info_str

    ! counter

    integer :: i, j
    integer :: i_pp_species, i_species
    integer :: i_atom, i_pp_atom, i_pp_basis_fns


    N = sbtgrid_N      ! N = 16384
    lnr0 = sbtgrid_lnr0
    lnrange = sbtgrid_lnrange
    lnk0 = sbtgrid_lnk0

    max_L = maxval(basis_l)


    moment_1 = 0.0d0
    moment_1 = 0.0d0    

    basis2row = (/(i, i=1,n_basis)/)        ! IDENTITY .. this needs to get modified if ( use_scalapack & use_packed_matrix) 
    basis2col = (/(i, i=1,n_pp_basis)/)         ! IDENTITY

    allocate(baswave(N,n_basis_fns))                 
    allocate(pseudowave(N,n_pp_basis_fns))            


!   start: import all pseudo wave fctns and basis wave fctns and
!          spline them on a finer-meshed grid (sbt_import_spline), 
!          and finally perform spherical Bessel transformation. 
!          (we need the fctns in k-space for the sbt_atomic_overlap routine!)
!
!          initially all wave fctns are available as u = r*f,
!          sbt_import_spline scales all functions u --> f,
!          logsbt_scale converts to f\tilde *r^{3/2} as needed for sbt_atomic_ovlp



    call sbt_atomic_field_transforms(N, lnr0, lnk0, lnrange, OVLP_TYPE_OVERLAP, &
   &                 n_basis_fns, basisfn_l, basisfn_species, basis_wave_spl, &
   &                 baswave, power_bias_int)



!TODO: sbt_atomic_field_transforms is not applicable for pseudowaves because loggridparameters are not able
!      to be parsed as variables

!-----------------------------------------------
    do i_pp_basis_fns = 1, n_pp_basis_fns

       i_pp_species = pp_basisfn_species(i_pp_basis_fns)

       ! spline -> sbtgrid
       call sbt_import_spline(N, pseudowave(:,i_pp_basis_fns), lnr0, lnrange, pp_basisfn_l(i_pp_basis_fns), &
       &                      n_points_pp_fn(i_pp_species), pseudo_chi_spl(:,:, i_pp_basis_fns), &
       &                      pp_r_grid_min(i_pp_species), pp_r_grid_inc(i_pp_species))



    end do



    ! ff == f(r)                     ! [assume full relative accuracy]
    call logsbt_scale(N, n_pp_basis_fns, pseudowave, lnr0, lnrange, 1.5d0)
    ! ff == f(r)*r**(3-power_bias)   ! [assume full relative accuracy]
    call logsbt_multi_driver(N, lnr0, lnk0, lnrange, 1.5d0, &
    &                        n_pp_basis_fns, pp_basisfn_l, pseudowave)


!
! --------- inquiry some array data needed for sbt_atomic_ovlp
!


     allocate(fnL_to_row(-max_pp_basis_L:max_pp_basis_L, max_n_basis_fnLsp, 0:max_pp_basis_L))


!------------ now we compute the overlaps <basis(O)|pp_basis(R_vec)>



    call initialize_fast_kernel(N, max_pp_basis_L + max_pp_basis_l, &
                                lnr0, lnk0, lnrange, power_bias_int)

    do i_atom = 1, n_pp_atoms

        i_species = species(i_atom)

        call get_fnL_to_rowcol(i_atom, &
                      n_pp_basis, pp_basis_atom, pp_basis_fn, pp_basis_l, pp_basis_m, &
                      max_pp_basis_l, max_n_pp_basis_fnLsp, Lsp2n_pp_basis_fnLsp(:,i_pp_species), & 
                      Lsp2_pp_basis_fn(:,:,i_pp_species), &
                   basis2row, fnL_to_row)


        do i_pp_atom = 1, n_pp_atoms

            i_pp_species = pp_species(i_pp_atom)

            call get_fnL_to_rowcol(i_pp_atom, &
                      n_pp_basis, pp_basis_atom, pp_basis_fn, pp_basis_l, pp_basis_m, &
                      max_pp_basis_l, max_n_pp_basis_fnLsp, Lsp2n_pp_basis_fnLsp(:,i_pp_species), & 
                      Lsp2_pp_basis_fn(:,:,i_pp_species), &
                      basis2col, pp_fnL_to_col)


!changed on 06/13/12: sbt_atomic_ovlp uses Rvec in the other orientation
!            Rvec = coords(:, i_atom) - pp_coords(:,i_pp_atom)
            Rvec = pp_coords(:,i_pp_atom) - coords(:, i_atom) 

            Rabs_sq = dot_product(Rvec,Rvec)

            Rvec = 0.0

!TODO:DB 10.11.: for huge systems saving basiswave_pp_overlap can become a bottle neck .. 
!                we should think about something to solve that, one day


            if (Rabs_sq < (maxval(outer_radius) + maxval(pp_outer_radius))**2) then

                call sbt_atomic_ovlp(N, lnr0, lnk0, lnrange, Rvec, 1D0, .false., &
                                 max_pp_basis_l, max_pp_basis_l, &
                                 max_n_pp_basis_fnLsp, max_n_pp_basis_fnLsp, &
                                 n_pp_basis_fns, pseudowave, pp_outer_radius, moment_2, & 
                                 Lsp2n_pp_basis_fnLsp(:,i_pp_species), Lsp2_pp_basis_fn(:,:,i_pp_species), fnL_to_row,&
                                 n_pp_basis_fns, pseudowave, pp_outer_radius, moment_2, & 
                                 Lsp2n_pp_basis_fnLsp(:,i_pp_species), Lsp2_pp_basis_fn(:,:,i_pp_species), &
                                 pp_fnL_to_col, basiswave_pp_overlap, times)


!            --> this provides the output overlap matrix entries
            endif
!          write(use_unit,*) 'basiswave_pp_overlap', basiswave_pp_overlap
        end do
    end do 

    call cleanup_fast_kernel()


    write(info_str,'(2X,A,F5.3,2X,A)') 'Extension of Kleinman-Bylander projectors:     ',&
                                  maxval(pp_outer_radius),'bohr'
    call localorb_info(info_str,use_unit,'(A)', OL_low) 

    write(info_str,'(2X,A,F6.3,2X,A)') 'Projection is zero for distance larger than   ',&
                                  maxval(outer_radius) + maxval(pp_outer_radius),'bohr'
    call localorb_info(info_str,use_unit,'(A)', OL_low) 

    write(info_str,*)
    call localorb_info(info_str,use_unit,'(A)', OL_low)


! ERASE AFTER DEBUGGING
!     if(myid==0)  write(use_unit,*) 'baswave_pp_overlap: maxval(basiswave_pp_overlap) =', maxval(basiswave_pp_overlap)
!     if(myid==0)  write(use_unit,*) 'baswave_pp_overlap: minval(basiswave_pp_overlap) =', minval(basiswave_pp_overlap)
  

  write(use_unit,*) 'basiswave_pp_overlap()', basiswave_pp_overlap


  if (allocated(baswave)) then
      deallocate(baswave)
  end if

  if (allocated(pseudowave)) then
      deallocate(pseudowave)
  end if


  stop

  end subroutine baswave_pp_overlap2

