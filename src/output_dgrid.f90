!  NAME: output_dgrid
!  SYNOPSIS
!!!
!!!!!  **   written by Alim Ormeci - MPI CPfS Dresden --- Aug. 2018
!!!!   **   revised by CC // Nov 2018
subroutine output_dgrid ()
  !  PURPOSE
  !  Print out necessary information to enable ELI, electron density, etc.
  !  calculations and their analysis through DGRID, a separate external 
  !  program package.
  !
  !  USES
  !
  use dimensions
  use localorb_io
  use runtime_choices 
  use physics
  use synchronize_mpi
  use mpi_tasks
  use scalapack_wrapper
  use pbc_lists
  use geometry
  use species_data
  use basis
  use grids
  use aims_memory_tracking
  implicit none

  !
  ! local variables
  integer                                    :: i, j, k, indx, isp, nmaxtotal
  real*8, dimension(:,:,:), allocatable      :: dummy_eigvec
  real*8, dimension(:), allocatable          :: EvecR
  complex*16, dimension(:,:,:), allocatable  :: dummy_eigvec_cmplx
  complex*16, dimension(:), allocatable      :: EvecI
  character(*), parameter :: func    = 'output_dgrid'
  character(*), parameter :: version = '1.0'
  character(80) :: info_txt
  logical warn_msg

  ! Greet the user
  call localorb_info(' DGrid output files are now being written...')

  ! Initialize the warning flag to false 
  ! This flag may be necessary in some hybrid functional calculations
  warn_msg = .false.

  ! Work mostly on Proc 0:
  if ( myid == 0 ) then

    ! Open file:
    open(unit=67,file='dgrid_aims.dat',ACTION='WRITE')
    write(67, '( a, a )')  'FHI-aims - DGrid Interface Version  ', version

    ! check if sum of species radial function numbers equal to max. radial fn. number
    ! if additional aux. basis functions are enabled in the control.in file as 
    ! suggested for methods using "resolution of identity", then we have to correct
    ! the value of the max. number of radial functions for DGrid to work correctly
    nmaxtotal = Sum( n_basis_fn_species(1:n_species) )
    if ( nmaxtotal /= n_max_basis_fns ) then
       warn_msg = .true.
    end if

    ! print out the values of the basic integer variables regarding the number of
    ! different elements (species), atoms, states, basis functions, spin etc.
    write(67,'( 8i6 )') n_species, n_atoms, n_centers,  n_states, n_basis, &
      & nmaxtotal, n_periodic, n_spin

    ! total number of k points; k-mesh; highest value of L used in the basis set
    write(67,'(5i6)') n_k_points, n_k_points_xyz(:), max_basis_L

    ! number of grid points for each species
    write(67,*) (n_grid(i), i = 1, n_species)

    ! coordinates of the atoms in the molecule/cluster or in the unit cell
    do i = 1, n_atoms
      j = species(i)
      write(67,'(i3,2x,3(f20.15,1x),1x,(a))') j,coords(:,i),trim(species_name(j))
    end do

    ! energy and chemical potential
    write(67,'('' Total energy in Ha units:'', 2x, f25.9)') total_energy
    write(67,'('' Chemical pot in Ha units:'', 2x, f25.9)') chemical_potential

    ! information regarding the basis functions
    do j = 1, n_basis_fns
      write(67,'( 3(i5), f20.15 )') basisfn_species(j), basisfn_l(j), &
        & basisfn_n(j), outer_radius(j)
    end do
    write(67,*) ( n_basis_fn_species(i), i = 1, n_species )

    ! radial grid information for each species: r_min and incremental value
    do i = 1, n_species
      write(67,*) r_grid_min(i), r_grid_inc(i)
    end do

    ! spline coefficients for each basis function
    do i = 1, n_basis_fns
      k = basisfn_species(i)
      do j = 1, n_grid(k)
        write(67,'( 4(es28.20) )') basis_wave_spl(:,j,i)
      end do
      write(67,*)   
    end do

    ! eigenvalues and their occupancies
    do k = 1, n_k_points
      do j = 1, n_spin
        do i = 1, n_states
          write(67,'(2(es28.20,4x))') KS_eigenvalue(i,j,k), occ_numbers(i,j,k)
        end do
      end do
      write(67,*)    
    end do
  end if ! End Proc 0 part

  ! CC: For eigvecs we need MPI, writes on Proc 0 are inside the loop
  ! Case with Real-Eigenvectors first:
  if (real_eigenvectors) then
     if ( myid == 0 ) then
        write(67,*) -5, '  case of real eigenvectors'
     end if
    ! Local, temporary real. Evec for output
    call aims_allocate( EvecR, n_basis, "DGRid_EvecR" ) 

    ! Scalapack Case
    if (use_scalapack) then
      ! Loop over kpoints, spins, states
      do k = 1, n_k_points
        do isp = 1, n_spin
          do i = 1, n_states
            ! Now sync and write on Proc 0
            call sync_single_eigenvec_scalapack(EvecR,i,isp,k)
            if ( myid == 0 ) then
              write(67,'( 2(es30.22,1x) )') (EvecR(j), j = 1, n_basis)
              write(67,*)
            end if
          end do
        end do
      end do !Finished
   
    else ! Lapack Case
      ! Full Evec Array
      call aims_allocate( dummy_eigvec, n_basis, n_states, n_spin, "DGRid_dummy_eigvec" )
      indx = 0
      do k = 1, n_k_points
        if (myid == MOD(k, n_tasks) .and. myid <= n_k_points ) then
          indx = indx + 1
          dummy_eigvec(:,:,:) = KS_eigenvector(:,:,:,indx)
        else
          dummy_eigvec = 0.0d0
        end if
        do isp = 1, n_spin
          do i = 1, n_states
            ! Now sync and write on Proc 0
            EvecR(:) = dummy_eigvec(:,i,isp)
            call sync_vector(EvecR,n_basis)
            if ( myid == 0 ) then
              write(67,'( 2(es30.22,1x) )') (EvecR(j), j = 1, n_basis)
              write(67,*)
            end if
          end do
        end do ! spin loop
      end do ! k-point loop
      call aims_deallocate( dummy_eigvec, "+DGRid_dummy_eigvec" )
    end if
    ! Deallocate array
    call aims_deallocate( EvecR, "+DGRid_EvecR" ) ! Needed?

  ! ...now come the complex eigenvectors
  else
     if ( myid == 0 ) then
        write(67,*)  5, '  case of complex eigenvectors'
     end if
    ! Local, temporary cmpl. Evec for output
    call aims_allocate( EvecI, n_basis, "DGRid_EvecI" )

    ! Scalapack Case
    if (use_scalapack) then
      ! Loop over kpoints, spins, states
      do k = 1, n_k_points
        do isp = 1, n_spin
          do i = 1, n_states
            ! Now sync and write on Proc 0
            call sync_single_eigenvec_scalapack_complex(EvecI,i,isp,k)
            if ( myid == 0 ) then
              write(67,'( 2(es30.22,1x) )') (EvecI(j), j = 1, n_basis)
              write(67,*)
            end if
          end do
        end do
      end do !Finished

    else ! Lapack Case
      ! Full Evec Array
      call aims_allocate( dummy_eigvec_cmplx, n_basis, n_states, n_spin, "DGRid_dummy_eigvec_cmplx" )
      indx = 0
      do k = 1, n_k_points
        if (myid == MOD(k, n_tasks) .and. myid <= n_k_points ) then
          indx = indx + 1
          dummy_eigvec_cmplx(:,:,:) = KS_eigenvector_complex(:,:,:,indx)
        else
          dummy_eigvec_cmplx = (0.0d0,0.0d0)
        end if
        do isp = 1, n_spin
          do i = 1, n_states
            ! Now sync and write on Proc 0
            EvecI(:) = dummy_eigvec_cmplx(:,i,isp)
            call sync_vector_complex(EvecI,n_basis)
            if ( myid == 0 ) then
              write(67,'( 2(es30.22,1x) )') (EvecI(j), j = 1, n_basis)
              write(67,*)
            end if
          end do
        end do ! spin loop
      end do ! k-point loop
      call aims_deallocate( dummy_eigvec_cmplx, "+DGRid_dummy_eigvec_cmplx" )
    end if
    ! Deallocate array
    call aims_deallocate( EvecI, "+DGRid_EvecI" )
  end if

  ! Periodic Systems // again work on Proc 0
  if ( n_periodic > 0 ) then
    ! real and reciprocal space lattice vectors
    if ( myid == 0 ) then 
      do i = 1, n_periodic
        write(67,'( 3f22.16)') lattice_vector(:,i)
      end do
      do i = 1, n_periodic
        write(67,'( 3f22.16)') recip_lattice_vector(:,i)
      end do
      ! list of k-points and their weights
      do i = 1, n_k_points
        write(67,'( 3(f19.15), 1x, f20.16)') k_point_list(i,:), k_weights(i)
      end do
    end if ! myid = 0 check
  end if

  ! Close file:
  close(67)

  ! Say Goodbye
  if ( warn_msg ) then
     write(info_txt, '(a, i6)') 'INFO: number of radial functions         = ', n_max_basis_fns 
     call localorb_info(info_txt)
     write(info_txt, '(a, i6)') '      differs from sum of species values = ', nmaxtotal
     call localorb_info(info_txt)
     call localorb_info('      This is only OK if auxil. basis functions have been used.')
     call localorb_info('      See file control.in!')
  end if
  call localorb_info('  Finished writing the DGrid output file (dgrid_aims.dat).')
  call localorb_info('')
  
end subroutine output_dgrid
