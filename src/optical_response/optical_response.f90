subroutine optical_response()

!  PURPOSE
!
!  This subroutine is called once after SCF and scaledZORA are complete
!  and computes the optical resonse of a given molecule.
!
!  ONLY usable for NON-periodic structures. (for now)
!
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!    This file (and optics_module) was written by Jan Kloppenburg
!  SOURCE

  use localorb_io
  use runtime_choices
  use dimensions
  use physics
  use prodbas
  use gw_para
  use hartree_fock
  use species_data
  use optics_module
  use free_atoms
  use timing

  implicit none

  integer :: errnum, max_pairs, k_cnt, count_electrons, i_spin, i_k_point

  ! Time total
  call get_timestamps(time_total_neutral_excitation, &
        clock_time_total_neutral_excitation)

  occ_states = 0
  unocc_states = 0

  do i_state=1, n_states
    count_electrons = 0d0
    do i_spin = 1, n_spin
      count_electrons = count_electrons +  occ_numbers(i_state,i_spin,n_k_points)
    end do
    if( count_electrons.le.1d-6 ) then
      unocc_states = unocc_states + 1
    else
      occ_states = occ_states + 1
    endif
  enddo

  max_pairs = occ_states * unocc_states

  occ_limit = 0
  unocc_limit = 0

  if(casida_reduce_matrix) then
    do i_state=1, occ_states
      if(ks_eigenvalue(i_state,n_spin,n_k_points).lt.casida_occ_limit) occ_limit = occ_limit+1
    enddo
    do i_state=occ_states+1, n_states
      if(ks_eigenvalue(i_state,n_spin,n_k_points).gt.casida_unocc_limit) unocc_limit = unocc_limit+1
    enddo
  else
    occ_limit = 1
    unocc_limit = 0
  endif

  if(.not.use_hartree_fock) then
    if(myid==0) then
      write(use_unit,'(2x,a)') '| Initializing the product basis'
    endif
    call initialize_prodbas()
    if(myid==0) then
      write(use_unit,*)
      write(use_unit,'(2x,a)') '| Product basis set up.'
    endif
  endif


  if(myid==0) then
    write(use_unit,*) ' '
    write(use_unit,'(2x,a)') '|'
    write(use_unit,'(2x,a)') '| Starting routine for neutral optical excitation energies.'
    write(use_unit,'(2x,a)') '|'
    write(use_unit,'(2X,a,i3.1,a,i4.1,a,a,i8.1,a)') '| We have ',occ_states,' occupied and ',&
      unocc_states,' virtual states, ','yielding ',occ_states*unocc_states*n_spin,' paired states.'
    if(casida_reduce_matrix) then
      write(use_unit,'(2x,a)') '|'
      write(use_unit,'(2x,a,f9.2,a,i4.1)') '| Number of occupied states removed ( E < ',casida_occ_limit,'Ha): ', occ_limit
      write(use_unit,'(2x,a,f9.2,a,i4.1)') '| Number of  virtual states removed ( E > ',casida_unocc_limit,'Ha): ', unocc_limit
      write(use_unit,'(2x,a)') '|'
      write(use_unit,'(2x,a,i4.1)') '| Remaining occupied states: ', occ_states-occ_limit
      write(use_unit,'(2x,a,i4.1)') '| Remaining virtual states: ', unocc_states-unocc_limit
      write(use_unit,'(2x,a)') '|'
      write(use_unit,'(2x,a,i14.1,a,i9.1,a)') '| The number of paired states was reduced from ', max_pairs, ' to ', &
        ((occ_states-occ_limit)*(unocc_states-unocc_limit)), ' .'
    endif
    if(num_excited_states.eq.-1) then
      write(use_unit,'(2x,a)') '| The number of states to be printed was not specified.'
      num_excited_states = min(20,occ_states*unocc_states*n_spin)
      write(use_unit,'(2x,a,i4.1)') '| Defaulting to ', num_excited_states
    elseif(occ_states*unocc_states*n_spin.lt.num_excited_states) then
      write(use_unit,'(2x,a)') '| The number of requested states is larger than the number of paired states.'
      write(use_unit,'(2x,a)') '| Reducing the printed states to the maximum.'
      num_excited_states = occ_states*unocc_states*n_spin
      write(use_unit,'(2x,a,i4.1,a)') '| There will be ', num_excited_states, ' printed.'
    else
      write(use_unit,'(2x,a,i4.1,a)') '| You requested ', num_excited_states, ' to be printed.'
    endif
    write(use_unit,'(2x,a)') '|'
  endif

  if(casida_reduce_matrix) then
    n_pair_states = (occ_states-occ_limit)*(unocc_states-unocc_limit)
    low_limit = occ_limit+1
    if(occ_limit.eq.0) low_limit=1
  else
    n_pair_states = occ_states * unocc_states
    low_limit = 1
    unocc_limit = 0
  endif

  ! Set up "mostly" square processor grid
  do num_pcol = nint(dsqrt(dble(n_tasks))),2,-1
    if(mod(n_tasks,num_pcol)==0) exit
  enddo
  num_prow = n_tasks/num_pcol

  ! Initialize Blacs grid
  my_casida_ctxt = mpi_comm_world
  call blacs_gridinit( my_casida_ctxt, 'C' , num_prow, num_pcol )
  call blacs_gridinfo( my_casida_ctxt, num_prow, num_pcol, my_prow, my_pcol )

  num_blocks = min(n_tasks,floor(dble(n_pair_states)/dble(n_tasks)))

  na_rows = numroc(n_pair_states, num_blocks, my_prow, 0, num_prow)
  na_cols = numroc(n_pair_states, num_blocks, my_pcol, 0, num_pcol)

  call descinit( casida_desc, n_pair_states, n_pair_states, &
                 num_blocks, num_blocks, 0, 0, my_casida_ctxt, na_rows, errnum )

  call blacs_gridinfo(casida_desc(2), nprow, npcol, myrow, mycol)

  allocate(coord_of_center(3),stat=errnum)
  call check_allocation(errnum,' coord_of_center allocate')
  allocate(atom_coord_wrt_center(3,n_atoms),stat=errnum)
  call check_allocation(errnum,' atom_coord_of_center allocate')
  call determine_center_of_molecule(coord_of_center, atom_coord_wrt_center)

  if(n_spin==1) then
    ! prepare for spin
    i_spin = 1
    ! prepare for k_points
    i_k_point = 1
    if(n_tasks>1) then
!    if(use_2d_corr) then ! JK: apperently, the internal use of this flag has changed,
                          !     I use n_tasks>1 now, that should be safe to determine
                          !     whether we are running on more than one core

      ndim1_o3KS = (n_states-1)/np1_o3ks + 1
      ndim2_o3KS = (n_states-1)/np2_o3ks + 1

      if (.not.allocated(ovlp_3ks)) then
         allocate(ovlp_3ks(n_basbas,ndim1_o3ks,ndim2_o3ks,n_spin), stat=errnum)
         call check_allocation(errnum, 'ovlp_3ks (2D)')
      endif
      ovlp_3ks = 0d0

      allocate(dipole_moment(occ_states,unocc_states,n_spin,3),stat=errnum)
      call check_allocation(errnum,'casida dipole_moment')
      call integrate_fxc_pairstates(occ_states, unocc_states, coord_of_center, &
              l_shell_max, KS_eigenvector(:,:,:,i_k_point), dipole_moment)

      if(neutral_excitation_tddft) then
        allocate(fxc_3ks(n_basbas,ndim1_o3ks,ndim2_o3ks,n_spin), stat=errnum)
        call check_allocation(errnum, 'fxc_3ks (2D)')
        allocate(pure_3ks(n_basbas,ndim1_o3ks,ndim2_o3ks,n_spin), stat=errnum)
        call check_allocation(errnum, 'pure_3ks (2D)')
        fxc_3ks = 0d0
        pure_3ks = 0d0
        call get_tddft_orbitals_2d()
        if(excited_mode_singlet.and.excited_mode_triplet) then
         excited_mode_singlet = .true.
         excited_mode_triplet = .false.
         call casida_tddft_2d()
         excited_mode_singlet = .false.
         excited_mode_triplet = .true.
         call casida_tddft_2d()
        else
         call casida_tddft_2d()
        endif
      endif

      if(neutral_excitation_rpa) then
        call get_ks_orbitals_2d()
        if(excited_mode_singlet.and.excited_mode_triplet) then
         excited_mode_singlet = .true.
         excited_mode_triplet = .false.
         call casida_rpa_2d()
         excited_mode_singlet = .false.
         excited_mode_triplet = .true.
         call casida_rpa_2d()
        else
         call casida_rpa_2d()
        endif
      endif

      if(neutral_excitation_tdhf) then
        call aims_stop('TDHF calculation in MPI mode not yet possible due to the lack of a non-Hermitean eigensolver')
       !call get_ks_orbitals_2d()
       !call casida_tdhf_2d()
      endif
    else

      if (.not.allocated(ovlp_3ks)) then
        allocate(ovlp_3ks(n_basbas,n_states,n_states,n_spin), stat=errnum)
        call check_allocation(errnum, 'ovlp_3ks (1D)')
      endif

      allocate(dipole_moment(occ_states,unocc_states,n_spin,3),stat=errnum)
      call check_allocation(errnum,'casida dipole_moment')
      call integrate_fxc_pairstates(occ_states, unocc_states, coord_of_center, &
              l_shell_max, KS_eigenvector(:,:,:,i_k_point), dipole_moment)

      if(neutral_excitation_tddft) then
       allocate(fxc_3ks(n_basbas,n_states,n_states,n_spin), stat=errnum)
       call check_allocation(errnum, 'fxc_3ks (1D)')
       allocate(pure_3ks(n_basbas,n_states,n_states,n_spin), stat=errnum)
       call check_allocation(errnum, 'pure_3ks (1D)')
       call get_tddft_orbitals_1d()
       if(excited_mode_singlet.and.excited_mode_triplet) then
        excited_mode_singlet = .true.
        excited_mode_triplet = .false.
        call casida_tddft()
        excited_mode_singlet = .false.
        excited_mode_triplet = .true.
        call casida_tddft()
       else
        call casida_tddft()
       endif
       deallocate(fxc_3ks)
       deallocate(pure_3ks)
      endif

      if(neutral_excitation_rpa) then
       call get_ks_orbitals_1d()
       if(excited_mode_singlet.and.excited_mode_triplet) then
        excited_mode_singlet = .true.
        excited_mode_triplet = .false.
        call casida_rpa()
        excited_mode_singlet = .false.
        excited_mode_triplet = .true.
        call casida_rpa()
       else
        call casida_rpa()
       endif
      endif

      if(neutral_excitation_tdhf) then
       call get_ks_orbitals_1d()
       if(excited_mode_singlet.and.excited_mode_triplet) then
        excited_mode_singlet = .true.
        excited_mode_triplet = .false.
        call casida_tdhf()
        excited_mode_singlet = .false.
        excited_mode_triplet = .true.
        call casida_tdhf()
       else
        call casida_tdhf()
       endif
      endif
    deallocate(ovlp_3ks)
    endif

  else

    call aims_stop('************   Spin UNrestricted calculations are not yet possible!')

  endif

  call get_times(time_total_neutral_excitation,&
        clock_time_total_neutral_excitation)

  if(myid==0) then
    write(use_unit,'(2x,a)') '|'
    write(use_unit,'(2x,a)') '| End of neutral optical excitation energy calculation routine.'
    write(use_unit,'(2x,a)') '|'
  endif

end subroutine
