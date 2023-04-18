!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Parent subroutine for performing magnetic response (MR)
!!  calculations. It is responsible for some initializations, printing
!!  output, and calling MR_core which does the computationally heavy
!!  parts. MR_core is called one or more times. Current functionality
!!  includes the NMR shieldings, J-couplings, magnetizability, and the
!!  electric field gradient. See MR_core.f90 for more details.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
!!  COMMENTS
!!
!!  See MR_global.f90 for the description of the main MR work
!!  variables.
!!
subroutine MR_main()

  use aims_memory_tracking,  only: push_current_memory
  use constants,             only: pi
  use DFPT,                  only: reset_DFPT_timers
  use fjson_datatype,        only: fjson_handle
  use fjson_rw,              only: fjson_finish_object
  use geometry,              only: coords, isotope_numbers, mr_custom_mus, &
       & mr_custom_spins, species
  use json_output,           only: get_json_handle
  use localorb_io,           only: localorb_multi
  use mpi_tasks,             only: myid
  use MR_global
  use MR_output
  use nuclear_data
  use integration,           only: cleanup_integration, initialize_integration
  use tools,                 only: newunit, start_wall, stop_wall, str, &
       & walltime_mat_conv
  use runtime_choices,       only: mr_gauge_origin, out_aims_json_log, &
       & output_sxml, use_scalapack
  use species_data,          only: species_name, species_m
  use synchronize_mpi_basic, only: sync_vector
  use timing,                only: clock_time_magnetic_response
  use types,                 only: dp

  implicit none

  ! Information about isotopes present in the system
  type(nucleus_t), allocatable :: MR_nuclei(:)
  real(dp) :: cputime_tmp1, cputime_tmp2
  ! Auxiliary variables for 'where' statements needed for the PGI
  ! compiler.
  character(3), allocatable :: names_tmp(:)
  real(dp), allocatable :: values_tmp(:)
  ! Other
  type(fjson_handle) :: json_log_handle
  integer :: i_atom, i_counter, i_dir, i_term

  ! STEP 1 - Prepare the environment based on input from control.in
  !          and geometry.in. Some of that was already done in
  !          check_consistency_of_keywords.
  allocate(c_atoms(size(active_nuclei)))
  do i_atom = 1, size(active_nuclei)
     c_atoms(i_atom) = 'atom '//str(active_nuclei(i_atom))//' ('// &
          & trim(species_name(species(active_nuclei(i_atom))))//')'
  end do
  ! Unless explicitly given in control.in, set the gauge origin to the
  ! center of mass.
  default_not_changed: if (sum(mr_gauge_origin) > 1d299) then
     gauge_origin = [(sum(species_m(species)*coords(i_dir,:)), i_dir=1,3)]/ &
          & sum(species_m(species))
  else
     gauge_origin = mr_gauge_origin
  end if default_not_changed

  ! STEP 2 - Allocate the main work variables
  call push_current_memory()
  call initialize_MR()
  call initialize_integration()
  call reset_DFPT_timers()
  walltime_mat_conv = 0

  ! STEP 3 - Print general information
  call localorb_multi('', &
       & '#############################################', &
       & '#                                           #', &
       & '# WELCOME TO MAGNETIC RESPONSE CALCULATIONS #', &
       & '#                                           #', &
       & '#############################################', &
       & '', &
       & 'Current functionality includes:', &
       & '  * Indirect spin-spin couplings (J-couplings)', &
       & '  * NMR shieldings', &
       & '  * Magnetizabilities', &
       & '  * Electric Field Gradient (quadrupolar coupling)', &
       & '--------------------------------------------------', &
       & '', format='(2x, a)')
  call initial_general_info()

  ! STEP 4 - Do the calculations
  over_which_terms: do i_term = 1, size(which_terms)
     active_term: if (which_terms(i_term)) then
        call localorb_multi( &
             & 'Next term to be processed:', &
             & ' * '//term_names(i_term), &
             & '===============================', &
             & '', format='(2x, a)')
        call start_wall(walltime_all)
        call cpu_time(cputime_tmp1)
        if (term_names(i_term) == 'Spin-dipole') then
           call localorb_multi( &
                & 'Spin-dipole calculation is split into two parts:', &
                & 'First, the xx, yy, zz terms are computed', &
                & '', format='(2x, a)')
           calc_spin_dipole_diag = .true.
           call MR_core(term_names(i_term), timers(:,i_term))
           call localorb_multi('Next, the xy, xz, yz terms are computed', &
                & '', format='(2x, a)')
           calc_spin_dipole_diag = .false.
           call MR_core(term_names(i_term), timers(:,i_term))
        else
           call MR_core(term_names(i_term), timers(:,i_term))
        end if
        call stop_wall(walltime_all)
        call cpu_time(cputime_tmp2)
        clock_time_magnetic_response = &
             & clock_time_magnetic_response + (cputime_tmp2 - cputime_tmp1)
        ! Until this point, if use_scalapack==.true., each task held
        ! its own data. Here we synchronize the response tensor over
        ! tasks.
        if (use_scalapack) then
           select case(term_names(i_term))
           case('Fermi contact')
              call sync_vector(FC_tensor,          size(FC_tensor))
           case('Paramagnetic spin-orbit')
              call sync_vector(PSO_tensor,         size(PSO_tensor))
           case('Spin-dipole')
              call sync_vector(spin_dipole_tensor, size(spin_dipole_tensor))
           case('Diamagnetic spin-orbit')
              call sync_vector(DSO_tensor,         size(DSO_tensor))
           case('Paramagnetic shielding')
              call sync_vector(shield_para_tensor, size(shield_para_tensor))
           case('Diamagnetic shielding')
              call sync_vector(shield_dia_tensor,  size(shield_dia_tensor))
           case('Paramagnetic magnetizability')
              call sync_vector(mag_para_tensor,    size(mag_para_tensor))
           case('Diamagnetic magnetizability')
              call sync_vector(mag_dia_tensor,     size(mag_dia_tensor))
           case('Electric field gradient')
              call sync_vector(EFG_tensor,         size(EFG_tensor))
           end select
        end if
     end if active_term
  end do over_which_terms
  ! If J-couplings were requested, also compute the dipolar couplings
  ! because these cost nothing.
  if (j_is_active) call compute_dipolar_couplings(dipolar_tensor)

  ! STEP 5 - Convert spin-spin couplings to SI units using data for
  !          the most abundant NMR active isotopes (unless overridden
  !          by keywords in geometry.in).
  allocate(MR_nuclei(size(active_nuclei)), names_tmp(size(active_nuclei)), &
       & values_tmp((size(active_nuclei))))
  MR_nuclei = get_nuclear_data(species_name(species(active_nuclei)), &
       & default_isotope_number(species_name(species(active_nuclei)), 'NMR'))
  ! If any isotopes, magnetic moments, or spins were specified in
  ! geometry.in, override the default values.
  where (isotope_numbers(active_nuclei) > 0) &
       & MR_nuclei = get_nuclear_data(species_name(species(active_nuclei)), &
       & isotope_numbers(active_nuclei))
  names_tmp = MR_nuclei%name(4:)
  values_tmp = MR_nuclei%spin
  where (mr_custom_spins(active_nuclei) >= 0d0)
     names_tmp = '*'
     values_tmp = mr_custom_spins(active_nuclei)
  end where
  MR_nuclei%spin = values_tmp
  values_tmp = MR_nuclei%mu
  where (mr_custom_mus(active_nuclei) < 1d299)
     names_tmp = '*'
     values_tmp = mr_custom_mus(active_nuclei)
  end where
  do i_counter = 1, size(names_tmp)
     MR_nuclei(i_counter)%name(4:) = names_tmp(i_counter)
  end do
  MR_nuclei%mu = values_tmp
  if (j_is_active) then
     call convert_spin_spin_units(FC_tensor)
     call convert_spin_spin_units(PSO_tensor)
     call convert_spin_spin_units(spin_dipole_tensor)
     call convert_spin_spin_units(DSO_tensor)
     call convert_spin_spin_units(dipolar_tensor)
  end if

  ! STEP 6 - Print results
  call print_header()
  ! Start magnetic response JSON record
  json_log_handle = get_json_handle()
  if (out_aims_json_log) call fjson_start_MR_object(json_log_handle)
  ! Output shieldings
  if (shield_is_active) call print_shieldings(json_log_handle)
  ! Output dipolar couplings and J-couplings
  if (j_is_active) call print_spin_spin_couplings(MR_nuclei, json_log_handle)
  ! Output magnetizability
  if (magnet_is_active) call print_magnet(json_log_handle)
  ! Output EFG
  if (which_terms(9)) then
     ! The defaults nuclear properties are different for EFG and NMR
     ! calculations. If any isotopes, magnetic moments, or spins were
     ! specified in geometry.in, override the default values.
     where (isotope_numbers(active_nuclei) < 0 .and. &
          & mr_custom_spins(active_nuclei) < 0d0 .and. &
          & mr_custom_mus(active_nuclei) > 1d299) &
          & MR_nuclei = get_nuclear_data(species_name(species(active_nuclei)), &
          & default_isotope_number(species_name(species(active_nuclei)), 'EFG'))
     call print_efg(MR_nuclei)
  end if
  ! Output memory report
  call memory_report()
  ! Output timings
  call timings()
  ! End magnetic response JSON record
  if (myid == 0 .and. out_aims_json_log) &
       & call fjson_finish_object(json_log_handle)
  ! STEP 7 - Produce an sxml file
  if (output_sxml .and. myid == 0) call create_sxml_file(MR_nuclei)

  ! STEP 8 - clean up
  call cleanup_MR()
  call cleanup_integration()

  ! STEP 9 - Closing words
  call localorb_multi('', &
       & '############################################', &
       & '#                                          #', &
       & '# DONE WITH MAGNETIC RESPONSE CALCULATIONS #', &
       & '#                                          #', &
       & '############################################', &
       & '', format='(2x,a)')

contains
  ! Convert dipolar and J-couplings from atomic units to SI.
  subroutine convert_spin_spin_units(tensor)
    real(dp), intent(in out) :: tensor(:,:,:,:)
    real(dp), parameter :: atomic_time = 2.418884326509d-17 ! s
    real(dp) :: gamma1, gamma2
    integer :: j_atom
    do i_atom = 1, size(active_nuclei)
       gamma1 = MR_nuclei(i_atom)%mu/MR_nuclei(i_atom)%spin
       do j_atom = i_atom+1, size(active_nuclei)
          gamma2 = MR_nuclei(j_atom)%mu/MR_nuclei(j_atom)%spin
          tensor(:,:,i_atom,j_atom) = &
               & gamma1*gamma2/2/pi*tensor(:,:,i_atom,j_atom)/atomic_time
       end do
    end do
  end subroutine convert_spin_spin_units
end subroutine MR_main
