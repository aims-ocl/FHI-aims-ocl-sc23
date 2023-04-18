!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Subroutines for printing output before and after a magnetic
!!  response calculation.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
module MR_output

  use aims_memory_tracking, only: pop_memory_estimate
  use constants,            only: bohr
  use DFPT,                 only: timestamps_H1, timestamps_Ha1, &
       & timestamps_rho1, walltime_Ha1_update
  use dimensions,           only: n_basis, n_full_points, &
       & n_hamiltonian_matrix_size, n_spin
  use fjson_datatype,       only: fjson_handle
  use fjson_rw,             only: fjson_finish_array, fjson_finish_object, &
       & fjson_start_name_array, fjson_start_name_object, fjson_start_object, &
       & fjson_write_name_value
  use geometry,             only: coords
  use generate_aims_uuid,   only: write_aims_uuid
  use load_balancing,       only: batch_perm, n_bp_integ
  use localorb_io,          only: localorb_multi
  use MR_global
  use mpi_tasks,            only: aims_stop, mpi_min_max, myid, n_tasks
  use nuclear_data,         only: nucleus_t
  use pbc_lists,            only: column_index_hamiltonian
  use runtime_choices,      only: out_aims_json_log, sxml_filename, &
       & use_load_balancing
  use scalapack_wrapper,    only: eigenvec
  use timing,               only: time_magnetic_response
  use tools,                only: newunit, str, walltime_mat_conv
  use types,                only: dp

  implicit none

  private
  public :: create_sxml_file, fjson_start_MR_object, initial_general_info, &
       & memory_report, print_efg, print_header, print_shieldings, &
       & print_spin_spin_couplings, print_magnet, timings

  character(256) :: info_str

contains

  subroutine initial_general_info()
    integer :: i_count
    has_active_nuclei: if (j_is_active .or. shield_is_active .or. &
         & efg_is_active) then
       call localorb_multi('Computing responses for:', format='(2x, a)')
       do i_count = 1, size(active_nuclei)
          write(info_str, '(4x, a14, " at (", 2(f7.3, ", "), f7.3, ")")') &
               & trim(c_atoms(i_count)), coords(:,active_nuclei(i_count))*bohr
          call localorb_multi(info_str)
       end do
       call localorb_multi( &
            & '-------------------------------------------------------', &
            & '', format='(2x, a)')
    end if has_active_nuclei
    if (calc_full_tensor) then
       call localorb_multi('Both diagonal and off-diagonal tensor elements are &
            &calculated.', format='(2x, a)')
    else
       call localorb_multi('Only the diagonal tensor components are &
            &calculated.', format='(2x, a)')
    end if
    call localorb_multi('The following terms are included:', format='(2x, a)')
    do i_count = 1, size(which_terms)
       if (which_terms(i_count)) &
            & call localorb_multi('* '//term_names(i_count), format='(3x, a)')
    end do
    call localorb_multi('---------------------------------------------------', &
         & '', format='(2x, a)')
  end subroutine initial_general_info

  subroutine print_header()
    call localorb_multi( &
         & 'Results for the magnetic response', &
         & '=================================', &
         & '', &
         & format='(2x, a)')
    if (no_giao) then
       write(info_str, '("(", 2(f7.3, ", "), f7.3, ")")') gauge_origin*bohr
       call localorb_multi( &
            & 'Not using GIAOs! Gauge origin is at '//trim(info_str), &
            & '---------------------------------------------------------------&
            &', '', format='(2x, a)')
    end if
  end subroutine print_header

  subroutine fjson_start_MR_object(json_log_handle)
    type(fjson_handle), intent(in out) :: json_log_handle
    if (myid /= 0) return
    call fjson_start_object(json_log_handle)
    call fjson_write_name_value(json_log_handle, 'record_type', &
         & 'magnetic response')
    call fjson_write_name_value(json_log_handle, 'using GIAOs', .not. no_giao)
    call fjson_write_name_value(json_log_handle, 'calculate full tensor', &
         & calc_full_tensor)
  end subroutine fjson_start_MR_object

  subroutine print_shieldings(json_log_handle)
    type(fjson_handle), intent(in out) :: json_log_handle
    real(dp) :: tmp_tensor(3,3), eigenvalues(3), workspace(9)
    integer :: i_atom, info
    call localorb_multi(&
         & 'NMR shielding tensors (ppm)', &
         & '---------------------------', '', format='(2x, a)')
    call localorb_multi( &
         & 'Definitions:', &
         & '  Anisotropy: sigma_z - (sigma_y + sigma_x)/2', &
         & '  Asymmetry:  (sigma_y - sigma_x)/(sigma_z - sigma_iso)', &
         & '  Span:       sigma_z - sigma_x', &
         & '  Skew:       3(sigma_y - sigma_iso)/(sigma_z - sigma_x)', &
         & '  sigma_x < sigma_y < sigma_z are the eigenvalues of the shielding &
         &tensor', &
         & '', format='(4x, a)')
    if (myid == 0 .and. out_aims_json_log) &
         & call fjson_start_name_array(json_log_handle, 'shielding tensors')
    do i_atom = 1, size(active_nuclei)
       call localorb_multi(trim(c_atoms(i_atom))//':', format='(3x, a)')
       tmp_tensor = 1d6*shield_dia_tensor(:,:,i_atom)
       write(info_str, '(es13.6e2)') trace(tmp_tensor)/3
       call localorb_multi('Diamagnetic:  '//trim(info_str)//' (isotropic)', &
            & format='(4x, a)')
       call print_tensor(tmp_tensor)
       tmp_tensor = 1d6*shield_para_tensor(:,:,i_atom)
       write(info_str, '(es13.6e2)') trace(tmp_tensor)/3
       call localorb_multi('Paramagnetic: '//trim(info_str)//' (isotropic)', &
            & format='(4x, a)')
       call print_tensor(tmp_tensor)
       tmp_tensor = 1d6*(shield_para_tensor(:,:,i_atom) + &
            & shield_dia_tensor(:,:,i_atom))
       write(info_str, '(es13.6e2)') trace(tmp_tensor)/3
       call localorb_multi('Total:        '//trim(info_str)// &
            & ' (isotropic)', format='(4x, a)')
       if (.not. calc_full_tensor) then
          ! Set all off-diagonal elements to zero to avoid ambiguity.
          tmp_tensor(1,2:) = 0d0
          tmp_tensor(2,1::2) = 0d0
          tmp_tensor(3,:2) = 0d0
       end if
       call print_tensor(tmp_tensor)
       json_shield: if (myid == 0 .and. out_aims_json_log) then
          call fjson_start_object(json_log_handle)
          call fjson_write_name_value(json_log_handle, 'atom', &
               & active_nuclei(i_atom))
          call json_print_tensor(json_log_handle, tmp_tensor, 'ppm')
          call fjson_write_name_value(json_log_handle, '1/3 trace', &
               & trace(tmp_tensor)/3)
       end if json_shield
       ! Find eigenvalues of the shielding tensor
       call dsyev('n', 'u', 3, tmp_tensor, 3, eigenvalues, workspace, 9, info)
       ! Reorder them so that sigma_x < sigma_y < sigma_z
       call heapsort(eigenvalues, size(eigenvalues))
       ! Print anisotropy, asymmetry, span, skew
       write(info_str, '(es13.6e2)') eigenvalues(3) - sum(eigenvalues(1:2))/2
       call localorb_multi('Anisotropy:   '//trim(info_str), &
            & format='(4x, a)')
       write(info_str, '(es13.6e2)') (eigenvalues(2) - eigenvalues(1))/ &
            & (eigenvalues(3) - sum(eigenvalues)/3)
       call localorb_multi('Asymmetry:    '//trim(info_str), &
            & format='(4x, a)')
       write(info_str, '(es13.6e2)') eigenvalues(3) - eigenvalues(1)
       call localorb_multi('Span:         '//trim(info_str), &
            & format='(4x, a)')
       write(info_str, '(es13.6e2)') 3*(eigenvalues(2) - sum(eigenvalues)/3)/ &
            & (eigenvalues(3) - eigenvalues(1))
       call localorb_multi('Skew:         '//trim(info_str), &
            & format='(4x, a)')
       call localorb_multi('')
       json_other: if (myid == 0 .and. out_aims_json_log) then
          call fjson_write_name_value(json_log_handle, 'anisotropy', &
               & eigenvalues(3) - sum(eigenvalues(1:2))/2)
          call fjson_write_name_value(json_log_handle, 'asymmetry', &
               & (eigenvalues(2) - eigenvalues(1))/ &
               & (eigenvalues(3) - sum(eigenvalues)/3))
          call fjson_write_name_value(json_log_handle, 'span', &
               & eigenvalues(3) - eigenvalues(1))
          call fjson_write_name_value(json_log_handle, 'skew', &
               & 3*(eigenvalues(2) - sum(eigenvalues)/3)/ &
               & (eigenvalues(3) - eigenvalues(1)))
          call fjson_finish_object(json_log_handle)
       end if json_other
    end do
    if (myid == 0 .and. out_aims_json_log) &
         & call fjson_finish_array(json_log_handle)
  end subroutine print_shieldings

  subroutine print_spin_spin_couplings(MR_nuclei, json_log_handle)
    type(nucleus_t), intent(in) :: MR_nuclei(:)
    type(fjson_handle), intent(in out) :: json_log_handle
    real(dp) :: tmp_tensor(3,3,size(active_nuclei),size(active_nuclei))
    real(dp) :: eigenvalues(3), workspace(9)
    integer :: i_term, i_atom, j_atom
    call localorb_multi(&
         & 'J-couplings (Hz)', &
         & '----------------', '', format='(2x, a)')
    do i_term = 1, 4
       if (.not. which_terms(i_term)) cycle
       select case(i_term)
       case(1); tmp_tensor = FC_tensor
       case(2); tmp_tensor = PSO_tensor
       case(3); tmp_tensor = spin_dipole_tensor
       case(4); tmp_tensor = DSO_tensor
       end select
       call localorb_multi(trim(term_names(i_term))//':', &
            & format='(3x, a)')
       if (calc_full_tensor .and. i_term == 3) & ! Spin-dipole comment
            & call localorb_multi('Off-diagonal elements not evaluated!', &
            & format='(4x, a)')
       do i_atom = 1, size(active_nuclei)
          do j_atom = i_atom+1, size(active_nuclei)
             call localorb_multi(trim(c_atoms(i_atom))//' - '// &
                  & trim(c_atoms(j_atom))//':', format='(4x, a)')
             call print_tensor(tmp_tensor(:,:,i_atom,j_atom))
          end do
       end do
       call localorb_multi('Isotropic:', format='(4x, a)')
       call print_couplings()
    end do
    call localorb_multi('Total:', format='(3x, a)')
    tmp_tensor = FC_tensor + PSO_tensor + DSO_tensor + spin_dipole_tensor
    call print_couplings()
    json_J: if (myid == 0 .and. out_aims_json_log) then
       call fjson_start_name_array(json_log_handle, 'J-couplings')
       do i_atom = 1, size(active_nuclei)
          do j_atom = i_atom+1, size(active_nuclei)
             call fjson_start_object(json_log_handle)
             call fjson_write_name_value(json_log_handle, 'atom 1', &
                  & active_nuclei(i_atom))
             call fjson_write_name_value(json_log_handle, 'isotope 1', &
                  & trim(adjustl(MR_nuclei(i_atom)%name)))
             call fjson_write_name_value(json_log_handle, 'atom 2', &
                  & active_nuclei(j_atom))
             call fjson_write_name_value(json_log_handle, 'isotope 1', &
                  & trim(adjustl(MR_nuclei(j_atom)%name)))
             call json_print_tensor(json_log_handle, &
                  & tmp_tensor(:,:,i_atom,j_atom), 'Hz')
             call fjson_write_name_value(json_log_handle, '1/3 trace', &
                  & trace(tmp_tensor(:,:,i_atom,j_atom))/3)
             call fjson_finish_object(json_log_handle)
          end do
       end do
       call fjson_finish_array(json_log_handle)
    end if json_J
    call localorb_multi(&
         & 'Dipolar coupling tensors (kHz)', &
         & '------------------------------', &
         & '', format='(2x, a)')
    dipolar_tensor = 1d-3*dipolar_tensor ! Output in kHz!
    tmp_tensor = 0d0
    do i_atom = 1, size(active_nuclei)
       do j_atom = i_atom+1, size(active_nuclei)
          call localorb_multi(trim(c_atoms(i_atom))//' - '// &
               & trim(c_atoms(j_atom))//':', format='(4x, a)')
          call print_tensor(dipolar_tensor(:,:,i_atom,j_atom))
          ! Find the eigenvalues of the dipolar tensor
          call dsyev('n', 'u', 3, dipolar_tensor(:,:,i_atom,j_atom), 3, &
               & eigenvalues, workspace, 9, i_term)
          ! The largest eigenvalue is the experimentally measured one
          tmp_tensor(1,1,i_atom,j_atom) = 3*maxval(eigenvalues)
       end do
    end do
    call print_couplings()
    json_dipolar: if (myid == 0 .and. out_aims_json_log) then
       call fjson_start_name_array(json_log_handle, 'dipolar couplings')
       do i_atom = 1, size(active_nuclei)
          do j_atom = i_atom+1, size(active_nuclei)
             call fjson_start_object(json_log_handle)
             call fjson_write_name_value(json_log_handle, 'atom 1', &
                  & active_nuclei(i_atom))
             call fjson_write_name_value(json_log_handle, 'atom 2', &
                  & active_nuclei(j_atom))
             call fjson_write_name_value(json_log_handle, 'value in kHz', &
                  & tmp_tensor(1,1,i_atom,j_atom)/3)
             call fjson_finish_object(json_log_handle)
          end do
       end do
       call fjson_finish_array(json_log_handle)
    end if json_dipolar
    call localorb_multi('Isotopes used for dipolar and J-couplings:', &
         & '', format='(2x, a)')
    do i_atom = 1, size(MR_nuclei)
       write(info_str, '(i5, ": ", a, " at (", 2(f7.3, ", "), &
            &f7.3, ")")') active_nuclei(i_atom), MR_nuclei(i_atom)%name, &
            & coords(:,active_nuclei(i_atom))*bohr
       call localorb_multi(info_str)
    end do
    call localorb_multi('')
  contains
    subroutine print_couplings()
      real(dp) :: J_value
      integer :: str_len, j_atom, i_atom
      character(:), allocatable :: info_str_tmp
      write(info_str, '(5x, '//str(size(active_nuclei)-1)//'(i10, 5x))') &
           & [(active_nuclei(j_atom), j_atom = 2, size(active_nuclei))]
      call localorb_multi(info_str)
      do i_atom = 1, size(active_nuclei)-1
         write(info_str, '(i5)') active_nuclei(i_atom)
         str_len = len(trim(info_str)) - 15
         over_l: do j_atom = 2, size(active_nuclei)
            str_len = str_len + 15
            if (j_atom > i_atom) then
               J_value = trace(tmp_tensor(:,:,i_atom,j_atom))/3
               ! For printing purposes put all values lower than 1d-99
               ! to zero and greater than 1d99 to infinity.
               if (abs(J_value) > 1d+99) J_value = J_value/0d0
               if (abs(J_value) < 1d-99) J_value = 0d0
               info_str_tmp = info_str(:str_len)
               write(info_str, '(a, es15.6e2)') info_str_tmp, J_value
            end if
         end do over_l
         write(info_str, '(a, i5)') trim(info_str)//'  ', active_nuclei(i_atom)
         call localorb_multi(info_str)
      end do
      call localorb_multi('')
    end subroutine print_couplings
  end subroutine print_spin_spin_couplings

  subroutine print_efg(MR_nuclei)
    type(nucleus_t), intent(in) :: MR_nuclei(:)
    real(dp) :: efg(3,3)
    real(dp) :: efg_principle_values(3)
    real(dp) :: workspace(9), efg_cq, efg_eta, efg_theta, efg_phi
    integer :: i_atom, info
    call localorb_multi(&
         & 'Electric field gradient tensor in a.u. and general axis system', &
         & '--------------------------------------------------------------', &
         & '', format='(2x, a)')
    do i_atom = 1, size(active_nuclei)
       call localorb_multi(trim(c_atoms(i_atom))//':', format='(3x, a)')
       call localorb_multi('', '     Total EFG tensor in a.u.', '')
       efg = EFG_tensor(:,:,i_atom)
       call print_tensor(efg)
       call localorb_multi('')
       efg = (efg+transpose(efg))/2d0
       ! Diagonalize matrix using LAPACK routine for a real,
       ! symmetric matrix - including the calculation of eigenvectors
       ! etc
       call dsyev('v','u', 3, efg, 3, efg_principle_values, workspace, 9, info)
       if (info /= 0) &
            & call aims_stop('ERROR in DSYEV', 'MagneticResponse/MR_main')
       ! Sort eigenvalues and -vectors according the absolute
       ! eigenvalues in ascending order
       call sort_eigenvalues(efg_principle_values, efg)
       call localorb_multi(&
            & '     EFG eigenvalues and -vectors in a.u.', '', format='(2x, a)')
       call print_eigenv(efg_principle_values, efg)
       call evaluate_quadrupolar_coupling(efg_principle_values, efg, &
            & efg_cq, efg_eta, efg_theta, efg_phi, MR_nuclei(i_atom)%q_in_barn)

       write(info_str, '(es13.4e2)') MR_nuclei(i_atom)%q_in_barn
       call localorb_multi('     nuclear Q /Barn:         '// &
            & trim(info_str), format='(4x, a)')
       write(info_str, '(es13.4e2)') efg_cq
       call localorb_multi('           C_Q /MHz :         '// &
            & trim(info_str), format='(4x, a)')
       if (abs(efg_eta + 1d0) < epsilon(1d0)) then
          call localorb_multi('           eta      :         N/A', &
               & format='(4x, a)')
       else
          write(info_str, '(es13.4e2)') efg_eta
          call localorb_multi('           eta      :         '// &
               & trim(info_str), format='(4x, a)')
       endif
       write(info_str, '(es13.4e2)') efg_theta
       call localorb_multi('     theta /rad (polar)  :         '// &
            & trim(info_str), format='(4x, a)')
       write(info_str, '(es13.4e2)') efg_phi
       call localorb_multi('     phi /rad (azimuthal):         '// &
            & trim(info_str), '', format='(4x, a)')
    end do
    call localorb_multi('Isotopes used for EFG:', '', format='(2x, a)')
    do i_atom = 1, size(MR_nuclei)
       write(info_str, '(i5, ": ", a, " at (", 2(f7.3, ", "), &
            &f7.3, ")")') active_nuclei(i_atom), MR_nuclei(i_atom)%name, &
            & coords(:,active_nuclei(i_atom))*bohr
       call localorb_multi(info_str)
    end do
    call localorb_multi('')
  end subroutine print_efg

  subroutine evaluate_quadrupolar_coupling(evals, evecs, cq, eta, &
       & theta, phi, nucl_q)
    real(dp) :: evals(3), evecs(3,3), cq, eta, theta, phi, nucl_q
    ! conversion of EFG from a.u. (Hartree / (e * a_0^2)) to SI
    !         e * V_zz * Q
    ! C_Q = ------------------
    !            h
    !
    !               Ha      4.35974*10^-18      J
    ! [e * V_zz] = ------ = ------------------  ---
    !               a_0^2   (5.29177*10^-11)^2  m^2
    !
    !       V_zz * Q[m^2 to Barn]
    !C_Q = ---------------------- * [Hz to MHz]
    !       h
    !
    !            V_zz * 10^-28 * 10^-6    J    *  m^2     *      MHz
    !        = ------------------------   ---------------------------
    !            6.626070*10^-34          m^2  * Barn  *  J*s  *  Hz
    !
    !                Ha                        a_0^2 * MHz
    !        = V_zz ----- * 2.3496471626*10e2  ----------- * Q Barn
    !               a_0^2                       Ha   * Barn
    cq = evals(3)*2.3496471626d2*nucl_q
    if (evals(3) < 0d0) then
       eta = -1d0
    else
       eta = (evals(1) - evals(2)) / evals(3)
    end if
    ! polar angle of eigenvector zz
    !           V_zz(z)
    ! acos( -------------- )
    !           |V_zz|
    !|V_zz| of normalized eigenvector = 1
    theta = acos(evecs(3,3))
    ! azimuthal angle of eigenvector zz
    !           V_zz(y)
    ! atan2( ------------- )
    !           V_zz(x)
    phi   = atan2(evecs(2,3), evecs(1,3))
  end subroutine evaluate_quadrupolar_coupling

  subroutine print_eigenv(evals, evecs)
    real(dp) :: evals(3), evecs(3,3)
    write(info_str, '(es13.4e2)') evals(1)
    call localorb_multi('     eigenvalue xx:         '// &
         & trim(info_str), format='(4x, a)')
    write(info_str, '(3(es15.4e2))') evecs(:,1)
    call localorb_multi('     eigenvector xx:        '// &
         & trim(info_str), format='(4x, a)')
    call localorb_multi('')
    write(info_str, '(es13.4e2)') evals(2)
    call localorb_multi('     eigenvalue yy:         '// &
         & trim(info_str), format='(4x, a)')
    write(info_str, '(3(es15.4e2))') evecs(:,2)
    call localorb_multi('     eigenvector yy:        '// &
         & trim(info_str), format='(4x, a)')
    call localorb_multi('')
    write(info_str, '(es13.4e2)') evals(3)
    call localorb_multi('     eigenvalue zz:         '// &
         & trim(info_str), format='(4x, a)')
    write(info_str, '(3(es15.4e2))') evecs(:,3)
    call localorb_multi('     eigenvector zz:        '// &
         & trim(info_str), format='(4x, a)')
    call localorb_multi('')
  end subroutine print_eigenv

  pure real(dp) function trace(x) result(y)
    real(dp), intent(in) :: x(:,:)
    integer :: i_counter
    y = sum([(x(i_counter,i_counter), i_counter=1, size(x,1))])
  end function trace

  subroutine print_tensor(T)
    real(dp), intent(in) :: T(3,3)
    ! Fortran note: VALUE attribute would conflict with DIMENSION attribute
    real(dp) :: T_loc(3,3)
    integer :: i_dir
    T_loc = T
    ! For printing purposes put all values lower than 1d-99 to zero
    ! and greater than 1d99 to infinity.
    where (abs(T_loc) > 1d+99) T_loc = T_loc/0d0
    where (abs(T_loc) < 1d-99) T_loc = 0d0
    do i_dir = 1, 3
       write(info_str, '(6x, 3(es15.6e2))') T_loc(i_dir,:)
       call localorb_multi(info_str)
    end do
  end subroutine print_tensor

  subroutine json_print_tensor(json_log_handle, T, unit)
    type(fjson_handle), intent(in out) :: json_log_handle
    real(dp), intent(in) :: T(3,3)
    character(*), intent(in) :: unit
    call fjson_write_name_value(json_log_handle, 'xx in '//unit, T(1,1))
    call fjson_write_name_value(json_log_handle, 'yy in '//unit, T(2,2))
    call fjson_write_name_value(json_log_handle, 'zz in '//unit, T(3,3))
    if (calc_full_tensor) then
       call fjson_write_name_value(json_log_handle, 'xy in '//unit, T(1,2))
       call fjson_write_name_value(json_log_handle, 'xz in '//unit, T(1,3))
       call fjson_write_name_value(json_log_handle, 'yx in '//unit, T(2,1))
       call fjson_write_name_value(json_log_handle, 'yz in '//unit, T(2,3))
       call fjson_write_name_value(json_log_handle, 'zx in '//unit, T(3,1))
       call fjson_write_name_value(json_log_handle, 'zy in '//unit, T(3,2))
    end if
  end subroutine json_print_tensor

  subroutine print_magnet(json_log_handle)
    type(fjson_handle), intent(in out) :: json_log_handle
    real(dp) :: tmp_tensor(3,3)
    ! Conversion factor for magnetizability:
    ! (e^2*hbar^2)/(m^3*c^2*alpha^2)*1d30
    associate(xi_au_to_SI => 78.910355719306921d0)
      ! Minus sign from the definition of chi
      mag_para_tensor = -xi_au_to_SI * mag_para_tensor
      mag_dia_tensor  = -xi_au_to_SI * mag_dia_tensor
    end associate
    call localorb_multi( &
         & 'Magnetizability (10^-30 J/T^2)', &
         & '------------------------------', &
         & '', format='(2x, a)')
    write(info_str, '(es13.6e2)') trace(mag_dia_tensor)/3
    call localorb_multi('Diamagnetic:  '//trim(info_str)//' (isotropic)', &
         & format='(4x, a)')
    call print_tensor(mag_dia_tensor)
    write(info_str, '(es13.6e2)') trace(mag_para_tensor)/3
    call localorb_multi('Paramagnetic: '//trim(info_str)//' (isotropic)', &
         & format='(4x, a)')
    call print_tensor(mag_para_tensor)
    write(info_str, '(es13.6e2)') trace(mag_para_tensor+mag_dia_tensor)/3
    call localorb_multi('Total:        '//trim(info_str)// &
         & ' (isotropic)', format='(4x, a)')
    tmp_tensor = mag_para_tensor + mag_dia_tensor
    if (.not. calc_full_tensor) then
       ! Set all off-diagonal elements to zero to avoid ambiguity.
       tmp_tensor(1,2:) = 0d0
       tmp_tensor(2,1::2) = 0d0
       tmp_tensor(3,:2) = 0d0
    end if
    call print_tensor(tmp_tensor)
    call localorb_multi('')
    json_magnet: if (myid == 0 .and. out_aims_json_log) then
       call fjson_start_name_object(json_log_handle, 'magnetizability')
       call json_print_tensor(json_log_handle, tmp_tensor, '10^-30 J/T^2')
       call fjson_write_name_value(json_log_handle, '1/3 trace', &
            & trace(tmp_tensor)/3)
       call fjson_finish_object(json_log_handle)
    end if json_magnet
  end subroutine print_magnet

  subroutine memory_report()
    real(dp) :: r_min, r_max
    integer :: i_minloc, i_maxloc, n_digits
    n_digits = int(log(real(n_tasks,8))/log(10d0),4)+1
    call localorb_multi( &
         & 'Memory report (for developers)', &
         & '==============================', &
         & '', format='(2x, a)')
    call localorb_multi( &
         & '  Maximum memory usage of magnetic response calculations', &
         & '                     min (task)          max (task)', &
         & format='(2x, a)')
    max_total_memory = pop_memory_estimate()
    call mpi_min_max(real(max_total_memory,8), r_min, i_minloc, r_max, i_maxloc)
    write(info_str, '(16x, 2(f11.2, " MB (",i'//str(n_digits)//',")"))') &
         & r_min/1d6, i_minloc, r_max/1d6, i_maxloc
    call localorb_multi(info_str)
    call localorb_multi('', format='(4x, 66("-"), a)')
    call localorb_multi( &
         & '  Maximum memory usage for the DFPT part separately', &
         & '                     min (task)          max (task)', &
         & format='(2x, a)')
    call mpi_min_max(real(max_dfpt_memory,8), r_min, i_minloc, r_max, i_maxloc)
    write(info_str, '(16x, 2(f11.2, " MB (",i'//str(n_digits)//',")"))') &
         & r_min/1d6, i_minloc, r_max/1d6, i_maxloc
    call localorb_multi(info_str)
    call localorb_multi('', format='(4x, 66("-"), a)')
    if (use_load_balancing) then
       call localorb_multi( &
            & '  batch_perm(n_bp_integ)%n_local_matrix_size -', &
            & '    characteristic size of packed matrices', &
            & '                     min (task)          max (task)', &
            & format='(2x, a)')
       call mpi_min_max(real(batch_perm(n_bp_integ)%n_local_matrix_size,8), &
            & r_min, i_minloc, r_max, i_maxloc)
       write(info_str, '(16x, 2(i14, " (",i'//str(n_digits)//',")"))') &
            & nint(r_min,8), i_minloc, nint(r_max,8), i_maxloc
       call localorb_multi(info_str)
       write(info_str, '(4x, a, 2(f11.2, " MB (",i'//str(n_digits)//',")"))') &
            & 'For real(8):', 8d-6*r_min, i_minloc, 8d-6*r_max, i_maxloc
       call localorb_multi(info_str)
    else
       call localorb_multi( &
            & '  n_hamiltonian_matrix_size - &
            &characteristic size of packed matrices', &
            & '                     min (task)          max (task)', &
            & format='(2x, a)')
       call mpi_min_max(real(n_hamiltonian_matrix_size,8), r_min, i_minloc, &
            & r_max, i_maxloc)
       write(info_str, '(16x, 2(i14, " (",i'//str(n_digits)//',")"))') &
            & nint(r_min,8), i_minloc, nint(r_max,8), i_maxloc
       call localorb_multi(info_str)
       write(info_str, '(4x, a, 2(f11.2, " MB (",i'//str(n_digits)//',")"))') &
            & 'For real(8):', 8d-6*r_min, i_minloc, 8d-6*r_max, i_maxloc
       call localorb_multi(info_str)
    end if
    call localorb_multi('', format='(4x, 66("-"), a)')
    call localorb_multi('  size(eigenvec) - &
         &characteristic size of 2D block cyclic matrices', &
         & '                      min (task)          max (task)', &
         & format='(2x, a)')
    call mpi_min_max(real(size(eigenvec),8), r_min, i_minloc, r_max, i_maxloc)
    write(info_str, '(16x, 2(i14, " (",i'//str(n_digits)//',")"))') &
         & nint(r_min,8), i_minloc, nint(r_max,8), i_maxloc
    call localorb_multi(info_str)
    write(info_str, '(4x, a, 2(f11.2, " MB (",i'//str(n_digits)//',")"))') &
         & 'For real(8):', 8d-6*r_min, i_minloc, 8d-6*r_max, i_maxloc
    call localorb_multi(info_str)
    call localorb_multi('', format='(4x, 66("-"), a)')
    call localorb_multi('  n_full_points - &
         &characteristic size of grid point based arrays', &
         & '                      min (task)          max (task)', &
         & format='(2x, a)')
    call mpi_min_max(real(n_full_points,8), r_min, i_minloc, r_max, i_maxloc)
    write(info_str, '(16x, 2(i14, " (",i'//str(n_digits)//',")"))') &
         & nint(r_min,8), i_minloc, nint(r_max,8), i_maxloc
    call localorb_multi(info_str)
    write(info_str, '(4x, a, 2(f11.2, " MB (",i'//str(n_digits)//',")"))') &
         & 'For real(8):', 8d-6*r_min, i_minloc, 8d-6*r_max, i_maxloc
    call localorb_multi(info_str)
    call localorb_multi('', format='(4x, 66("-"), a)')
    call localorb_multi('Global dimensions (size in MB shown for real(8))', &
         & format='(4x, a)')
    write(info_str, '(6x, a, 3(i0, a), f10.2, a)') &
         & '(n_basis x n_basis) = (', n_basis, ' x ', n_basis, ') = ', &
         & int(n_basis,8)*int(n_basis,8), ' = ', &
         & 8d-6*real(n_basis,8)*real(n_basis,8), ' MB'
    call localorb_multi(info_str)
    call localorb_multi('', format='(4x, 66("-"), a)')
  end subroutine memory_report

  subroutine timings()
    character(:), allocatable :: fmt
    integer :: wall_maxloc, count_rate
    real(dp) :: r_min, r_max, wall_tmp, wall_max
    real(dp) :: walltime_sum ! Sum of individual terms
    integer :: i_minloc, i_maxloc, i_term, n_digits
    n_digits = int(log(real(n_tasks,8))/log(10d0),4)+1
    call system_clock(count_rate=count_rate)
    fmt = '(4x, a, es11.4e2, " s (",i'//str(n_digits)// &
         & ',") | ", 3(es11.4e2, " s (",i'//str(n_digits)//',")"))'
    call localorb_multi('', &
         & 'Timings (for developers)', &
         & '========================', &
         & format='(2x, a)')
    write(info_str, '(a, 4x, '//str(n_digits)//'x, "|", 14x, a)') &
         & 'wall time', 'cpu time'
    call localorb_multi(info_str, format='(36x, a)')
    call localorb_multi('', format='(33x, '//str(50+3*n_digits)//'("-"), a)')
    write(info_str, '("max (task)", '//str(n_digits)// &
         & 'x, "    |    min (task)", '//str(n_digits)// &
         & 'x, 6x, "max (task)", a)')
    call localorb_multi(info_str, format='(35x, a)')
    ! Non-self-consistent parts of the first-order Hamiltonian (H1)
    walltime_sum = 0
    do i_term = 1, size(timers,2)
       if (which_terms(i_term)) then
          call mpi_min_max(timers(3,i_term), r_max=wall_max, &
               & i_maxloc=wall_maxloc)
          call mpi_min_max(timers(4,i_term), r_min, i_minloc, r_max, i_maxloc)
          write(info_str, fmt) term_names(i_term), wall_max, wall_maxloc, &
               & r_min, i_minloc, r_max, i_maxloc
          call localorb_multi(info_str)
          walltime_sum = walltime_sum + wall_max
       end if
    end do
    ! Self-consistent part of the first-order Hamiltonian (H1)
    call mpi_min_max(timestamps_H1(4), r_max=wall_max, i_maxloc=wall_maxloc)
    call mpi_min_max(timestamps_H1(3), r_min, i_minloc, r_max, i_maxloc)
    write(info_str, fmt) 'First-order xc              ', &
         & wall_max, wall_maxloc, r_min, i_minloc, r_max, i_maxloc
    call localorb_multi(info_str)
    walltime_sum = walltime_sum + wall_max
    ! rho1, grad rho1 update
    call mpi_min_max(timestamps_rho1(4), r_max=wall_max, i_maxloc=wall_maxloc)
    call mpi_min_max(timestamps_rho1(3), r_min, i_minloc, r_max, i_maxloc)
    write(info_str, fmt) 'First-order density         ', &
         & wall_max, wall_maxloc, r_min, i_minloc, r_max, i_maxloc
    call localorb_multi(info_str)
    walltime_sum = walltime_sum + wall_max
    ! First-order Hartree potential
    if (n_spin == 2) then
       call mpi_min_max(timestamps_Ha1(4), r_max=wall_max, i_maxloc=wall_maxloc)
       call mpi_min_max(timestamps_Ha1(3), r_min, i_minloc, r_max, i_maxloc)
       write(info_str, fmt) 'First-order Hartree         ', &
            & wall_max, wall_maxloc, r_min, i_minloc, r_max, i_maxloc
       call localorb_multi(info_str)
       walltime_sum = walltime_sum + wall_max
       wall_tmp = real(walltime_Ha1_update,8)/count_rate
       call mpi_min_max(wall_tmp, r_max=wall_max, i_maxloc=wall_maxloc)
       write(info_str, fmt) 'First-order Ha update       ', wall_max, &
            & wall_maxloc
       call localorb_multi(info_str)
       walltime_sum = walltime_sum + wall_max
    end if
    ! Matrix matrix multiplications
    wall_tmp = real(walltime_mul,8)/count_rate
    call mpi_min_max(wall_tmp, r_max=wall_max, i_maxloc=wall_maxloc)
    write(info_str, fmt) 'Mat-mat multiplications     ', wall_max, wall_maxloc
    call localorb_multi(info_str)
    walltime_sum = walltime_sum + wall_max
    ! packed <-> block cyclic
    wall_tmp = real(walltime_mat_conv,8)/count_rate
    call mpi_min_max(wall_tmp, r_max=wall_max, i_maxloc=wall_maxloc)
    write(info_str, fmt) 'packed <-> block cyclic     ', wall_max, wall_maxloc
    call localorb_multi(info_str)
    walltime_sum = walltime_sum + wall_max
    ! Calls to MR_core
    call localorb_multi('', format='(4x, '//str(79+3*n_digits)//'("-"), a)')
    time_magnetic_response = real(walltime_all,8)/count_rate
    call mpi_min_max(time_magnetic_response, r_max=wall_max, &
         & i_maxloc=wall_maxloc)
    write(info_str, fmt) 'Total                       ', wall_max, wall_maxloc
    call localorb_multi(info_str)
    ! Difference between total and individual wall times.
    write(info_str, '(4x, a, es11.4e2, " s")') &
         & 'Total minus individual terms', wall_max - walltime_sum
    call localorb_multi(info_str)
  end subroutine timings

  !!  FUNCTION
  !!
  !!  Sorts the eigenvalues in ascending order corresponding to their
  !!  ABSOLUTE value
  !!  |evals(3)| > |evals(2)| > |evals(1)|
  !!  and sorts the eigenvectors accordingly.
  !!  INPUT:    evals = real(3) array, real (& unsorted) eigenvalues
  !!            evecs = real(3,3) array, real (& unsorted) eigenvectors
  !!  OUTPUT:   evals = real(3) array, real & sorted eigenvalues
  !!            evecs = real(3,3) array, real & sorted eigenvectors
  !!
  subroutine sort_eigenvalues(evals, evecs)
    real(dp), intent(in out) :: evals(3), evecs(3,3)
    real(dp) :: e_combined(5,3)
    e_combined(1,:)   = abs(evals)
    e_combined(2,:)   = evals
    e_combined(3:5,:) = evecs
    call heapsort_general(e_combined, 5, 3, 1)
    evecs = e_combined(3:5,:)
    evals = e_combined(2,:)
  end subroutine sort_eigenvalues

  !!  FUNCTION
  !!
  !!  Creates an sxml file. The structure of the sxml file is
  !!  specified in Biternas et al., J. Mag. Res. 240, 124 (2014).
  !!
  subroutine create_sxml_file(MR_nuclei)
    type(nucleus_t), intent(in) :: MR_nuclei(:)
    character(:), allocatable :: name
    character(10) :: coords_str(3)
    ! See the output sxml file to understand id.
    integer :: id
    character(256) :: info_str
    integer :: i_atom, j_atom, io_stat, sxml_unit
    real(dp), allocatable :: shieldings(:,:,:), J_couplings(:,:,:,:)
    allocate(J_couplings(3,3,size(active_nuclei),size(active_nuclei)))
    allocate(shieldings(3,3,size(active_nuclei)))
    shieldings = 1d6*(shield_para_tensor + shield_dia_tensor)
    J_couplings = FC_tensor + PSO_tensor + DSO_tensor + spin_dipole_tensor
    open(newunit(sxml_unit), file=trim(sxml_filename), &
         & action='write', form='formatted', iostat=io_stat, iomsg=info_str)
    if (io_stat /= 0) call aims_stop(trim(info_str), &
         & 'MagneticResponse/MR_main::create_sxml_file')
    write(sxml_unit, '(a)') '<?xml version="1.0" encoding="utf-8"?>'
    write(sxml_unit, '(a)') '<spin_system>'
    ! Spin system
    do i_atom = 1, size(active_nuclei)
       name = MR_nuclei(i_atom)%name
       name = trim(name(4:))//adjustl(name(:2))
       info_str = '<spin id="'//str(i_atom)//'" isotope="'//trim(name)//'" >'
       write(sxml_unit, '(2x, a)') trim(info_str)
       write(coords_str, '(f10.3)') coords(:,active_nuclei(i_atom))*bohr
       coords_str = adjustl(coords_str)
       write(sxml_unit, '(4x, a)') '<coordinates x="'//trim(coords_str(1))//&
            &'" y="'//trim(coords_str(1))//'" z="'//trim(coords_str(1))//'" />'
       write(sxml_unit, '(2x, a)') '</spin>'
    end do
    id = 1
    ! Shieldings
    if (shield_is_active) then
       do i_atom = 1, size(active_nuclei)
          info_str = '<interaction kind="shielding" id="'//str(id)// &
               & '" units="ppm" spin_a="'//str(i_atom)//'" reference="">'
          write(sxml_unit, '(2x, a)') trim(info_str)
          info_str = '<tensor xx="'//str(shieldings(1,1,i_atom))//'" xy="'// &
               & str(shieldings(1,2,i_atom))//'" xz="'// &
               & str(shieldings(1,3,i_atom))//'"'
          write(sxml_unit, '(4x, a)') trim(info_str)
          info_str = '        yx="'//str(shieldings(2,1,i_atom))//'" yy="'// &
               & str(shieldings(2,2,i_atom))//'" yz="'// &
               & str(shieldings(2,3,i_atom))//'"'
          write(sxml_unit, '(4x, a)') trim(info_str)
          info_str = '        zx="'//str(shieldings(3,1,i_atom)) //'" zy="'// &
               & str(shieldings(3,2,i_atom))//'" zz="'// &
               & str(shieldings(3,3,i_atom))//'" />'
          write(sxml_unit, '(4x, a)') trim(info_str)
          write(sxml_unit, '(2x, a)') '</interaction>'
          id = id + 1
       end do
    end if
    ! J-couplings
    if (j_is_active) then
       do i_atom = 1, size(active_nuclei)
          do j_atom = i_atom+1, size(active_nuclei)
             info_str = '<interaction kind="jcoupling" id="'//str(id)// &
                  & '" units="Hz" spin_a="'//str(i_atom)//'" spin_b="'// &
                  & str(j_atom)//'" >'
             write(sxml_unit, '(2x, a)') trim(info_str)
             info_str = '<scalar iso="'// &
                  & str(trace(J_couplings(:,:,i_atom,j_atom))/3)//'" />'
             write(sxml_unit, '(4x, a)') trim(info_str)
             write(sxml_unit, '(2x, a)') '</interaction>'
             id = id + 1
          end do
       end do
    end if
    write(sxml_unit, '(a)') '</spin_system>'
    close(sxml_unit)
  end subroutine create_sxml_file

  ! subroutine create_mag_files()
  !   use dimensions,         only: n_spin, n_states_k(1), n_basis
  !   use tools,              only: newunit
  !   use physics,            only: KS_eigenvalue, KS_eigenvector, occ_numbers
  !   !  use scalapack_wrapper,  only: eigenvec

  !   integer :: calph_unit, cbeta_unit, &
  !        & ealph_unit, ebeta_unit, &
  !        & doca_unit, docb_unit
  !   integer :: i_state, i_basis

  !   inquire(file='CALPH.771', exist=exists)
  !   if (exists) then
  !      open(newunit(calph_unit), file='CALPH.771')
  !      close(calph_unit, status='delete')
  !   end if

  !   inquire(file='CBETA.772', exist=exists)
  !   if (exists) then
  !      open(newunit(cbeta_unit), file='CBETA.772')
  !      close(cbeta_unit, status='delete')
  !   end if

  !   inquire(file='EALPH.771', exist=exists)
  !   if (exists) then
  !      open(newunit(ealph_unit), file='EALPH.771')
  !      close(ealph_unit, status='delete')
  !   end if

  !   inquire(file='EBETA.772', exist=exists)
  !   if (exists) then
  !      open(newunit(ebeta_unit), file='EBETA.772')
  !      close(ebeta_unit, status='delete')
  !   end if

  !   inquire(file='DOCA.771', exist=exists)
  !   if (exists) then
  !      open(newunit(doca_unit), file='DOCA.771')
  !      close(doca_unit, status='delete')
  !   end if

  !   inquire(file='DOCB.772', exist=exists)
  !   if (exists) then
  !      open(newunit(docb_unit), file='DOCB.772')
  !      close(docb_unit, status='delete')
  !   end if

  !   if (n_spin == 1) then
  !      open(newunit(calph_unit), file='CALPH.771',action='write',&
  !           & form='formatted')
  !      write(unit=calph_unit,fmt=*) &
  !           & [((KS_eigenvector(i_basis, i_state, 1, 1), &
  !           & i_basis = 1, n_basis), &
  !           & i_state = 1, n_states_k(1))]
  !      close(calph_unit)

  !      open(newunit(cbeta_unit), file='CBETA.772',action='write',&
  !           & form='formatted')
  !      write(unit=cbeta_unit,fmt=*) &
  !           & [((KS_eigenvector(i_basis, i_state, 1, 1), &
  !           & i_basis = 1, n_basis), &
  !           & i_state = 1, n_states_k(1))]
  !      close(cbeta_unit)

  !      open(newunit(ealph_unit), file='EALPH.771',action='write',&
  !           & form='formatted')
  !      write(unit=ealph_unit,fmt=*) &
  !           & [(KS_eigenvalue(i_state, 1, 1), i_state = 1, n_states_k(1))]
  !      close(ealph_unit)

  !      open(newunit(ebeta_unit), file='EBETA.772',action='write',&
  !           & form='formatted')
  !      write(unit=ebeta_unit,fmt=*) &
  !           & [(KS_eigenvalue(i_state, 1, 1), i_state = 1, n_states_k(1))]
  !      close(ebeta_unit)

  !      open(newunit(doca_unit), file='DOCA.771',action='write',&
  !           & form='formatted')
  !      write(unit=doca_unit,fmt=*) &
  !           & [(occ_numbers(i_state, 1, 1)*0.5, i_state = 1, n_states_k(1))]
  !      close(doca_unit)

  !      open(newunit(docb_unit), file='DOCB.772',action='write',&
  !           & form='formatted')
  !      write(unit=docb_unit,fmt=*) &
  !           & [(occ_numbers(i_state, 1, 1)*0.5, i_state = 1, n_states_k(1))]
  !      close(docb_unit)

  !   else if (n_spin == 2) then
  !      ! spin-polarized case

  !      open(newunit(calph_unit), file='CALPH.771',action='write',&
  !           & form='formatted')
  !      write(unit=calph_unit,fmt=*) &
  !           & [((KS_eigenvector(i_basis, i_state, 1, 1), &
  !           & i_basis = 1, n_basis), &
  !           & i_state = 1, n_states_k(1))]
  !      close(calph_unit)

  !      open(newunit(cbeta_unit), file='CBETA.772',action='write',&
  !           & form='formatted')
  !      write(unit=cbeta_unit,fmt=*) &
  !           & [((KS_eigenvector(i_basis, i_state, 2, 1), &
  !           & i_basis = 1, n_basis), &
  !           & i_state = 1, n_states_k(1))]
  !      close(cbeta_unit)

  !      open(newunit(ealph_unit), file='EALPH.771',action='write',&
  !           & form='formatted')
  !      write(unit=ealph_unit,fmt=*) &
  !           & [(KS_eigenvalue(i_state, 1, 1), i_state = 1, n_states_k(1))]
  !      close(ealph_unit)

  !      open(newunit(ebeta_unit), file='EBETA.772',action='write',&
  !           & form='formatted')
  !      write(unit=ebeta_unit,fmt=*) &
  !           & [(KS_eigenvalue(i_state, 2, 1), i_state = 1, n_states_k(1))]
  !      close(ebeta_unit)

  !      open(newunit(doca_unit), file='DOCA.771',action='write',&
  !           & form='formatted')
  !      write(unit=doca_unit,fmt=*) &
  !           & [(occ_numbers(i_state, 1, 1)*0.5, i_state = 1, n_states_k(1))]
  !      close(doca_unit)

  !      open(newunit(docb_unit), file='DOCB.772',action='write',&
  !           & form='formatted')
  !      write(unit=docb_unit,fmt=*) &
  !           & [(occ_numbers(i_state, 2, 1)*0.5, i_state = 1, n_states_k(1))]
  !      close(docb_unit)
  !   end if
  ! end subroutine create_mag_files
end module MR_output
