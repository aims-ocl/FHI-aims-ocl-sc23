!****s* FHI-aims/output_KS_dos
!  NAME
!   output_KS_dos
!  SYNOPSIS

subroutine output_KS_dos &
     ( KS_eigenvalue, occ_numbers, chemical_potential, n_electrons, &
       filename_shifted, filename_raw )

! PURPOSE
! The subroutine prints to file the Kohn-Sham density of states, derived from an advanced version of 
! output_mulliken.f90.
!
! Note that the present version, same as the partial density of states, does NOT include
! support scalapack. This means scalapack type of eigenvectors. If normal form of eigenvectors
! are created from scalapack routine then this this can be called.
!
! Note also that the output of the SOC perturbed DOS is not found here,
! because this routine is called at the end of scf_solver(), before the SOC
! perturbed eigenvalues have been calculated.  Much as I (WPH) hate doing so,
! this forces me to fork the DOS outputting for this case, and it can be found in
! output_KS_dos_soc()
!
!  USES

  use runtime_choices,    only : dos_alpha, dos_n_en_points, &
                                 dos_low_energy, dos_high_energy
  use constants,          only : hartree, one_over_sqrt2
  use dimensions,         only : n_states, n_spin, n_k_points, spin_degeneracy
  use mpi_tasks,          only : myid, n_tasks
  use localorb_io,        only : use_unit, OL_norm, localorb_info
  use arch_specific,      only : arch_erf
  use pbc_lists,          only : k_weights
  use synchronize_mpi,    only : sync_vector 
  use timing,             only : get_times, get_timestamps, &
                                 tot_time_band_dos, tot_clock_time_band_dos
  use generate_aims_uuid, only : write_aims_uuid

  implicit none

!  ARGUMENTS

  real*8, dimension(n_states, n_spin, n_k_points)               :: KS_eigenvalue 
  real*8, dimension(n_states, n_spin,n_k_points)                :: occ_numbers
  real*8:: chemical_potential
  real*8  :: n_electrons
  character(len=*) :: filename_shifted
  character(len=*) :: filename_raw

!  INPUTS
!   o KS_eigenvalue -- Kohn-Sham eigenvalues
!   o occ_numbers -- occupations weights of the Kohn-Sham eigenstates
!   o n_electrons -- total number of electrons in the system
!  OUTPUT
!   o chemical_potential -- chemical potential
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

  ! local variables
  real*8  :: diff_electrons
  real*8  :: de,  E1, E2
  real*8  :: en
  real*8, dimension(:,:), allocatable :: KS_dos

  ! counters
  integer :: i_spin, i_state, i_counter
  integer :: i_k_point, i_e
  character*200 :: info_str

  real*8 :: dos_time = 0.d0
  real*8 :: clock_dos_time = 0.d0

  call get_timestamps(dos_time, clock_dos_time)

  !!!VB: I removed the call to check_norm , for now.
  !      * For one thing, I do not see how occ_numbers are used in this routine.
  !      * Then, the Fermi level may be spin dependent if a spin constraint is used - the 
  !      results below are then wrong.
  !      * The occupation numbers must be updated whenever the KS eigenvalues change,
  !        right away, not here
  !      * Finally, if the KS eigenvalues changed, the Fermi level should itself be updated. The
  !        routine below only recomputes occupation numbers for a fixed Fermi level, the Fermi level
  !        does not change.
  !
  ! calculate position of Fermi level, among other things ... 
  ! i_counter = 0
  ! call check_norm_p0( chemical_potential, KS_eigenvalue, n_electrons, occ_numbers, diff_electrons, i_counter)

  if(.not. allocated (KS_dos)) allocate(KS_dos(n_spin,dos_n_en_points)) 
  
  write(info_str,'(2X,A)') 'Calculating total density of states ...'
  call localorb_info(info_str, use_unit)
  write(info_str,'(2X,A,F10.6,A)') '| Chemical potential is ', chemical_potential* hartree,' eV.'
  call localorb_info(info_str, use_unit)

  de= (dos_high_energy - dos_low_energy)/dble(dos_n_en_points-1)
  KS_dos =0d0 
  do i_k_point = 1, n_k_points
     if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
        do i_spin = 1, n_spin
           do i_e = 1, dos_n_en_points
              en= dos_low_energy + dble(i_e-1)*de   
              E1 = en - de/2d0
              E2 = en + de/2d0
              do i_state = 1, n_states, 1
                 KS_dos(i_spin,i_e) =  KS_dos (i_spin,i_e) + spin_degeneracy * k_weights(i_k_point)&
                      *(arch_erf((E2-(KS_eigenvalue(i_state,i_spin,i_k_point))*hartree)*one_over_sqrt2/dos_alpha)&
                      -arch_erf((E1-(KS_eigenvalue(i_state,i_spin,i_k_point))*hartree)*one_over_sqrt2/dos_alpha)&
                      )/(2d0*dE)
              enddo
           enddo
        end do
     end if
  enddo
  
  call sync_vector(KS_dos(1,1),n_spin*dos_n_en_points)

  if(myid.eq.0)then
     ! arrange output for DOS: different for spin polarized and normal calculations
     write(info_str,'(2X,A,A)') '| writing DOS (shifted by electron chemical potential) to file ', &
          trim(filename_shifted)
     call localorb_info(info_str, use_unit)
     write(info_str,'(2X,A,A)') '| writing DOS (raw data) to file ', &
          trim(filename_raw)
     call localorb_info(info_str, use_unit)

     open(88, file=trim(filename_shifted))
     open(89, file=trim(filename_raw))
     if (n_spin.eq.2) then
        write(88,'(2A)') '# total density of states output by FHI-aims, ',&
             'for a spin polarized system' 
        write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', &
             chemical_potential*hartree, ' eV'
        call write_aims_uuid(info_str)
        write(88,'(A,2X,A)') '#', info_str

        write(88,'(A)') '#    Energy (eV)      DOS(spin up)   DOS(spin down)'
        write(89,'(2A)') '# total density of states output by FHI-aims, ',&
             'for a spin polarized system' 
        write(89,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
        call write_aims_uuid(info_str)
        write(89,'(A,2X,A)') '#', info_str
        write(89,'(A)') '#    Energy (eV)      DOS(spin up)   DOS(spin down)'
        do i_e = 1, dos_n_en_points ,1
           en = dos_low_energy + dble(i_e-1)*de 
           write(88,'(3F20.8)') en-chemical_potential*hartree, KS_dos(1,i_e), KS_dos(2,i_e)
           write(89,'(3F20.8)') en, KS_dos(1,i_e), KS_dos(2,i_e)
        enddo
     else
        write(88,'(2A)') '# total density of states output by FHI-aims '
        write(88,'(A,F15.6,A)') '# The energy reference for this output is the chemical potential, mu = ', &
             chemical_potential*hartree, ' eV'
        write(88,'(A)') '#    Energy (eV)      DOS '
        write(89,'(2A)') '# total density of states output by FHI-aims '
        write(89,'(A,F15.6,A)') '# The energy reference for this output is the vacuum level'
        write(89,'(A)') '#    Energy (eV)      DOS '
        do i_e = 1, dos_n_en_points ,1
           en = dos_low_energy + dble(i_e-1)*de 
           write(88,'(3F20.8)') en-chemical_potential*hartree, KS_dos(1,i_e)
           write(89,'(3F20.8)') en, KS_dos(1,i_e)
        enddo
        
     end if
     close(88)
     close(89)

  end if

  if (allocated(KS_dos)) deallocate(KS_dos)

  ! I'm (WPH) choosing here not to output the time, as this should be a fast
  ! calculation, but for the sake of completeness I'm keeping track of the
  ! timing anyway
  call get_times(dos_time, clock_dos_time, tot_time_band_dos, tot_clock_time_band_dos, .true.)

end subroutine output_KS_dos
!******	
