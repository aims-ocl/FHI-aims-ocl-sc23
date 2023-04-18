!****s* FHI-aims/output_KS_dos_tetrahedron
!  NAME
!   output_KS_dos_tetrahedron
!  SYNOPSIS

subroutine output_KS_dos_tetrahedron &
     ( KS_eigenvalue, occ_numbers, chemical_potential, n_electrons, &
       filename_shifted, filename_raw )

! PURPOSE
! The subroutine prints to file the Kohn-Sham density of states with tetrahedron method.
! derived from output_KS_dos.f90
!
!  USES

  use runtime_choices,    only : dos_alpha, dos_n_en_points, &
                                 dos_low_energy, dos_high_energy, &
                                 n_k_points_xyz_nosym, &
                                 k_points_offset
  use constants,          only : hartree, one_over_sqrt2
  use dimensions,         only : n_states, n_spin, n_k_points, &
                                 spin_degeneracy, ik2irred_map
  use mpi_tasks,          only : myid, n_tasks, mpi_comm_global
  use localorb_io,        only : use_unit, OL_norm, localorb_info
  use arch_specific,      only : arch_erf
  use pbc_lists,          only : k_weights
  use synchronize_mpi,    only : sync_vector 
  use timing,             only : get_times, get_timestamps, &
                                 tot_time_band_dos, tot_clock_time_band_dos
  use geometry,           only : lattice_vector
  use generate_aims_uuid, only : write_aims_uuid
  use tetrahedron_integration, only : ltidos

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
!   Yi Yao
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
!    Yi Yao added April. 2019
!  SOURCE

  ! local variables
  real*8  :: diff_electrons
  real*8  :: de,  E1, E2
  real*8  :: en
  real*8, dimension(:,:), allocatable :: KS_dos
  real*8  :: energies(dos_n_en_points)
  real*8  :: eigs(n_states,n_k_points_xyz_nosym(1),n_k_points_xyz_nosym(2),n_k_points_xyz_nosym(3),n_spin)
  real*8  :: r_x, r_y, r_z

  ! counters
  integer :: i_spin, i_state, i_counter
  integer :: i_k_point, i_e
  character*200 :: info_str
  integer :: i_x, i_y, i_z

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
  
  ! assert n_periodic == 3

  if(.not. allocated (KS_dos)) allocate(KS_dos(n_spin,dos_n_en_points)) 
  
  write(info_str,'(2X,A)') 'Calculating total density of states ...'
  call localorb_info(info_str, use_unit)
  write(info_str,'(2X,A,F10.6,A)') '| Chemical potential is ', chemical_potential* hartree,' eV.'
  call localorb_info(info_str, use_unit)

  !i_k_point = 0
  do i_x = 1, n_k_points_xyz_nosym(1)
    do i_y = 1, n_k_points_xyz_nosym(2)
        do i_z = 1, n_k_points_xyz_nosym(3)
          !i_k_point = i_k_point + 1
          !r_x = dble(i_x-1) / dble(n_k_points_xyz_nosym(1)) + k_points_offset(1)
          !r_y = dble(i_y-1) / dble(n_k_points_xyz_nosym(2)) + k_points_offset(2)
          !r_z = dble(i_z-1) / dble(n_k_points_xyz_nosym(3)) + k_points_offset(3)
          i_k_point = ik2irred_map(i_x,i_y,i_z)
          !if(myid.eq.0) then
          !  write(info_str,'(2X,I)') i_k_point
          !  call localorb_info(info_str, use_unit)
          !end if
          do i_spin = 1, n_spin
            eigs(:,i_x,i_y,i_z,i_spin) = KS_eigenvalue(:,i_spin,i_k_point) * hartree
          end do
        enddo
    enddo
  enddo
  de= (dos_high_energy - dos_low_energy)/dble(dos_n_en_points-1)
  do i_e = 1, dos_n_en_points
    energies(i_e) = dos_low_energy + de * (i_e - 1)
  enddo
  do i_spin = 1, n_spin
    call ltidos(lattice_vector, n_k_points_xyz_nosym(1), n_k_points_xyz_nosym(2), &
                n_k_points_xyz_nosym(3), n_states, &
                eigs(:,:,:,:,i_spin), dos_n_en_points, energies, KS_dos(i_spin,:), mpi_comm_global)
  enddo
  KS_dos = KS_dos * spin_degeneracy

  !KS_dos =0d0 
  !do i_k_point = 1, n_k_points
  !   if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
  !      do i_spin = 1, n_spin
  !         do i_e = 1, dos_n_en_points
  !            en= dos_low_energy + dble(i_e-1)*de   
  !            E1 = en - de/2d0
  !            E2 = en + de/2d0
  !            do i_state = 1, n_states, 1
  !               KS_dos(i_spin,i_e) =  KS_dos (i_spin,i_e) + spin_degeneracy * k_weights(i_k_point)&
  !                    *(arch_erf((E2-(KS_eigenvalue(i_state,i_spin,i_k_point))*hartree)*one_over_sqrt2/dos_alpha)&
  !                    -arch_erf((E1-(KS_eigenvalue(i_state,i_spin,i_k_point))*hartree)*one_over_sqrt2/dos_alpha)&
  !                    )/(2d0*dE)
  !            enddo
  !         enddo
  !      end do
  !   end if
  !enddo
  
  !call sync_vector(KS_dos(1,1),n_spin*dos_n_en_points)

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

end subroutine output_KS_dos_tetrahedron
!******	
