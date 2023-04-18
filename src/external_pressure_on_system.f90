!****h* FHI-aims/external_pressure_on_system
!  NAME
!   external_pressure 
!  SYNOPSIS
module external_pressure_on_system

!  PURPOSE
!    This module contains routines for applying external pressure to system
!  USES
implicit none

contains

  !! FK: Write summary for external pressure
  subroutine external_pressure_summary(tot_energy, ext_pressure, volume, stress)
    use localorb_io
    use constants
    implicit none

    real*8,                   intent(in)           :: tot_energy   ! has unit Ha
    real*8,                   intent(in)           :: ext_pressure ! has unit Ha/bohr**3
    real*8,                   intent(in)           :: volume       ! has unit bohr**3
    real*8,dimension(1:3,1:3),intent(in), optional :: stress       ! has unit Ha/bohr**3

    ! local
    logical       :: flag_stress = .false.
    real*8        :: int_pressure
    character*120 :: info_str

    if (present(stress)) flag_stress = .true.

    write (info_str,'(2X,A)') &
      "Summary for the variables related to the applied external pressure:"
    call localorb_info(info_str)
    write (info_str,'(2X,A,F15.8,A)') &
      "| Cell volume                   : ", volume * (bohr**3), " A**3"
    call localorb_info(info_str)
    write (info_str,'(2X,A,F15.8,A)') &
      "| External pressure             : ", ext_pressure * hartree/(bohr**3), " eV/A**3"
    call localorb_info(info_str)

    if (flag_stress) then
      int_pressure = (stress(1,1) + stress(2,2) + stress(3,3))/(-3.0d0)
      write (info_str,'(2X,A,F15.8,A)') &
        "| Internal pressure             : ", int_pressure * hartree/(bohr**3), " eV/A**3"
      call localorb_info(info_str)
      write (info_str,'(2X,A,F15.8,A)') &
        "| Total pressure inside system  : ", (int_pressure - ext_pressure) * hartree/(bohr**3), " eV/A**3"
      call localorb_info(info_str)
    end if

    write (info_str,'(2X,A)') &
      "Since the calculation is done at T = 0K, the Gibbs free energy is G = E_tot + p_ext * V:"
    call localorb_info(info_str)
    write (info_str,'(2X,A,E23.15,A)') &
      "| Gibbs free energy             : ", (tot_energy + ext_pressure*volume) * hartree, " eV"
    call localorb_info(info_str)
    call localorb_info('')

  end subroutine external_pressure_summary

  !! FK: Write notice that output geometry has external pressure applied
  subroutine external_pressure_notice(ext_pressure)
    use localorb_io
    use constants
    implicit none

    real*8,intent(in) :: ext_pressure ! has unit Ha/bohr**3

    ! local
    character*120 :: info_str

    write (info_str,'(2X,A,E10.3,A,E10.3,A)') &
       "Reminder: External pressure of ", ext_pressure * hartree/(bohr**3), &
       " eV/A**3 (", ext_pressure * hartree/(bohr**3)*giga_pascal , &
       " GPa) is applied to system."
    call localorb_info(info_str)

  end subroutine external_pressure_notice

  !! FK: Add the external pressure term p*V to energy
  subroutine add_pV_to_energy(energy, ext_pressure, volume)
    implicit none

    real*8,intent(inout) :: energy
    real*8,intent(in)    :: ext_pressure
    real*8,intent(in)    :: volume

    energy = energy + ext_pressure * volume

  end subroutine add_pV_to_energy

  !! FK: Add the external pressure to diagonal of stress tensor
  subroutine add_p_to_stress(stress, ext_pressure)
    implicit none

    real*8,dimension(1:3,1:3),intent(inout) :: stress
    real*8,                   intent(in)    :: ext_pressure

    ! local
    integer :: i_temp

    do i_temp = 1, 3, 1
      stress(i_temp,i_temp) = stress(i_temp,i_temp) + ext_pressure
    end do

  end subroutine add_p_to_stress

end module external_pressure_on_system
