!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Provides the nuclear magnetic properties of the most abundant
!!  isotopes used in magnetic response experiments.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
module nuclear_data

  use tools, only: str
  use types, only: dp

  implicit none

  private
  public :: nucleus_t
  public :: default_isotope_number, get_nuclear_data

  type :: nucleus_t
     character(6) :: name
     real(dp) :: spin, mu, q_in_barn
  end type nucleus_t

  ! Magnetic data taken from
  !   N.J. Stone, Atomic Data and Nuclear Data Tables 90, 75 (2005)
  ! Magnetic moments are given in bohr magnetons in that paper. Here
  ! they have been converted to a.u. by multiplying with mu_N in
  ! a.u. = 5.44617021352d-4/2 (physics.nist.gov/cuu/Constants, 2016).
  ! Nuclear quadrupole moments in Barn (=10^-28 m^2) taken from
  !   Pekka Pyykkoe, Molecular Physics 106, 1965 (2008)
  integer, parameter :: N_NUCLEI = 21
  type(nucleus_t), parameter :: nuclei(N_NUCLEI) = [ &
       & nucleus_t(' H-1  ', 1/2d0,  7.605160997d-4, 0d0), &
       & nucleus_t(' H-2  ', 1/1d0,  2.334877269d-4, 2.860d-3), &
       & nucleus_t('He-3  ', 1/2d0, -5.793357356d-4, 0d0), &
       & nucleus_t('Li-7  ', 3/2d0,  8.867527865d-4, -40.1d-3), &
       & nucleus_t('Be-9  ', 3/2d0, -3.206247543d-4, 52.88d-3), &
       & nucleus_t(' B-11 ', 3/2d0,  7.321419777d-4, 40.59d-3), &
       & nucleus_t(' C-11 ', 3/2d0, -2.625054043d-4, 33.27d-3), &
       & nucleus_t(' C-13 ', 1/2d0,  1.912727111d-4, 0d0), &
       & nucleus_t(' N-14 ', 1/1d0, -1.099475566d-4, 20.44d-3), &
       & nucleus_t(' N-15 ', 1/2d0, -7.711473126d-5, 0d0), &
       & nucleus_t(' O-17 ', 5/2d0, -5.156951344d-4, -25.58d-3), &
       & nucleus_t(' F-19 ', 1/2d0,  7.158631298d-4, -94.2d-3), &
       & nucleus_t('Ne-21 ', 3/2d0, -1.802129554d-4, 101.55d-3), &
       & nucleus_t('Na-23 ', 3/2d0,  6.038501132d-4, 104d-3), &
       & nucleus_t('Mg-25 ', 5/2d0, -2.329463155d-4, 199.420d-3), &
       & nucleus_t('Al-27 ', 5/2d0,  9.916133206d-4, 146.610d-3), &
       & nucleus_t('Si-29 ', 1/2d0, -1.512101929d-4, 0d0), &
       & nucleus_t(' P-31 ', 1/2d0,  3.081443107d-4, 0d0), &
       & nucleus_t(' S-33 ', 3/2d0,  1.753179921d-4, -67.813d-3), &
       & nucleus_t('Cl-35 ', 3/2d0,  2.238033666d-4, -81.6580d-3), &
       & nucleus_t('Ar-39 ', 7/2d0, -3.155511022d-4, 0d0) &
       & ]

  ! Default isotope numbers depending on the type of calculation.
  ! These are in general different for NMR and EFG.
  type :: isotope_t
     character(2) :: name
     integer :: isotope_NMR, isotope_EFG
  end type isotope_t

  integer, parameter :: N_DEFAULT_NUCLEI = 18
  type(isotope_t), parameter :: default_isotopes(N_DEFAULT_NUCLEI) = [ &
       & isotope_t(' H',  1,  2), & ! First column NMR, second column EFG
       & isotope_t('He',  3,  3), &
       & isotope_t('Li',  7,  7), &
       & isotope_t('Be',  9,  9), &
       & isotope_t(' B', 11, 11), &
       & isotope_t(' C', 13, 11), &
       & isotope_t(' N', 15, 14), &
       & isotope_t(' O', 17, 17), &
       & isotope_t(' F', 19, 19), &
       & isotope_t('Ne', 21, 21), &
       & isotope_t('Na', 23, 23), &
       & isotope_t('Mg', 25, 25), &
       & isotope_t('Al', 27, 27), &
       & isotope_t('Si', 29, 29), &
       & isotope_t(' P', 31, 31), &
       & isotope_t(' S', 33, 33), &
       & isotope_t('Cl', 35, 35), &
       & isotope_t('Ar', 39, 39) &
       & ]

contains

  !!  FUNCTION
  !!
  !!  Returns nuclear data corresponding to isotope 'isotope' of
  !!  element 'name'. See MR_main() for usage.
  !!
  !!  INPUTS
  !!
  !!  name - name of the element, e.g., 'He'.
  !!  isotope - isotope number, e.g., '13' for C-13.
  !!
  elemental type(nucleus_t) function get_nuclear_data(name, isotope) result(y)
    character(*), intent(in) :: name
    integer, intent(in) :: isotope
    integer :: i_atom
    do i_atom = 1, size(nuclei)
       if (trim(name)//'-'//str(isotope) == &
            & trim(adjustl(nuclei(i_atom)%name))) then
          y = nuclei(i_atom)
          return
       end if
    end do
    ! Default if not found
    y = nucleus_t(' ---', 1d0, 0d0, 0d0)
  end function get_nuclear_data

  !!  FUNCTION
  !!
  !!  Unless specified in geometry.in with the 'isotope' keyword,
  !!  returns the default isotope number depending on the element. See
  !!  MR_main() for usage.
  !!
  !!  INPUTS
  !!
  !!  name - name of the element, e.g., 'He'.
  !!  calc_type - For now, 'NMR' or 'EFG'.
  !!
  elemental integer function default_isotope_number(name, calc_type) result(y)
    character(*), intent(in) :: name, calc_type
    integer :: i_atom
    do i_atom = 1, size(default_isotopes)
       if (trim(name) == trim(adjustl(default_isotopes(i_atom)%name))) then
          if (calc_type == 'NMR' .or. calc_type == 'nmr') then
             y = default_isotopes(i_atom)%isotope_NMR
          else if (calc_type == 'EFG' .or. calc_type == 'efg') then
             y = default_isotopes(i_atom)%isotope_EFG
          end if
          return
       end if
    end do
    y = -1
  end function default_isotope_number
end module nuclear_data
