!****s* FHI-aims/renormalize_density
!  NAME
!    renormalize_density
!  SYNOPSIS

subroutine renormalize_density(density, density_gradient, kinetic_density, &
                               density_dimension, partition_tab, n_electrons)

  !  PURPOSE
  !
  !     Takes an input electron density and renormalizes it to exactly n_electrons electrons
  !     as integrated on the normal 3D integration grid.
  !     The partition_tab on input can be anything, including the hartree_partition_tab,
  !     which is the partitioning used for the electrostatic potential. 
  !
  !  USES

  use grids
  use dimensions
  use runtime_choices
  use localorb_io
  use synchronize_mpi_basic
  use geometry, only : total_initial_charge
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: density_dimension
  real*8, intent(INOUT) :: density(density_dimension, n_full_points)
  real*8, intent(INOUT) :: density_gradient(3, density_dimension, n_full_points)
  real*8, intent(INOUT) :: kinetic_density(density_dimension, n_full_points)
  real*8, intent(IN) :: partition_tab(n_full_points)
  real*8, intent(IN) :: n_electrons

  !  INPUTS
  !    o density_dimension -- n_spin for rho, 1 for free_rho_superpos
  !    o density           -- Electron density which will be renormalized
  !    o density_gradient  -- Corresponding gradient to be renormalized
  !    o partition_tab -- Integration weights
  !    o n_electrons   -- REAL*8 (!) number of electrons that should be in the system
  !                       according to input file control.in
  !  OUTPUTS
  !    o density and density_gradient will be modified
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2012).
  !  SOURCE

  integer :: info
  character*150 :: info_str
  character(*), parameter :: func = 'renormalize_density'

  integer :: i_my_batch, i_index, i_spin
  integer :: i_full_point, i_atom
  real*8 :: weight, this_density
  real*8 :: electron_count
  real*8 :: norm

  write(info_str, '(2X,A)') &
    'Renormalizing the density to the exact electron count on the 3D integration grid.' 
  call localorb_info(info_str, use_unit, '(A)', OL_norm)

  if (n_electrons.gt.0.d0) then

    electron_count = 0.d0

    ! simplest possible loop over the integration grid

    i_full_point = 0
    do i_my_batch = 1, n_my_batches
      do i_index = 1, batches(i_my_batch)%size
        i_full_point = i_full_point + 1
        weight = partition_tab(i_full_point)
        this_density = 0.d0
        do i_spin = 1, density_dimension, 1
           this_density = this_density+density(i_spin, i_full_point)
        enddo
        
        electron_count &
        & = electron_count + weight * this_density

      end do
    end do

    call sync_real_number(electron_count)

    if (electron_count.gt.0.d0) then
      norm = n_electrons / electron_count
    else
      write(info_str, '(2X,A,F18.10)') &
        '*** Error : The actual electron count is', electron_count
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      write(info_str, '(2X,A)') &
        '*** We will continue without renormalizing anything.'
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      norm = 1.d0
    end if
 
    write(info_str, '(2X,A,F18.10)') &
      '| Formal number of electrons (from input files) : ', n_electrons   
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
    write(info_str, '(2X,A,F18.10)') &
      '| Integrated number of electrons on 3D grid     : ', electron_count
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
    write(info_str, '(2X,A,F18.10)') &
      '| Charge integration error                      : ', electron_count-n_electrons
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
    write(info_str, '(2X,A,F18.10)') &
      '| Normalization factor for density and gradient : ', norm
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
  
    ! Corner case. In the case of an initial density that does not have the
    ! same charge of the system, we can not safely renormalize.
    ! As s.c.f. cycles pass, the density will gradually approach the correct
    ! norm, and then we can start to renormalize.
    !
    ! We hardwire that the initial and formal charge are equal if they differ
    ! by less that 10^-5 electrons (also in read_geo,f90!), and that we will
    ! only start to renormalize once the renormalization affects fewer than
    ! one in a thousand electrons.
    !
    if ( ( dabs(total_initial_charge - charge) .gt. 1.d-5 ) & 
         .and. ( dabs(norm - 1.d0) .gt. 1.d-3 ) ) then

      write(info_str, '(2X,A)') &
        '* The proposed renormalization factor deviates too far from 1.'
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      write(info_str, '(2X,A)') &
        '* Since the initial charge in geometry.in was not equal to the formal charge,'
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      write(info_str, '(2X,A)') &
        '* it is not yet safe to renormalize the charge. If this message persists'
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      write(info_str, '(2X,A)') &
        '* even in the final s.c.f. cycle, your integration accuracy may not be good enough.'
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      
    else
      ! hopefully it is safe to renormalize

      ! Now renormalize the actual density and gradient
      density = density*norm
      if (use_density_gradient) then
        density_gradient = density_gradient*norm
        if (use_meta_gga) then
          kinetic_density = kinetic_density*norm
        end if
      endif

    end if

  else ! zero or fewer electrons?

      write(info_str, '(2X,A,F18.10)') &
        '*** Error : The formal number of electrons is n_electrons = ', n_electrons
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      write(info_str, '(2X,A)') &
        '*** We will continue without renormalizing anything.'
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      norm = 1.d0

  end if

end subroutine renormalize_density
!******
!****s* FHI-aims/renormalize_free_density
!  NAME
!    renormalize_density
!  SYNOPSIS

subroutine renormalize_free_density(density, density_gradient, density_dimension, partition_tab, n_electrons)

  !  PURPOSE
  !
  !     WARNING! THIS ROUTINE SHOULD ONLY BE USED TO RENORMALIZE free_rho_superpos,
  !     AS IT ALSO RENORMALIZES THE FREE ATOM SPLINES.
  !
  !     Takes an input electron density and renormalizes it to exactly n_electrons electrons
  !     as integrated on the normal 3D integration grid.
  !     The partition_tab on input can be anything, including the hartree_partition_tab,
  !     which is the partitioning used for the electrostatic potential. 
  !
  !     This subroutine must be added specifically to handle free_rho_superpos,
  !     which is multiplied by a factor 4*pi compared to the normal rho.
  !
  !  USES

  use constants
  use grids
  use dimensions
  use runtime_choices
  use localorb_io
  use synchronize_mpi_basic
  use free_atoms
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: density_dimension
  real*8, intent(INOUT) :: density(density_dimension, n_full_points)
  real*8, intent(INOUT) :: density_gradient(3, density_dimension, n_full_points)
  real*8, intent(IN) :: partition_tab(n_full_points)
  real*8, intent(IN) :: n_electrons

  !  INPUTS
  !    o density_dimension -- n_spin for rho, 1 for free_rho_superpos
  !    o density           -- Electron density which will be renormalized
  !    o density_gradient  -- Corresponding gradient to be renormalized
  !    o partition_tab -- Integration weights
  !    o n_electrons   -- REAL*8 (!) number of electrons that should be in the system
  !                       according to input file control.in
  !  OUTPUTS
  !    o density and density_gradient will be modified
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2012).
  !  SOURCE

  integer :: info
  character*150 :: info_str
  character(*), parameter :: func = 'renormalize_free_density'

  integer :: i_my_batch, i_index, i_spin
  integer :: i_full_point, i_atom
  real*8 :: weight, this_density
  real*8 :: electron_count
  real*8 :: norm

  write(info_str, '(2X,A)') &
    'Renormalizing the free-atom superposition density to the exact electron count on the 3D integration grid.' 
  call localorb_info(info_str, use_unit, '(A)', OL_norm)

  ! only normalize if there are more than zero electrons in the neutral density!
  if ( (n_electrons+charge).gt.0.d0 ) then

    electron_count = 0.d0

    ! simplest possible loop over the integration grid

    i_full_point = 0
    do i_my_batch = 1, n_my_batches
      do i_index = 1, batches(i_my_batch)%size
        i_full_point = i_full_point + 1
        weight = partition_tab(i_full_point)
        this_density = 0.d0
        do i_spin = 1, density_dimension, 1
           this_density = this_density+density(i_spin, i_full_point)
        enddo
        
        electron_count &
        & = electron_count + weight * pi4_inv * this_density

      end do
    end do

    call sync_real_number(electron_count)

    ! in this case, the free atom density must be renormalized to the formally neutral system,
    ! i.e., n_electrons - (-charge)

    if (electron_count.gt.0.d0) then
      norm = (n_electrons + charge) / electron_count
    else
      write(info_str, '(2X,A,F18.10)') &
        '*** Error : The actual electron count is', electron_count
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      write(info_str, '(2X,A)') &
        '*** We will continue without renormalizing anything.'
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      norm = 1.d0
    end if

    write(info_str, '(2X,A,F18.10)') &
      '| Formal number of electrons (from input files) : ', n_electrons  
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
    write(info_str, '(2X,A,F18.10)') &
      '| Integrated number of electrons on 3D grid     : ', electron_count
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
    write(info_str, '(2X,A,F18.10)') &
      '| Charge integration error                      : ', & 
      electron_count- ( n_electrons + charge ) 
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
    write(info_str, '(2X,A,F18.10)') &
      '| Normalization factor for density and gradient : ', norm
    call localorb_info(info_str, use_unit, '(A)', OL_norm)

    ! Now renormalize the actual density and gradient
    density = density*norm
    if (use_density_gradient) then
      density_gradient = density_gradient*norm
    end if

    ! Since this here is the free-atom version ONLY, we also renormalize the free-atom density splines
    renormalized_free_rho_spl = free_rho_spl * norm
    if (use_density_gradient) then
      renormalized_free_drho_dr_spl = free_drho_dr_spl * norm
    end if

  else ! zero or fewer electrons in neutral free-atom superposition ?

      write(info_str, '(2X,A)') &
        '*** Error : The formal number of electrons in the free-atom superposition density '
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      write(info_str, '(2X,A,F18.10)') &
        '*** is n_electrons+charge = ', & 
       n_electrons+charge
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      write(info_str, '(2X,A)') &
        '*** We will continue without renormalizing anything.'
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      norm = 1.d0

  end if

end subroutine renormalize_free_density
!******
!****s* FHI-aims/renormalize_initial_density
!  NAME
!    renormalize_initial_density
!  SYNOPSIS

subroutine renormalize_initial_density& 
   (density, density_gradient, density_dimension, partition_tab, n_electrons)

  !  PURPOSE
  !
  !  Only for the initial density. In a charged system, the formal charge of the
  !  initial density is known exactly.
  !
  !  USES

  use grids
  use dimensions
  use runtime_choices
  use localorb_io
  use synchronize_mpi_basic
  use geometry, only : total_initial_charge
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: density_dimension
  real*8, intent(INOUT) :: density(density_dimension, n_full_points)
  real*8, intent(INOUT) :: density_gradient(3, density_dimension, n_full_points)
  real*8, intent(IN) :: partition_tab(n_full_points)
  real*8, intent(IN) :: n_electrons

  !  INPUTS
  !    o density_dimension -- n_spin for rho, 1 for free_rho_superpos
  !    o density           -- Electron density which will be renormalized
  !    o density_gradient  -- Corresponding gradient to be renormalized
  !    o partition_tab -- Integration weights
  !    o n_electrons   -- REAL*8 (!) number of electrons that should be in the system
  !                       according to input file control.in
  !  OUTPUTS
  !    o density and density_gradient will be modified
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2012).
  !  SOURCE

  integer :: info
  character*150 :: info_str
  character(*), parameter :: func = 'renormalize_initial_density'

  integer :: i_my_batch, i_index, i_spin
  integer :: i_full_point, i_atom
  real*8 :: weight, this_density
  real*8 :: electron_count
  real*8 :: norm

  write(info_str, '(2X,A)') &
    'Renormalizing the initial density to the exact electron count on the 3D integration grid.' 
  call localorb_info(info_str, use_unit, '(A)', OL_norm)

  if (n_electrons.gt.0.d0) then

    electron_count = 0.d0

    ! simplest possible loop over the integration grid

    i_full_point = 0
    do i_my_batch = 1, n_my_batches
      do i_index = 1, batches(i_my_batch)%size
        i_full_point = i_full_point + 1
        weight = partition_tab(i_full_point)
        this_density = 0.d0
        do i_spin = 1, density_dimension, 1
           this_density = this_density+density(i_spin, i_full_point)
        enddo
        
        electron_count &
        & = electron_count + weight * this_density

      end do
    end do

    call sync_real_number(electron_count)

    if (electron_count.gt.0.d0) then
      norm = ( n_electrons + charge - total_initial_charge ) / electron_count
    else
      write(info_str, '(2X,A,F18.10)') &
        '*** Error : The actual electron count is', electron_count
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      write(info_str, '(2X,A)') &
        '*** We will continue without renormalizing anything.'
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      norm = 1.d0
    end if
 
    write(info_str, '(2X,A,F18.10)') &
      '| Initial density: Formal number of electrons (from input files) : ', & 
      n_electrons + charge - total_initial_charge
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
    write(info_str, '(2X,A,F18.10)') &
      '| Integrated number of electrons on 3D grid     : ', electron_count
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
    write(info_str, '(2X,A,F18.10)') &
      '| Charge integration error                      : ', & 
      electron_count-( n_electrons + charge - total_initial_charge ) 
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
    write(info_str, '(2X,A,F18.10)') &
      '| Normalization factor for density and gradient : ', norm
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
  
    ! Now renormalize the actual density and gradient
    density = density*norm
    if (use_density_gradient) then
      density_gradient = density_gradient*norm
    end if

  else ! zero or fewer electrons?

      write(info_str, '(2X,A,F18.10)') &
        '*** Error : The formal number of electrons is n_electrons = ', n_electrons
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      write(info_str, '(2X,A)') &
        '*** We will continue without renormalizing anything.'
      call localorb_info(info_str, use_unit, '(A)', OL_norm)
      norm = 1.d0

  end if

end subroutine renormalize_initial_density
!******
