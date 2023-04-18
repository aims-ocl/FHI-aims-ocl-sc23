!****f* FHI-aims/atom_sphere_wrapper
!*  NAME
!*    atom_sphere_wrapper
subroutine atom_sphere_wrapper(n_grid, r_grid, species_name, flag_xc, hybrid_coeff, hse_omega, &
     nuclear_charge, n_core, n_valence, n_max_ind_fns, n_shell, l_shell, s_shell, occup_shell, n_max_grid, &
     n_l_channels, n_spin, &
     cut_atom, cutoff_type, r_cut, scale_cutoff, w_cutoff, &
     eigenval, wave, wave_deriv, density, density_deriv, density_2nd_deriv, atom_pot, atom_pot_es, &
     kinetic, psi_at_0)
!*  SYNOPSIS
!*    call atom_sphere_wrapper(...)
!*  FUNCTION
!*    Provides a wrapper around atom_sphere's atomsolver subroutine to abstract out its idiosyncrasies 
!*    and prevent further cluttering of get_free_atom and its derivatives
!*  USES
  use libxcModule, only : libxc_functionals_init, libxc_functionals_end 
  use spline,      only : val_spline, val_spline_deriv, val_spline_2nd_deriv, cubic_spline
  use runtime_choices, only: flag_rel, REL_atomic_zora
  use mpi_tasks,   only : aims_stop
  use localorb_io, only : localorb_info
  implicit none
!*  INPUTS
!*    o n_grid 
!*    o r_grid 
!*    o species_name 
!*    o flag_xc 
!*    o hybrid_coeff 
!*    o hse_omega
!*    o nuclear_charge
!*    o n_core
!*    o n_valence
!*    o n_shell
!*    o l_shell
!*    o occup_shell
!*    o n_max_ind_fns
!*    o n_max_grid
!*    o n_l_channels
!*    o n_spin
!*    o cut_atom
!*    o cutoff_type
!*    o r_cut
!*    o scale_cutoff
!*    o w_cutoff
!*  OUTPUTS
!*    o eigenval
!*    o wave
!*    o wave_deriv
!*    o density
!*    o density_deriv
!*    o density_2nd_deriv
!*    o atom_pot
!*    o atom_pot_es
!*    o kinetic
!*    o psi_at_0
!*  AUTHORS
!*    William Huhn
!*  NOTES
!*    o This code is not thread-safe, as a random guess is needed to initialize
!*      unoccupied states.  This will lead to errors if not accounted for, as
!*      different ranks will have different basis functions.  To get around
!*      this, synchronization of outputs should occur after this routine is called.
!*    o This subroutine is a heavily modified version of the main subroutine for
!*      atom_sphere.  The parts that look Fortran 95-y are mine, the parts that 
!*      look FORTRAN 77-y are theirs.
!*    o The suffix "as" on variable names refers to quantities evaluated on 
!*      atom_sphere's rr grid, but using aim's definitions (for example, wave = rr * psi)
!*    o The suffix "test" on variable names refers to quantities splined onto
!*      aim's r_grid
!*      
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is subject
!*    the terms and conditions of the respective license agreement."
!* SOURCE

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ARGUMENTS TO/FROM AIMS !
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These are the variables that will be passed to/from aims
  integer, intent(in) :: n_grid
  real*8, intent(in), dimension(n_grid) :: r_grid
  character*2, intent(in) :: species_name
  integer, intent(in) :: flag_xc
  real*8, intent(in) :: hybrid_coeff
  real*8, intent(in) :: hse_omega
  real*8, intent(in) :: nuclear_charge 
  integer, intent(in) :: n_core
  integer, intent(in) :: n_valence
  integer, intent(in) :: n_max_ind_fns
  integer, intent(in), dimension(n_max_ind_fns) :: n_shell
  integer, intent(in), dimension(n_max_ind_fns) :: l_shell
  integer, intent(in), dimension(n_max_ind_fns) :: s_shell
  real*8,  intent(in), dimension(n_max_ind_fns) :: occup_shell
 integer, intent(in) :: n_max_grid
 ! Note: aims combines the l-channels and spin channels into one
 ! variable n_channels.  Well, it claims to.  It actually always
 ! sets n_channels = 1, and then passes n_spin in its place when
 ! it needs spin-polarized values.
 ! While I understand why these two quantities would be combined
 ! into one index for a potential, they are quite different things.
 ! I am pretty sure that, when l-channel-dependent potentials are
 ! implemented in aims, this scheme will be quickly abandonded.
 ! Hell, for all intents and purposes, it already has.
 ! In this interface, I divide the two concepts.
 integer, intent(in) :: n_l_channels
 integer, intent(in) :: n_spin
 logical, intent(in) :: cut_atom
 integer, intent(in) :: cutoff_type
 real*8, intent(in) :: r_cut
 real*8, intent(in) :: scale_cutoff
 real*8, intent(in) :: w_cutoff
 ! Output variables, needed by aims to construct the basis set
 real*8, intent(out), dimension(n_max_ind_fns) :: eigenval
 real*8, intent(out), dimension(n_max_grid, n_max_ind_fns) :: wave
 real*8, intent(out), dimension(n_max_grid, n_max_ind_fns) :: wave_deriv
 real*8, intent(out), dimension(n_max_grid, n_max_ind_fns) :: kinetic
 ! Wavefunction value at the nucleus
 real*8, intent(out) :: psi_at_0(n_max_ind_fns)
 real*8, intent(out), dimension(n_max_grid, n_spin) :: density
 real*8, intent(out), dimension(n_max_grid, n_spin) :: density_deriv
 ! Note:  I (WPH) don't trust the values for density_2nd_deriv coming out of
 ! either sratom or the splined interpolation of atom_sphere near the nucleus.
 ! That being said, atom_sphere looks closer to being correct (sratom's appears
 ! to be garbage) and the variable in question doesn't look to be used anywhere
 real*8, intent(out), dimension(n_max_grid, n_spin) :: density_2nd_deriv
 real*8, intent(out), dimension(n_max_grid, n_l_channels*n_spin) :: atom_pot
 ! I refuse to give the electrostatic and Hartree potentials an l-channel/spin
 ! dependence, even though the main code does.
 real*8, intent(out), dimension(n_max_grid) :: atom_pot_es
  !!!!!!!!!!!!!!!!!!!!!!!!
  ! VARIABLES I'VE ADDED !
  !!!!!!!!!!!!!!!!!!!!!!!!
  integer :: i_orbital_1, i_orbital_2, i_l, i_n, i_grid_1, i_grid_2, i_interp, i_interp2, i_index, i_spin
  real*8  :: spin_polarization_error
  logical :: spin_polarization_error_out
  ! The conversion between aims' orbital ordering and atom_sphere's.  See the section of the code 
  ! constructing this for more details.
  integer, allocatable :: as_to_aims_orb_order(:)
  real*8 :: log_coeff ! The logrithmic factor for r_grid, inferred from r_grid itself
  real*8 :: x, r_grid_0, temp1, temp2
  ! Whether to use a modified form of aim's r_grid or generate it from atom_sphere's loggrid and
  ! loggrid4 functions with their hardwired parameters
  logical, parameter :: use_atom_sphere_hardwired_grid = .false.
  !!!!!! NOTE !!!!!!!!
  ! We've found that the accuracy of the splined second derivative of the wavefunction
  ! depends on the number of points added to the beginning of the radial grid; the more 
  ! points we add, the closer the second derivative comes to the expected value.
  ! This is not a minor issue, either;  we need accurate second derivatives to
  ! accurately calculate kinetic, as close to the nucleus the (non-relativistic)
  ! kinetic term consists of two large numbers that mostly cancel one another.
  ! My current belief is that this is due to the shift of the radial grid needs by
  ! atom_sphere in order that the first point lies on the origin;  the more points we add 
  ! to the beginning of the radial grid, the smaller that shift will be, and the more the 
  ! atom_sphere rr radial grid resembles aims' r_grid, reducing interpolation error.
  ! 
  ! 100 pre-pended points works, but this is far too many points and may cause convergence
  ! issues in some cases in atom_sphere (according to the authors)!  Something needs to be
  ! done to avoid this interpolation in the first place, likely a rewrite of atom_sphere to 
  ! support aims' radial grid)
  integer :: n_prepend_points = 100, n_postpend_points = 3

  character*300 :: info_str
  ! The following are quantities constructed on the atom_sphere rr grid, but using the appropriate
  ! re-definition to make them compatible with aims' expected form (e.g. wave = rr * psi and so on )
  real*8, allocatable :: wave_as(:,:,:,:)
  real*8, allocatable :: kinetic_as(:,:,:,:)
  real*8, allocatable :: density_deriv_as(:,:)
  real*8, allocatable :: density_2nd_deriv_as(:,:)
  real*8, allocatable :: atom_pot_es_as(:)
  real*8, allocatable :: atom_pot_hart_as(:)
  real*8, allocatable :: atom_pot_xc_as(:,:)
  real*8, allocatable :: atom_pot_as(:,:)
  ! Variables related to interpolation of quantities onto aim's r_grid
  real*8, allocatable :: wave_test(:,:)
  real*8, allocatable :: wave_deriv_test(:,:)
  real*8, allocatable :: wave_2nd_deriv_test(:,:)
  real*8, allocatable :: kinetic_test(:,:)
  real*8, allocatable :: density_test(:,:)
  real*8, allocatable :: density_deriv_test(:,:)
  real*8, allocatable :: density_2nd_deriv_test(:,:)
  real*8, allocatable :: atom_pot_es_test(:)
  real*8, allocatable :: atom_pot_hart_test(:)
  real*8, allocatable :: atom_pot_xc_test(:,:)
  real*8, allocatable :: atom_pot_test(:,:)
  real*8, allocatable :: spl_param(:,:)
  real*8 :: r_output
  real*8 :: convert_radial_grid_to_as_log_grid ! external function, see bottom of file
  real*8 :: didr, di2dr2
  real*8, parameter   :: TOLERANCE = 1.0d-6

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ARGUMENTS FOR ATOMSOLVER !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The following are the arguments for atom_sphere's atomsolver subroutine.
  ! I have deliberately kept their naming convention, F77-y as it is, for easy
  ! comparison to the original source code.

  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! HARD CODED PARAMETERS !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: lmax = 3
  integer, parameter :: lmaxp = 2*lmax
  integer, parameter :: nprinx = 7
  ! integer, parameter :: nrad = 532
  integer :: nrad ! This WAS hard-coded to be 532; here it can either be set to 532 or based upon aims' grid
  integer, parameter :: norbmax = 60
  integer, parameter :: idsx = 1
  integer, parameter :: nspol = 2
  character*1 :: il(5)
  integer, parameter :: ncovmax = 5
  real*8 :: epstres
  logical :: pspcalcul
  logical :: kinonly

  !!!!!!!!!!!!!!!!!!!!!!!!
  ! SPECIFIED IN INPUT() !
  !!!!!!!!!!!!!!!!!!!!!!!!
  ! Input parameters originally read in by the subroutine input() from the file atom.dat,
  ! here provided by aims.
  ! Organized by order of appearance in atom.dat
  character*2 :: nameat
  integer :: iXC
  real*8 :: hfmix
  real*8 :: omega
  character*1 :: ispp ! "Treatment of the electrons: Do not change this (spin-polarized)."
  real*8 :: znuc  ! "Znuc of the atom"
  integer :: ncov ! "number of covalent radius on which the norm conservation is imposed. At present set it to 1."
  real*8 :: rcov(ncovmax) ! "the cut-off radius"
  real*8, parameter :: rprb =0.0d0 ! " a parameter for  the external parabolic potential to be applied. For AE cases 
                                   !   set it to 0.d0"
                                   ! Since we have modified atom_sphere to use aims' confining potential, turn the
                                   ! external parabolic potential off
  integer :: lo(norbmax)
  integer :: norb
  integer :: noae(norbmax)
  real*8 :: so(norbmax)
  real*8 :: zo(norbmax)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! INFERRED FROM USER INPUT !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical :: screened
  real*8, allocatable :: screenexkernel(:,:,:)   ! represents the action of the screened Coulomb operator
  integer :: nspin
  integer :: no(norbmax)
  integer :: nprin(lmax+1,nspol)    ! Store n,l,s of pseudo orbitals

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ATOM-SPECIFIC INPUT TO ATOM_SPHERE !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real*8, allocatable :: psi(:,:,:,:) ! radial wavefunctions (input and output of WF optimizationpsi
  real*8, allocatable :: as_kinetic(:,:,:,:) ! radial wavefunctions (input and output of WF optimizationpsi
  real*8 :: occup(nprinx,nspol,lmax+1) ! occupation number
  real*8 :: shift(nprinx,lmax+1) ! shifts used for preconditoning (CRITCAL FOR FAST CONVERGENCE!)

  !!!!!!!!!!!!!!!!!!!!!!
  ! RADIAL GRID ARRAYS !
  !!!!!!!!!!!!!!!!!!!!!!
  real*8, allocatable :: rr(:)    !radial grid
  real*8, allocatable :: rw(:,:)    ! radial weights rw(*,1) pure widths, rw(*,2) multiplied with r^2, rw(*,3) inverse radial weight for GGA's
  real*8, allocatable :: rr4(:) !denser version for kinetic energy
  real*8, allocatable :: rw4(:,:)  !denser version for kinetic energy

  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! OUTPUT OF ATOM_SPHERE !
  !!!!!!!!!!!!!!!!!!!!!!!!!
  real*8 :: eval(nprinx,nspol,lmax+1) ! eigenvalues of the Kohn Sham Hamiltonian (Fock matrix)
  real*8 :: chrg(nprinx,nspol,lmax+1,ncovmax) ! Charge of  the Kohn Sham Orbitals
  real*8 :: residue(nprinx,nspol,lmax+1)  ! residue of KS Orbitals
  real*8 :: ekin_orb(nprinx,nspol,lmax+1)  ! Kinetic Energy of KS Orbitals
  real*8 :: epot_orb(nprinx,nspol,lmax+1)  ! Potential Energy of KS Orbitals for znuc 
  real*8, allocatable :: potloc(:)    ! local potential of pseudopotential
  real*8, allocatable :: pothart(:)   ! Hartree potential
  real*8, allocatable :: potxc(:,:)  ! exchange correlation potential
  real*8, allocatable :: pothf(:,:) ! Hartree Fock potential
  real*8, allocatable :: pot0(:,:,:)  ! potential used for preconditioning
  real*8 :: ccleb(lmaxp+1,lmax+1,lmax+1) !contracted (overl m) Clebsch coefficients
  real*8, allocatable :: grad(:,:,:,:) ! total gradient
  real*8, allocatable :: gradp(:,:,:,:,:)  ! histopry [sic] list of preconditioner gradients for DIIS
  real*8, allocatable :: grads(:,:,:,:)  ! gradient of overlap
  real*8, allocatable :: gradhf(:,:,:,:) ! gradient form Hartree Fock
  real*8, allocatable :: hhp(:,:,:,:) ! Hamiltonian matrice for preconditioning
  real*8, allocatable :: ssp(:,:) ! overlap matrice for preconditioning
  real*8 :: epslag(nprinx,nprinx,nspol,lmax+1)  ! Lagrange multiplier matrix
  real*8, allocatable :: rhotot(:)    ! total charge density
  real*8, allocatable :: rhoud(:,:) ! spinup, down charge density
  real*8, allocatable :: rhoij(:)   ! cross charge densities for Hartree Fock
  real*8, allocatable :: wrky(:,:)  !Work array
  real*8 :: adiis(idsx+1,idsx+1,3),ipiv(idsx+1),rdiis(idsx+1) ! work arrays for DIIS
  real*8, allocatable :: ppois(:,:,:)  !auxiliary array for solving Poisson's equation
  real*8 :: etotal
  real*8, allocatable :: psid(:,:,:,:,:)   !WF history  used uniquely in the WF optimization part
  real*8 :: time(5)

  !!!!!!!!!!!!!!!!!!!!!!
  ! INTERNAL VARIABLES !
  !!!!!!!!!!!!!!!!!!!!!!
  ! (e.g., all variable not entering into atomsolver)
  character(80) :: errmsg
  character*200 :: label
  integer :: nomax(lmax+1) 
  integer :: nomin(lmax+1) 
  ! The following variables were implicit (this doesn't include variables passed
  !    into atomsolver which were defined implicitly)
  ! Fun fact: the FORTRAN 77 standard is that implicit variables beginning with
  !           letters i-n are integers, and all others are real
  !           The authors changed it from real to real*8, though 
  real*8 :: t0
  real*8 :: fourpi
  integer :: lmaxtemp
  integer :: lmaxoccup
  integer :: iorb
  integer :: ierr
  integer :: l
  integer :: nprinxtemp
  integer :: i
  integer :: np
  integer :: isp
  integer :: iprin
  real*8 :: fcr
  real*8 :: zeff
  integer :: ll
  real*8 :: si
  real*8 :: oup
  real*8 :: odown
  real*8 :: ss
  integer :: irad
  integer :: nn
  real*8 :: zz
  integer :: j
  real*8 :: t

  write(info_str, '(A)')
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  "************************************"
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  "* STARTING ATOM_SPHERE CALCULATION *"
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  "************************************"
  call localorb_info( info_str )
  write(info_str, '(A)')
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  "          KNOWN BUGS/ISSUES         "
  call localorb_info( info_str )
  write(info_str, '(A)')
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  " * Does not support relativistic corrections"
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  " - This is a limitation of the present status of integration into aims.  Programmer"
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  "   time, mostly."
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  " * Only works for atoms with Z=1 to Z=36"
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  " - This is a limitation in atom_sphere itself.  Looking at the code, I (WPH)"
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  "   am guessing that it has to do with their angular integration grids."
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  " * Confined basis elements are not supported."
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  " - This limitation arises from these basis elements requiring free_potential be"
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  "   valid, which will not be true for hybrid/HF calculations.  However, these basis"
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  "   elements only show up in the standard tiered basis set for f-shell metals, which"
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  "   A) have Z > 36 and B) are avoided like the plague in DFT, so it's a minor limitation."
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  " * Iterative Hirshfeld analysis is not supported."
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  " - I don't understand how this subroutine works, I don't have the time to understand"
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  "   it, and they're not paying me to understand it.  You're on your own, kid!"
  call localorb_info( info_str )
  write(info_str, '(A)')  ""
  call localorb_info( info_str )

  if ( n_l_channels.gt.1 ) then
    call aims_stop("Error:  the atom_sphere wrapper does not support outputting &
         &l-channel-dependent potentials.  Exiting.")
  end if

  spin_polarization_error_out = .false.
  ! Set up the grids based on whether we want to use aims' r_grid (with suitable
  ! shift and additional appended points) or the atom_sphere hardwired grid
  if (use_atom_sphere_hardwired_grid) then
    nrad = 532
  else
    nrad = n_grid + n_prepend_points + n_postpend_points
  end if

  allocate( screenexkernel(nrad,nrad,lmaxp+1) )   ! represents the action of the screened Coulomb operator
  allocate( psi(nrad,nspol,nprinx,lmax+1) ) ! radial wavefunctions (input and output of WF optimizationpsi
  allocate( as_kinetic(nrad,nspol,nprinx,lmax+1) ) ! radial wavefunctions (input and output of WF optimizationpsi
  allocate( rr(nrad) )    !radial grid
  allocate( rw(nrad,3) )    ! radial weights rw(*,1) pure wigths, rw(*,2) mulitplied with r^2, rw(*,3) inverse radial weight for GGA's
  allocate( rr4(4*nrad) ) !denser version for kinetic energy
  allocate( rw4(4*nrad,3) )  !denser version for kinetic energy
  allocate( potloc(nrad) )    ! local potential of pseudopotential
  allocate( pothart(nrad) )   ! Hartree potential
  allocate( potxc(nrad,nspol) )  ! exchange correlation potential
  allocate( pothf(nrad,lmaxp+1) ) ! Hartree Fock potential
  allocate( pot0(nrad,nprinx,lmax+1) )  ! potential used for preconditioning
  allocate( grad(nrad,nspol,nprinx,lmax+1) )  ! total gradient
  allocate( gradp(nrad,nspol,nprinx,lmax+1,idsx) )  ! histopry [sic] list of preconditioner gradients for DIIS
  allocate( grads(nrad,nspol,nprinx,lmax+1) )  ! gradient of overlap
  allocate( gradhf(nrad,nspol,nprinx,lmax+1) ) ! gradient form Hartree Fock
  allocate( hhp(2,nrad,nprinx,lmax+1) ) ! Hamiltonian matrice for preconditioning
  allocate( ssp(2,nrad) ) ! overlap matrice for preconditioning
  allocate( rhotot(nrad) )    ! total charge density
  allocate( rhoud(nrad,nspol) ) ! spinup, down charge density
  allocate( rhoij(nrad) )   ! cross charge densities for Hartree Fock
  allocate( wrky(nrad,nprinx) )  !Work array
  allocate( ppois(nrad,-(lmaxp+1):lmaxp,2) )  !auxiliary array for solving Poisson's equation
  allocate( psid(nrad,nspol,nprinx,lmax+1,idsx) )   !WF history  used uniquely in the WF optimization part

  allocate( as_to_aims_orb_order(2*n_valence) )

  allocate( wave_as(nrad, nspol, nprinx, lmax+1) )
  allocate( kinetic_as(nrad, nspol, nprinx, lmax+1) )
  allocate( density_deriv_as(nrad, n_spin) )
  allocate( density_2nd_deriv_as(nrad, n_spin) )
  allocate( atom_pot_es_as(nrad) )
  allocate( atom_pot_hart_as(nrad) )
  allocate( atom_pot_xc_as(nrad, n_spin) )
  allocate( atom_pot_as(nrad, n_spin) )

  ! Note that wave and its derivative quantities do not have an
  ! n_spin index, as these quantities implicitly include
  ! the spin channel in n_max_ind_fns
  allocate( wave_test(n_max_grid, n_max_ind_fns) )
  allocate( wave_deriv_test(n_max_grid, n_max_ind_fns) )
  allocate( wave_2nd_deriv_test(n_max_grid, n_max_ind_fns) )
  allocate( kinetic_test(n_max_grid, n_max_ind_fns) )
  allocate( density_test(n_grid, n_spin) )
  allocate( density_deriv_test(n_grid, n_spin) )
  allocate( density_2nd_deriv_test(n_grid, n_spin) )
  allocate( atom_pot_es_test(n_grid) )
  allocate( atom_pot_hart_test(n_grid) )
  allocate( atom_pot_xc_test(n_grid, n_spin) )
  allocate( atom_pot_test(n_grid, n_spin) )

  allocate( spl_param(4,nrad) )

  if (use_atom_sphere_hardwired_grid) then
    !!write(info_str,'(2X, A)') "Running atom_sphere with the atom_sphere hardwirded grid."
    call localorb_info( info_str )

    log_coeff = log(1.05d0) 
    r_grid_0 = 1.0d-3    
    call loggrid (nrad,r_grid_0, exp(log_coeff), rr, rw)
    call loggrid4(nrad,r_grid_0, exp(log_coeff), rr4, rw4) !! REFINED grid FOR K.E.
  else
    ! BEGIN AIM'S GRID !
    write(info_str,'(2X,A)') "Running atom_sphere with aims' (modified) r_grid."
    call localorb_info( info_str )
    write(info_str,'(2X,A,I5)') "Number prepended points:     ", n_prepend_points
    call localorb_info( info_str )
    write(info_str,'(2X,A,I5)') "Number postpended points:    ", n_postpend_points
    call localorb_info( info_str )

    ! We use the weight modifications provided with the original atom_sphere code, which are performed on the
    ! first 4 points and the last 4 point of the radial grid
    if (nrad .lt. 8) then
      call aims_stop("Error detected in atom_sphere_wrapper:  the radial grid provided has less than 8 points! &
                     &Either your input files are badly misformed, or you hit a bug in the code. Exiting.",'')
    end if

    ! To convert r_grid into atom_sphere's rr, we need knowledge of
    ! the grid construction; not only do we need its values, but we also 
    ! need the derivatives and a way to interpolate consistently between 
    ! grid points 
    ! The following code presumes that the aims radial grid is a logarithmetic 
    ! grid set up by
    !     r_grid(i) = r0*exp(A*(i-1)) = exp(A) * r_grid(i-1)
    ! which, at this time in aims, is a pretty safe bet.  In this case, the 
    ! derivatives/radial weights become
    !     rw(i) = dr/di = A*r_grid(i)
    ! We need the logarithmetic coefficient A, which can be inferred as
    log_coeff = log(r_grid(2)/r_grid(1))
    ! This does have a variable associated with it within aims, but this makes
    ! the math more transparent and cuts down on interface cruft.

    ! atom_sphere generates its radial grid by using the aforementioned
    ! logarithmetic grid, similar to aims, but it then shifts the grid back so that
    ! the first element is the origin
    !   rr(i) = r_grid(i) - r_grid(1)
    ! which, defining r_grid(1) = r_grid_0 and taking into account the points
    ! prepended to aims' radial grid, will be
    r_grid_0 = r_grid(1) * exp( -1.0d0*n_prepend_points*log_coeff )
    ! Once we have this, we have everything we need to modify aims' r_grid to
    ! fit with atom_sphere's rr

    ! Add the points at the beginning of the radial grid, if desired
    do i_grid_2 = 1, n_prepend_points
      rr(i_grid_2)   = r_grid_0 * exp( log_coeff*(i_grid_2-1) ) - r_grid_0
      rw(i_grid_2,1) = (rr(i_grid_2) + r_grid_0) * log_coeff
      rw(i_grid_2,2) = rr(i_grid_2)**2 * rw(i_grid_2,1)
      rw(i_grid_2,3) = 1.0d0 / rw(i_grid_2,1)
    end do

    ! Add points from r_grid to rr and apply the shift
    do i_grid_2 = n_prepend_points+1, nrad-n_postpend_points
      rr(i_grid_2) = r_grid( i_grid_2-n_prepend_points ) - r_grid_0
      rw(i_grid_2,1)  = r_grid( i_grid_2-n_prepend_points ) * log_coeff
      rw(i_grid_2,2)  = rr(i_grid_2)**2 * rw(i_grid_2,1)
      rw(i_grid_2,3)  = 1.0d0 / rw(i_grid_2,1)

      ! As an error check, make sure that r_grid really is a logarithmetic
      ! grid by recomputing the logarithmetic coefficient A at each point 
      ! and ensuring that it matches the reference value log_coeff
      if (i_grid_2 .gt. n_prepend_points+1) then
        if (abs(1.0d0 - log(r_grid(i_grid_2-n_prepend_points)/r_grid(i_grid_2-n_prepend_points-1))/log_coeff ) &
             .gt. TOLERANCE) then
          call aims_stop("Error detected in atom_sphere_wrapper:  radial grid provided was not &
               &logarithmetic!  This should not happen in aims (unless a new kind of grid was added &
               &after this function was written.)  Exiting.",'')
        end if
      end if
    end do

    ! Add points to the end of the radial grid
    do i_grid_2 = nrad-n_postpend_points+1, nrad
      rr(i_grid_2)   = r_grid_0 * exp( log_coeff*(i_grid_2-1) ) - r_grid_0
      rw(i_grid_2,1) = (rr(i_grid_2) + r_grid_0) * log_coeff
      rw(i_grid_2,2) = rr(i_grid_2)**2 * rw(i_grid_2,1)
      rw(i_grid_2,3) = 1.0d0 / rw(i_grid_2,1)
    end do

    ! Construct rr4 and rw4, the radial grid interpolated by a factor of 4
    ! for evaluation of the kinetic energy in atom_sphere
    do i_grid_1 = 1, nrad
      do i_grid_2 = 1, 4
        ! This requires we use the assumed logarithmetic form for the radial grid
        i_index = 4*(i_grid_1-1) + i_grid_2
        x = (i_grid_1-1)+0.25d0*(i_grid_2-1)
        rr4(i_index) = r_grid_0 * exp( log_coeff*x ) - r_grid_0
        rw4(i_index,1) = 0.25d0 * (rr4(i_index) + r_grid_0) * log_coeff ! grid is four times denser
        rw4(i_index,2) = rr4(i_index)**2 * rw4(i_index,1)
        rw4(i_index,3) = 0.25d0 / rw4(i_index,1)
      end do
    end do

    ! For reasons I don't yet understand, atom_sphere needs to have the weights
    ! altered at the end point
    rw(1,2)=rw(1,2)*17.d0/48.d0
    rw(2,2)=rw(2,2)*59.d0/48.d0
    rw(3,2)=rw(3,2)*43.d0/48.d0
    rw(4,2)=rw(4,2)*49.d0/48.d0
    rw(nrad-0,2)=rw(nrad-0,2)*17.d0/48.d0
    rw(nrad-1,2)=rw(nrad-1,2)*59.d0/48.d0
    rw(nrad-2,2)=rw(nrad-2,2)*43.d0/48.d0
    rw(nrad-3,2)=rw(nrad-3,2)*49.d0/48.d0
    rw4(1,2)=rw4(1,2)*17.d0/48.d0
    rw4(2,2)=rw4(2,2)*59.d0/48.d0
    rw4(3,2)=rw4(3,2)*43.d0/48.d0
    rw4(4,2)=rw4(4,2)*49.d0/48.d0
    rw4(4*nrad-0,2)=rw4(4*nrad-0,2)*17.d0/48.d0
    rw4(4*nrad-1,2)=rw4(4*nrad-1,2)*59.d0/48.d0
    rw4(4*nrad-2,2)=rw4(4*nrad-2,2)*43.d0/48.d0
    rw4(4*nrad-3,2)=rw4(4*nrad-3,2)*49.d0/48.d0
    ! END AIM'S GRID !
  end if

  write(info_str,'(2X,A,F25.17)') "r_grid_0:                    ", r_grid_0
  call localorb_info( info_str )
  write(info_str,'(2X,A,F15.8)') "log_coeff:                   ", log_coeff
  call localorb_info( info_str )
  write(info_str,'(2X,A,F25.17)') "r_grid(1):                   ", r_grid(1) 
  call localorb_info( info_str )
  write(info_str,'(2X,A,F25.17)') "r_grid(2):                   ", r_grid(2) 
  call localorb_info( info_str )
  write(info_str,'(2X,A,F15.8)') "First non-origin grid point: ", rr(2)
  call localorb_info( info_str )
  write(info_str,'(2X,A,F15.8)') "r_grid(n_grid):              ", r_grid(n_grid) 
  call localorb_info( info_str )
  write(info_str,'(2X,A,F15.8)') "Last grid point:             ", rr(nrad)
  call localorb_info( info_str )
 
  epstres=1.d-12  ! convergence threshold
  pspcalcul=.false.
  kinonly=.false.
  time=0.d0
  call cpu_time(t0)

  il(1) = 's'
  il(2) = 'p'
  il(3) = 'd'
  il(4) = 'f'
  il(5) = 'g'

  fourpi=16.d0*atan(1.d0)
  screened = .false. 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The calculation now has set up everything it needs to continue      !
  ! So now we perform reverse transcription and inject our RNA into the !
  ! code.  To understand what is going on, look at input.f90, the       !
  ! original subroutine for reading in atom.dat for atom_sphere         !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nameat = species_name
  znuc = nuclear_charge
  ! Set up details of the XC functional
  hfmix = 0.0d0
  omega = 0.0d0

  ! This nomenclature for specifying the functional used comes from a combination of 
  ! an interface layer between libxc and a third piece of code (ABINIT, possibly) that 
  ! was bolted onto atom_sphere and the ability of libxc to mix-and-match exchange and 
  ! correlation pieces of different ! functionals, requiring that certain functionals 
  ! set up both seperately.
  ! If you look at the libxc module libxc_func_m, where the xc functionals are defined by
  ! flags, this will make sense.  Some functional flags have XC in their name, corresponding 
  ! to both exchange and correlation combined, and some functionals only have X or 
  ! C in their name, signifying that you need to initalize both piece seperately.
    
  ! The first three numbers correspond to the exchange contribution to the
  ! functional.  The second three numbers correspond to the correlation
  ! contribution.  When the first three numbers are zero, that signifies
  ! that both exchange and correlation are combined, and the functional
  ! only needs to be initalized once.
  ! The negative sign is just there, because YOLO.

  ! LDAs
  if      (flag_xc .eq.  3) then
     ! PZ-LDA, added by WPH
     iXC=-001009
  else if (flag_xc .eq.  8) then
     ! PW-LDA, added by WPH
     iXC=-001012
  else if (flag_xc .eq. 15) then
     ! VWN5
     iXC=-001007
  else if (flag_xc .eq. 16) then
     ! VWN (RPA), added by WPH
     iXC=-001008
  ! GGAs
  else if (flag_xc .eq.  6) then  
     ! PBE 
     iXC=-101130
  else if (flag_xc .eq. 17) then
     ! PBEsol
     iXC=-116133 
  else if (flag_xc .eq. 22) then
     ! xPBE
     iXC=-123136
  ! Hybrids
  else if (flag_xc .eq.  1) then
     ! PBE0 
     iXC=-000406
     hfmix = hybrid_coeff
  else if (flag_xc .eq. 10) then  
     ! B3LYP 
     iXC=-000402
     hfmix = hybrid_coeff
  else if (flag_xc .eq.  7) then  
     ! HSE06
     ! NOTE:  In conversations with the developer of atom_sphere, it may be
     ! necessary to increase the cutoff radius of the screened exchange kernel
     ! constructed in create_screenexkernel by uncommenting the creation of
     ! "call  radgrid(nradxx,rr,rw)" in that function and increasing nradxx so that
     ! the auxillary grid goes out at least 50 a.u. past the main grid
     ! (i.e. nradxx must be larger than nrad,) because otherwise some 
     ! of the convolved charge will be lost and the ES potentials will be bad.
     iXC=-000428
     hfmix = hybrid_coeff
     omega = hse_omega
  else if (flag_xc .eq.  0) then
     ! Note:  When running with Hartree-Fock, atom_sphere never calls libxc
     ! functions, and you can run the code without having linked aginst
     ! libxc.  
     ! Not sure how useful that could be, but hey, it's there.
     ! Hartree-Fock
     iXC=100
     hfmix = 1.0d0 
  else if (flag_xc .le. -10) then
     ! We are already using LibXC for our functional choice in FHI-aims, 
     ! so lets set iXC appropriately
     iXC = flag_xc/10
  else
    write(info_str,'(2X,A)') "Error in atom_sphere_wrapper:  XC selected is currently not supported."
    call localorb_info( info_str )
    write(info_str,'(2X,A)') &
        "For many functionals, adding support is as simple as adding two lines into the atom_sphere wrapper!"
    call localorb_info( info_str )
    call aims_stop("Error detected in atom_sphere_wrapper:  XC selected is currently not supported.",'')
  end if
  ispp = 's'   ! Non-relativistic, spin-polarized calculations.  Contrary to what the comments in 
               ! atom_sphere state (which I have preserved, just in case,) this is NOT a relativistic 
               ! calculation.  
               ! Because aims uses non-spin-polarized basis elements (even for spin-polarized calculations,) 
               ! we equally occupy both spin channels of the atom to make the calculation
               ! non-spin-polarized.
  ! ispp = 'n' ! Non-relativistic, non-spin-polarized calculations.  Not implemented in atom_sphere.
  ! ispp = 'r' ! Relativistic, non-spin-polarized calculations.  Not implemented in atom_sphere.
  ! The following has to do with pseudopotentials; we turn this off for all-electron calculations
  ncov = 1
  rcov = 2.23d0

  ! This part of input.f90 is so utterly incomprehensible that I had to guess 
  ! what is going on by printing out variables and mapping them onto aims'.
  ! Given some of the comments in that code, whoever was in charge was also 
  ! confused by it.

  ! aims sorts orbitals first by the l quantum number, then by n, whereas 
  ! atom_sphere does the exact opposite.  So we have to resort them.  On
  ! top of that, we need to account for aims' different spin indexing
  ! of the orbitals; aims tacks on the down channel (should it be included)
  ! at the end of the array, whereas atom_sphere interweaves the two
  ! channels.  The two cases are sufficiently different that a hard-coded fork in
  ! spin channels is needed, inelegant as it is.
  i_l = 0
  i_n = 1
  norb = 0

  if (n_spin.eq.1) then
    do i_orbital_1 = 1, 2*n_valence   ! indexes atom_aphere orbitals (both up and down)
      do i_orbital_2 = 1, n_valence   ! indexes aims orbitals
        if (n_shell(i_orbital_2) .eq. i_n .and. l_shell(i_orbital_2) .eq. i_l) then
          norb = norb + 2
          noae(norb-1)   =  n_shell(i_orbital_2) 
          lo(norb-1)     =  l_shell(i_orbital_2)
          noae(norb)     =  n_shell(i_orbital_2)
          lo(norb)       =  l_shell(i_orbital_2)

          if (ispp .eq. 's') then
            so(norb-1)     =  0.5d0
            zo(norb-1)     =  occup_shell(i_orbital_2) / 2.0d0
            so(norb)       = -0.5d0
            zo(norb)       =  occup_shell(i_orbital_2) / 2.0d0
          else
            call aims_stop('You are attempting to run atom_sphere with an unsupported spin treatment.  Exiting.')
          end if
          ! Don't ask me why atom_sphere does this, it just does
          if ( zo(norb-1) .eq. 0.0d0 ) zo(norb-1) = 1.0d-20
          if ( zo(norb) .eq. 0.0d0 )   zo(norb)   = 1.0d-20
          ! Save the ordering, since we'll need to reverse it when exporting
          ! the results of atom_sphere back out to aims
          as_to_aims_orb_order(norb-1) = i_orbital_2
          as_to_aims_orb_order(norb)   = i_orbital_2
          ! Since we've found the desired orbital, go to the next l-orbital in
          ! this n-shell and reset the counter
          exit
        end if
      end do
          
      i_l = i_l + 1

      if (i_l .ge. i_n) then
        ! We've found all orbitals for this n-shell, move on to the next
        ! n-shell
        i_n = i_n + 1
        i_l = 0
      end if
    end do 
  else
    do i_orbital_1 = 1, 2*n_valence     ! indexes atom_aphere orbitals (both up and down, but interweaved)
      do i_orbital_2 = 1, n_valence     ! indexes aims orbitals, only the spin up channel
        if (n_shell(i_orbital_2) .eq. i_n .and. l_shell(i_orbital_2) .eq. i_l) then
          ! The first n_max_ind_fns/2 states for aims should be spin-up, and the
          ! last n_max_ind_fns/2 states should have identical orbital spin numbers
          ! but be spin down
          if ( ( n_shell(i_orbital_2) .ne. n_shell(i_orbital_2+n_valence) ) .or. &
               ( l_shell(i_orbital_2) .ne. l_shell(i_orbital_2+n_valence) ) .or. &
               ( s_shell(i_orbital_2) .ne. 1 ) .or. ( s_shell(i_orbital_2+n_valence) .ne. 2 ) ) then
            call aims_stop('atom_sphere_wrapper did not find the expected layout of quantum numbers. Exiting.')
          end if

          norb = norb + 2
          noae(norb-1)   =  n_shell(i_orbital_2)
          lo(norb-1)     =  l_shell(i_orbital_2)
          noae(norb)     =  n_shell(i_orbital_2+n_valence)
          lo(norb)       =  l_shell(i_orbital_2+n_valence)

          if (ispp .eq. 's') then
            so(norb-1)     =  0.5d0
            zo(norb-1)     =  occup_shell(i_orbital_2)
            so(norb)       = -0.5d0
            zo(norb)       =  occup_shell(i_orbital_2+n_valence)
          else
            call aims_stop('You are attempting to run atom_sphere with an unsupported spin treatment.  Exiting.')
          end if
          ! Don't ask me why atom_sphere does this, it just does
          if ( zo(norb-1) .eq. 0.0d0 ) zo(norb-1) = 1.0d-20
          if ( zo(norb) .eq. 0.0d0 )   zo(norb)   = 1.0d-20
          ! Save the ordering, since we'll need to reverse it when exporting
          ! the results of atom_sphere back out to aims
          as_to_aims_orb_order(norb-1) = i_orbital_2
          as_to_aims_orb_order(norb)   = i_orbital_2+n_valence
          ! Since we've found the desired orbital, go to the next l-orbital in
          ! this n-shell and reset the counter
          exit
        end if
      end do
      
      i_l = i_l + 1

      if (i_l .ge. i_n) then
        ! We've found all orbitals for this n-shell, move on to the next
        ! n-shell
        i_n = i_n + 1
        i_l = 0
      end if
    end do 
  end if 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! From here until the final post-processing, this code is ripped from atomsphere.f90 !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! (though with some output lines removed, output funneled through localorb_info,
  !  and whitespace changed to the way I like it, because I can)

  lmaxtemp=0
  lmaxoccup=0
  do iorb=1,norb
    if (zo(iorb).gt.1.d-10)lmaxoccup=max(lo(iorb),lmaxoccup)
    lmaxtemp=max(lo(iorb),lmaxtemp)
  end do
  if(lmaxtemp.gt.lmax)ierr=3   !Stefan modified +1 -> -1
  write(errmsg,*)'array dimension problem:lmax'

  nomin=100
  nomax=0
  do iorb=1,norb
    nomin(lo(iorb)+1)=min(noae(iorb),nomin(lo(iorb)+1))
    nomax(lo(iorb)+1)=max(noae(iorb),nomax(lo(iorb)+1))
    no(iorb)=noae(iorb)
  enddo

  nprin=0
  do iorb=1,norb
    no(iorb)=noae(iorb)-nomin(lo(iorb)+1)+1
    if (so(iorb)==0.5d0) then
      nprin(lo(iorb)+1,1)=max(no(iorb),nprin(lo(iorb)+1,1))
    else if (so(iorb)==-0.5d0) then
      nprin(lo(iorb)+1,2)=max(no(iorb),nprin(lo(iorb)+1,2))
    end if
  end do

  do l=0,lmaxtemp
    nprinxtemp=(nomax(l+1)-nomin(l+1))+1
    if (nprinxtemp.gt.nprinx)ierr=3
    write(errmsg,*) 'array dimension problem: nprinxtemp'
  end do
  do i=1,ncov
    call detnp(nrad,rr,rcov(i),np)
    rcov(i)=rr(np)
  end do

  occup=0.d0 !! occup(nspol,nprinx,lmax+1)
  do iorb=1,norb
    if (so(iorb).gt.0.1d0)isp=1
    if (so(iorb).lt.0.1d0)isp=2
    l=lo(iorb)
    iprin=no(iorb)
    occup(iprin,isp,l+1)=zo(iorb)
    if (zo(iorb)==0.d0) then
      occup(iprin,isp,l+1)=1.d-20
    end if
  end do

!!  SETTING UP GUESS WAVE FUNCTION and SHIFT VALUES
  fcr=0.9d0
  zeff=znuc
  do iorb=1,norb,2
    ll=lo(iorb)
    iprin=no(iorb)
    si=so(iorb)
    isp=1
    if (si.lt.0.1d0)isp=2
    oup=occup(iprin,1,ll+1)
    odown=occup(iprin,2,ll+1)
    zeff=zeff-oup-odown
    if (zeff.lt.1.d0) then !! SS
      call random_number(ss) ! WPH:  Leads to a thread-non-safe situation!  This is necessary for
                             !       the initial guess for unoccupied states, however.
      !ss = 0.5              !       Will make the code (more) thread-safe, but at the risk of 
                             !       poor initial guess for unoccupied states.  
      zeff=ss*.5d0+.499d0 !! SG
    endif
    do irad=1,nrad
      psi(irad,1,iprin,ll+1)=exp(-fcr*zeff*rr(irad))*rr(irad)**ll
      psi(irad,2,iprin,ll+1)=exp(-fcr*zeff*rr(irad))*rr(irad)**ll
    end do
    shift(iprin,ll+1) = max((fcr*zeff)**2,1.d0)         
  end do

!!  SETTING UP LIBXC
  nspin=1
  if (ispp.eq.'s')nspin=2
  if (iXC.eq.100) then
    omega=0.d0;hfmix=1.d0 !! HF activation
    screened=.false.
  else
    call libxc_functionals_init(iXC,nspin,omega,hfmix)
  end if
  if (omega.gt.0.d0) then
    screened=.true.
    call create_screenexkernel(nrad,rr,rw,lmax,lmaxp,omega,screenexkernel)
  end if
  write(info_str,'(2X,A,F15.8,A,F15.8)')" Hfmix=",hfmix," Omega=",omega
  call localorb_info( info_str )

 ! Perform atom_sphere's SCF cycle
 ! (Rundong) atomsphere adopts an input guessed wavefunction psi, and thus generates charge density, energy terms (kinetic, Hartree, XC)
 ! at the beginning of the SCF loop. The information of the energy/potential terms are saved in the form of gradient, 
 ! as gradient (the array "grad") denotes: (the operater + potential) times psi.
 ! The grad for kinetic energy is generated through subroutine applykinpot. For ZORA case, this term is different (see subroutine applyzorakinpot).
 ! The Lagrange multiplier matrix seems to help to minimize the energy and thus to optimize the wavefunction. After that,
 ! the wavefunction is updated, and the SCF loop goes to the next iteration.
  if(flag_rel.eq.REL_atomic_zora)then
    call zorasolver(nameat,nprinx,nprin,nrad,lmaxp,lmax,norb,norbmax,idsx,nspin,nspol,occup,eval,chrg, &
       residue,ekin_orb,epot_orb,noae,no,il,lo,zo,so,ncovmax,ncov,potloc,pothart,potxc,pothf,pot0,ccleb,grad,&
       gradp,grads,gradhf,hhp,ssp,epslag,shift,rcov,rprb,znuc,psi,psid,as_kinetic,rr,rw,rr4,rw4,screened,&
       omega,hfmix,rhotot,rhoud,rhoij,screenexkernel,ppois,wrky,adiis,etotal,epstres,pspcalcul,kinonly,time,&
       cut_atom,cutoff_type,r_cut,scale_cutoff,w_cutoff)
  else
    call atomsolver(nameat,nprinx,nprin,nrad,lmaxp,lmax,norb,norbmax,idsx,nspin,nspol,occup,eval,chrg, &
       residue,ekin_orb,epot_orb,noae,no,il,lo,zo,so,ncovmax,ncov,potloc,pothart,potxc,pothf,pot0,ccleb,grad,&
       gradp,grads,gradhf,hhp,ssp,epslag,shift,rcov,rprb,znuc,psi,psid,as_kinetic,rr,rw,rr4,rw4,screened,&
       omega,hfmix,rhotot,rhoud,rhoij,screenexkernel,ppois,wrky,adiis,etotal,epstres,pspcalcul,kinonly,time,&
       cut_atom,cutoff_type,r_cut,scale_cutoff,w_cutoff)
  endif
 
  call libxc_functionals_end()
  call cpu_time(t)
  time(1)=t-t0
  write(info_str,'(A)')
  call localorb_info( info_str )
  write(info_str,'(2X,A)') ' Simple Time Profiling'
  call localorb_info( info_str )
  write(info_str,'(2X,A)') ' _____________________'
  call localorb_info( info_str )
  write(info_str,'(A)')
  call localorb_info( info_str )
  write(info_str,'(2X,A,F9.3,A)')' Atomsolver  overall runtime ',time(2),' seconds'
  call localorb_info( info_str )
  write(info_str,'(2X,A,F9.3,A)')' Atomsolver  actual  runtime ',time(2)-time(4),' seconds'
  call localorb_info( info_str )
  write(info_str,'(2X,A,F9.3,A)')' Libxc       overall runtime ',time(3),' seconds'
  call localorb_info( info_str )
  if (hfmix.gt.0.d0) then
    write(info_str,'(2X,A,F9.3,A)')' Hatree-Fock overall runtime ',time(5),' seconds'
    call localorb_info( info_str )
  end if
  write(info_str,'(2X,A)')' ___________________________________________________________'
  call localorb_info( info_str )
  write(info_str,'(2X,A,F9.3,A)')'                    CPU time ',t-t0,' seconds'
  call localorb_info( info_str )
  write(info_str,'(2X,A)')' ___________________________________________________________'
  call localorb_info( info_str )
  write(info_str,'(2X,A)')'                                              finished'
  call localorb_info( info_str )
  write(info_str,'(A)')
  call localorb_info( info_str )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! atom_sphere has finished,  Now we convert its results into aims format, and !
  ! export the basis functions and derived quantities back out to aims          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! First we construct those quantites that are not stored by atom_sphere, but will
  ! be needed by aims, on atom_sphere's rr
  ! We only evaluate the functions on grid points where they are defined
  do i_grid_2 = 1, nrad
    if (i_grid_2 .ne. 1) then
      atom_pot_hart_as(i_grid_2) = pothart(i_grid_2)/rw(i_grid_2,2)
      atom_pot_es_as(i_grid_2) =  -1.0d0*nuclear_charge/rr(i_grid_2) + pothart(i_grid_2)/rw(i_grid_2,2)
    end if
  end do

  do i_spin = 1, n_spin
    do i_grid_2 = 1, nrad
      if (i_grid_2 .ne. 1) then
        ! We currently only use the semi-local part of the xc potential (at least that's what I think
        ! potxc is, anyways.)  This will give the wrong values for hybrids and HF!
        atom_pot_xc_as(i_grid_2,i_spin)   = potxc( i_grid_2, i_spin )/rw( i_grid_2, 2 ) 
        atom_pot_as(i_grid_2,i_spin) = atom_pot_xc_as(i_grid_2,i_spin) + atom_pot_es_as(i_grid_2)
      end if
    end do
  end do

  ! We also check to make sure spontaneous spin polarization hasn't occured,
  ! since we want spin-non-polarized basis functions but are using a
  ! spin-polarized atomic solver
  psi_at_0(:n_valence) = 0d0
  do i_orbital_1 = 1, norb
    ss=so(i_orbital_1)
    nn=no(i_orbital_1)
    ll=lo(i_orbital_1)
    zz=zo(i_orbital_1)
    if (ss.gt.0.d0)isp=1
    if (ss.lt.0.d0)isp=2

    ! This check will be done twice.  Oh well.
    ! This check is not perfect; for example, one spin channel could go to
    ! zero but the other spin channel approaches a finite number.  It should
    ! be rewritten!
    if (n_spin.eq.1) then
      do i_grid_2 = 1, nrad
        if ( psi(i_grid_2,1,nn,ll+1).gt.TOLERANCE .and. psi(i_grid_2,2,nn,ll+1).gt.TOLERANCE ) then
          temp1 = abs(1.0d0-psi(i_grid_2,2,nn,ll+1)/psi(i_grid_2,1,nn,ll+1))
          if ( temp1.gt.TOLERANCE ) then
            spin_polarization_error = max(temp1, spin_polarization_error)
            spin_polarization_error_out = .true.
          end if
        end if
      end do
    end if

    wave_as(1,isp,nn,ll+1) = rr(1)*( psi(1,isp,nn,ll+1) )
    do i_grid_2 = 2, nrad
      wave_as(i_grid_2,isp,nn,ll+1) = rr(i_grid_2)*( psi(i_grid_2,isp,nn,ll+1) )
      kinetic_as(i_grid_2,isp,nn,ll+1) = rr(i_grid_2)*( as_kinetic(i_grid_2,isp,nn,ll+1) )/rw(i_grid_2,2)
    end do
    kinetic_as(1,isp,nn,ll+1) = kinetic_as(2,isp,nn,ll+1)

    ! Only the s orbitals may have nonzero values at the nucleus.
    if (ll == 0) then
       i_index = as_to_aims_orb_order(i_orbital_1)
       ! Because psi_at_0 is non-spin-polarized, whenever i_index refers
       ! to spin channel 2 (i_index > n_valence), we map it back to spin
       ! channel 1.
       if (i_index > n_valence) i_index = i_index - n_valence
       ! Divide by 2 because we average over the spin channels
       psi_at_0(i_index) = psi_at_0(i_index) + psi(1,isp,nn,ll+1)/2
    end if
  end do

  if (spin_polarization_error_out) then
    write(info_str, '(2X,A,F15.8)') "Error in atom_sphere_wrapper:  maximum spin polarization of ", spin_polarization_error
    call localorb_info( info_str ) 
    write(info_str, '(2X,A,F15.8)') "has arisen in basis functions, exceeding the tolerance of ", TOLERANCE
    call localorb_info( info_str ) 
    write(info_str, '(2X,A)')       "The code currently averages the wavefunctions over spin channels, so we will"
    call localorb_info( info_str ) 
    write(info_str, '(2X,A)')       "continue as before.  There is the possibility that this is a numerical artifact,"
    call localorb_info( info_str ) 
    write(info_str, '(2X,A)')       "in which case this message can be ignored.  However, if this spin polarization is"
    call localorb_info( info_str ) 
    write(info_str, '(2X,A)')       "significant, there is something wrong with your calculation!"
    call localorb_info( info_str ) 
  end if

  ! Now spline the functions on atom_sphere's rr, so that we may interpolate
  ! them on aims' r_grid

  ! Density and its derivatives
  do i_spin = 1, n_spin 
    if (n_spin.eq. 1) then
      call cubic_spline( rhotot, nrad, spl_param )
    else
      call cubic_spline( rhoud( 1, i_spin ), nrad, spl_param )
    end if

    do i_grid_1 = 1, n_grid
      ! The factor of fourpi is a simple conversion factor between aims and atomsphere
      density_test(i_grid_1,i_spin) = fourpi*val_spline( convert_radial_grid_to_as_log_grid(r_grid(i_grid_1), &
                                        r_grid_0, log_coeff), spl_param, nrad )
      ! val_spline_deriv returns the derivative with respect to the logarithmetic
      ! grid coordinate, that is, the i entering into
      !     r(i) = r0*exp(A*(i-1)) - r0
      ! (remember that even though we are using splines to interpolate the function on points
      !  on aims' radial grid, we did the spline on atom_sphere's logarithmetic grid, so that 
      !  is why we are using its defining equation here)
      ! But what we want is the derivative with respect to the radial grid coordinate, that is,
      ! r_grid.  This conversion is a straight forward application of the chain rule, with the 
      ! above formula yielding 
      !     di/dr = (dr/di)^-1 = (A*(r+r_0))^-1
      didr = 1.0d0/(log_coeff*(r_grid(i_grid_1) + r_grid_0))
      density_deriv_test(i_grid_1,i_spin) = fourpi * didr * &
           val_spline_deriv( convert_radial_grid_to_as_log_grid(r_grid(i_grid_1), r_grid_0, log_coeff), spl_param, nrad )
      ! More chain rule, I don't feel like typing it up, seriously it's like
      ! freshman-level math
      ! That comment will be really embarassing when someone discovers that I did
      ! it wrong, won't it be?
      di2dr2 = -1.0d0 / (log_coeff*(r_grid(i_grid_1) + r_grid_0)**2.0d0)
      density_2nd_deriv_test(i_grid_1,i_spin) = fourpi * di2dr2 * &
           val_spline_deriv( convert_radial_grid_to_as_log_grid(r_grid(i_grid_1),  r_grid_0, log_coeff), spl_param, nrad ) + &
           fourpi * didr**2 * &
           val_spline_2nd_deriv( convert_radial_grid_to_as_log_grid(r_grid(i_grid_1),  r_grid_0, log_coeff), spl_param, nrad )
    end do
  end do  

  ! The atomic potentials
  ! It is necessary to include the first point (i.e. the origin) in the spline to get agreement with sratom's output
  ! This confuses me; shouldn't the potential being set to 0 (actually infinite) at the origin mess with the fit?
  ! The Hartree and electrostatic contributions
  call cubic_spline( atom_pot_hart_as, nrad, spl_param )
  do i_grid_1 = 1, n_grid
    atom_pot_hart_test(i_grid_1) = val_spline( convert_radial_grid_to_as_log_grid(r_grid(i_grid_1), &
                                      r_grid_0, log_coeff), spl_param, nrad )
    atom_pot_es_test(i_grid_1) = -1.0d0*nuclear_charge/r_grid(i_grid_1) + atom_pot_hart_test(i_grid_1)
  end do
  ! The XC potential (which is garbage for HF/hybrids!) and the effective potential
  do i_spin = 1, n_spin
    call cubic_spline( atom_pot_xc_as( 1, i_spin ), nrad, spl_param )
    do i_grid_1 = 1, n_grid
      atom_pot_xc_test( i_grid_1, i_spin ) = val_spline( convert_radial_grid_to_as_log_grid(r_grid(i_grid_1), &
                                        r_grid_0, log_coeff), spl_param, nrad )
      atom_pot_test( i_grid_1, i_spin ) = atom_pot_es_test(i_grid_1) + atom_pot_xc_test( i_grid_1, i_spin )
    end do
  end do

  ! Finally, the eigenvalues and wave functions
  do i_orbital_1 = 1, norb
    ss=so(i_orbital_1)
    nn=no(i_orbital_1)
    ll=lo(i_orbital_1)
    zz=zo(i_orbital_1)
    if (ss.gt.0.d0)isp=1
    if (ss.lt.0.d0)isp=2

    ! If aims uses a non-spin-polarized atom, make sure we only
    !  calculate quantities for one spin channel of the orbital
    if ( (n_spin.eq.1) .and. (isp.eq.2) ) then
      continue
    end if

    i_index = as_to_aims_orb_order(i_orbital_1)

    ! We average the wavefunction over the spin channels when no spin
    ! polarization is used, since formally the atomic solver should not 
    ! introduce spin-polarization given a non-spin-polarized initial 
    ! condition, and any difference should be a numerical artifact.  
    ! Practically, however, there could be a bug or the calculation may
    ! not be converged.

    ! The eigenvalues
    if (n_spin.eq.1) then
      eigenval(i_index) = 0.5d0 * ( eval(nn,1,ll+1) + eval(nn,2,ll+1) )
    else
      eigenval(i_index) = eval(nn,isp,ll+1)
    end if

    ! The wavefunctions
    if (n_spin.eq.1) then
      do i_grid_1 = 1, n_grid
        wave_as(i_grid_1,1,nn,ll+1) = 0.5d0 * ( wave_as(i_grid_1,1,nn,ll+1) + wave_as(i_grid_1,2,nn,ll+1) )
      end do
    end if

    call cubic_spline( wave_as(1,isp,nn,ll+1), nrad, spl_param )
    do i_grid_1 = 1, n_grid
      wave_test(i_grid_1, i_index) = val_spline( convert_radial_grid_to_as_log_grid(r_grid(i_grid_1), &
                                     r_grid_0, log_coeff), spl_param, nrad )
      wave_deriv_test(i_grid_1, i_index) = &
           val_spline_deriv( convert_radial_grid_to_as_log_grid(r_grid(i_grid_1), r_grid_0, log_coeff), spl_param, nrad )
      wave_2nd_deriv_test(i_grid_1, i_index) = &
           val_spline_2nd_deriv( convert_radial_grid_to_as_log_grid(r_grid(i_grid_1), r_grid_0, log_coeff), spl_param, nrad )

      wave_2nd_deriv_test(i_grid_1, i_index) = wave_2nd_deriv_test(i_grid_1, i_index) &
           - log_coeff * wave_deriv_test(i_grid_1, i_index)
      wave_2nd_deriv_test(i_grid_1, i_index) = &
           wave_2nd_deriv_test(i_grid_1, i_index) / (log_coeff * (r_grid(i_grid_1)+r_grid_0))**2.d0
      
      didr = 1.0d0/(log_coeff*(r_grid(i_grid_1)+r_grid_0))
      wave_deriv_test(i_grid_1, i_index) = didr * wave_deriv_test(i_grid_1, i_index)
    end do

    ! The kinetic density
    if (n_spin.eq.1) then
      do i_grid_1 = 1, n_grid
        kinetic_as(i_grid_1,1,nn,ll+1) = 0.5d0 * ( kinetic_as(i_grid_1,1,nn,ll+1) + kinetic_as(i_grid_1,2,nn,ll+1) )
      end do
    end if

    call cubic_spline( kinetic_as(1,isp,nn,ll+1), nrad, spl_param )
    do i_grid_1 = 1, n_grid
      kinetic_test(i_grid_1, i_index) = val_spline( convert_radial_grid_to_as_log_grid(r_grid(i_grid_1), &
                                     r_grid_0, log_coeff), spl_param, nrad )
    end do
  end do

  ! For debugging purposes, write out output variables
!  do i_grid_1 = 1, n_grid
!    write(use_unit,*)
!    write(use_unit,*) "Grid point     ", i_grid_1, " at ", r_grid(i_grid_1)
!    write(use_unit,*) "Dens           ", density(i_grid_1,1), density_test(i_grid_1)  
!    write(use_unit,*) "Dens deriv     ", density_deriv(i_grid_1,1), density_deriv_test(i_grid_1)  
!    write(use_unit,*) "Dens 2nd deriv ", density_2nd_deriv(i_grid_1,1), density_2nd_deriv_test(i_grid_1)  
!    write(use_unit,*) "Atom pot       ", atom_pot(i_grid_1,1), atom_pot_test(i_grid_1)  
!    write(use_unit,*) "Atom pot es    ", atom_pot_es(i_grid_1), atom_pot_es_test(i_grid_1)  
!    do i_orbital_1 = 1, n_valence
!      write(use_unit,*) "Wave          ", i_orbital_1, wave(i_grid_1,i_orbital_1), wave_test(i_grid_1,i_orbital_1)
!      write(use_unit,*) "Wave     Deriv", i_orbital_1, wave_deriv(i_grid_1,i_orbital_1), wave_deriv_test(i_grid_1,i_orbital_1)
!      write(use_unit,*) "Wave 2nd Deriv", i_orbital_1, "NA", wave_2nd_deriv_test(i_grid_1,i_orbital_1)
!    end do
!    write(use_unit,*)
!  end do

!  ! Update the output variables
  do i_grid_1 = 1, n_grid
    atom_pot_es(i_grid_1) = atom_pot_es_test(i_grid_1)  
  end do

  do i_spin = 1, n_spin
    do i_grid_1 = 1, n_grid
      density(i_grid_1,i_spin) = density_test(i_grid_1,i_spin)
      density_deriv(i_grid_1,i_spin) = density_deriv_test(i_grid_1,i_spin)
      density_2nd_deriv(i_grid_1,i_spin) = density_2nd_deriv_test(i_grid_1,i_spin)
      atom_pot(i_grid_1,i_spin) = atom_pot_test(i_grid_1,i_spin) 
      do i_orbital_1 = 1, n_max_ind_fns
        wave( i_grid_1, i_orbital_1 ) = wave_test( i_grid_1, i_orbital_1)
        wave_deriv( i_grid_1, i_orbital_1 ) = wave_deriv_test( i_grid_1, i_orbital_1 )
        kinetic( i_grid_1, i_orbital_1 ) = kinetic_test( i_grid_1, i_orbital_1 )
      end do
    end do
  end do

  write(info_str, '(A)')
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  "*************************************"
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  "* COMPLETED ATOM_SPHERE CALCULATION *"
  call localorb_info( info_str )
  write(info_str, '(2X,A)')  "*************************************"
  call localorb_info( info_str )
  write(info_str, '(A)')
  call localorb_info( info_str )

end subroutine atom_sphere_wrapper
!******

! This function converts a point on the radial grid to the logarithmetic
! coordinate used by atom_sphere, that is, it inverts
!   r = r_0*exp(log_coeff*(i-1))-r_0
! This inversion is needed for splining atom_sphere's results back on to aims'
! grid
real*8 function convert_radial_grid_to_as_log_grid(r_input, r_grid_0, log_coeff) result(i_as)
  real*8, intent(in)  :: r_input
  real*8, intent(in)  :: r_grid_0
  real*8, intent(in)  :: log_coeff

  i_as = 1.0d0 + 1.0d0/log_coeff * log(1.0d0 + r_input/r_grid_0)

end function
