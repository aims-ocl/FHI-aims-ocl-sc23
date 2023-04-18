!****s* FHI-aims/hf_postproc
!  NAME
!    hf_postproc
!  SYNOPSIS

subroutine hf_postproc()

  !  PURPOSE
  !
  !     Calculate Hartree-Fock exchange perturbatively on top of a DFT
  !     calculation.
  !
  !     In principle, this subroutine could also be used for perturbative
  !     hybrid functionals.
  !
  !  USES

  use localorb_io
  use runtime_choices
  use dimensions
  use physics
  use mpi_tasks, only: check_allocation
  implicit none

  !  ARGUMENTS
  !    none
  !  INPUTS
  !    none
  !  OUTPUTS
  !    none
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
  !    Release version, FHI-aims (2011).
  !  SOURCE

  real*8, allocatable :: my_hamiltonian(:,:)
  real*8 :: my_fock_energy, my_en_xc, my_en_x, my_total_energy
  real*8 :: my_nlx_dc
  real*8 :: save_hybrid_coeff
  integer :: info
  character*150 :: info_str
  character(*), parameter :: func = 'hf_postproc'


  allocate(my_hamiltonian(n_basis*(n_basis+1)/2, n_spin), stat=info)
  call check_allocation(info, 'hamiltonian', func)
  my_hamiltonian = 0.d0
  my_fock_energy = 0.d0
  my_en_xc = 0.d0

  save_hybrid_coeff = hybrid_coeff
  hybrid_coeff = 1.d0
  call get_hf_hamiltonian(-1, &
  &                KS_eigenvalue, KS_eigenvector, n_electrons, occ_numbers, &
  &                my_hamiltonian, my_fock_energy, my_en_xc)
  hybrid_coeff = save_hybrid_coeff

  ! JW: For LDA/GGA, the double counting correction is in -en_pot_xc and the
  ! actual xc energy in en_xc.  For nonlocal exchange, the double counting
  ! correction amounts to Tr(KD) with the exchange matrix K and the density
  ! matrix D, so that Tr(KD) is positive.  The actual exchange energy is then
  ! -0.5*Tr(KD).  Unfortunately, en_xc contains the sum of these two terms and
  ! en_pot_xc gets nothing of it.
  my_en_x = -my_en_xc
  if (use_hartree_fock) then
     my_nlx_dc = fock_energy   ! Contains hybrid_coeff
  else
     my_nlx_dc = 0.d0
  end if

  my_total_energy = total_energy + my_nlx_dc - en_xc + my_en_x

  call localorb_info('')
  write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
  & "| Exchange energy               :", &
  & -my_en_x, " Ha", (-my_en_x)* hartree, " eV"
  call localorb_info(info_str)
  write (info_str,'(2X,A,1X,F20.8,A,1X,F20.8,A)') &
  & "| Total EX@DFT energy           :", &
  & my_total_energy, " Ha", my_total_energy*hartree, " eV"
  call localorb_info(info_str)

end subroutine hf_postproc
!******
