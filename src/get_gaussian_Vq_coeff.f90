!----------------------------------------------------------------------------
!****s* FHI-aims/get_gaussian_Vq_coeff
!  NAME
!    get_gaussian_Vq_coeff
!  SYNOPSIS

subroutine get_gaussian_Vq_coeff(qvec, recip_lattice_vector, chi, n_Dvec, Dvec, &
&                          gamma, max_L1, max_L2, max_n_LM, Vq_coeff)

  !  PURPOSE
  !
  !    Given the Gaussian orbitals
  !
  !       f_lm(r) = N_l e^{-gamma r^2/2} Y_lm(r)
  !
  !    normalized to unity multipole moment
  !
  !       1 == p_f := (2l+1)^-1 \int_0^\infty dr r^2 N_l e^{-gamma r^2/2}
  !
  !    by
  !
  !       N_l = gamma^(l+1.5) / [sqrt(pi/2) (2l-1)!!],
  !
  !    calculate the q-space matrix elements
  !
  !       V_q = \sum_R V(R) e^{-iqR} = Vc(q)+v^(1)/q + v^(2)/q^2
  !
  !    with
  !
  !       V(R) = \iint d^3r d^3r' f_lm(r) f_l'm'(r-R) / |r-r'|.
  !
  !    Thus Vc(q) is regular at q=0, v^(1) and v^(2) are the expansion
  !    coefficients for the 1/q and 1/q^2 term.
  !
  !  USES

  use bravais
  use constants, only: img_unit, pi
  use triple_Y
  use mpi_tasks, only: check_allocation
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: qvec(3)
  real*8, intent(IN) :: recip_lattice_vector(3,3)
  real*8, intent(IN) :: chi
  integer, intent(IN) :: n_Dvec
  real*8, intent(IN) :: Dvec(3, n_Dvec)
  real*8, intent(IN) :: gamma
  integer, intent(IN) :: max_L1, max_L2
  integer, intent(IN) :: max_n_LM
  complex*16, intent(OUT) :: Vq_coeff(max_n_LM, max_n_LM, n_Dvec)

  !  INPUTS
  !    o qvec -- q-point
  !    o recip_lattice_vector -- reciprocal lattice
  !    o chi -- Singularity corrected value of 1/q^2 at q=0.
  !    o n_Dvec -- Number of atomic distance vectors.
  !    o Dvec -- Atomic distance vectors.
  !    o gamma -- Gaussian decay coefficient in exp(-gamma*r**2/2.).
  !    o max_L1, max_L2 -- Maximum L for which to calculate interaction.
  !    o max_n_LM -- Array dimensions of Vq;
  !                  max_n_LM >= (max(max_L1, max_L2) + 1)**2
  !  OUTPUTS
  !    o Vq -- q-dependent Coulomb interaction matrix.
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  real*8 :: triple_Yr(-max_L1:max_L1, -max_L2:max_L2, 2)
  integer :: Ms(-max_L1:max_L1, -max_L2:max_L2, 2)
  integer :: max_L, n_LM1, n_LM2, n_LM, L1, L2, LL, M1, M2, MM, i_MM
  integer :: i_LM1, i_LM2, i_LLMM
  integer :: n_add_ksq, i_add_ksq, this_max_L
  real*8 :: recip_cell_volume  , exponent
  real*8 :: kpower
  real*8 :: norm_1, norm_2
  logical :: exclude_Gamma
  complex*16, allocatable :: lattice_sum(:,:)
  complex*16 :: phase, cfac
  integer :: info
  real*8, external :: double_factorial
  character(*), parameter :: func = 'get_gaussian_Vq_coeff'

  max_L = max_L1 + max_L2
  n_LM1 = (max_L1 + 1)**2
  n_LM2 = (max_L2 + 1)**2
  n_LM = (max_L + 1)**2
  n_add_ksq = min(max_L1, max_L2)
  call get_cell_volume(recip_lattice_vector, recip_cell_volume)

  allocate(lattice_sum(n_LM, n_Dvec), stat=info)
  call check_allocation(info, 'lattice_sum', func)

  ! --- Perform lattice sum

  exponent = 1d0 / gamma   ! 1/2gamma + 1/2gamma
  exclude_Gamma = .false.
!  do i_add_ksq = 0, n_add_ksq
!     power = -2.d0 + 2*i_add_ksq
!     power =  2*i_add_ksq
!     this_max_L = max_L - 2*i_add_ksq
     this_max_L = max_L 
!     exclude_Gamma = (all(abs(qvec) < 1d-10) .and. i_add_ksq == 0)
    call gaussian_lattice_gamma(recip_lattice_vector, -qvec, &
    &                         exponent, this_max_L, n_LM, &
    &                         exclude_Gamma, n_Dvec, Dvec, &
    &                         lattice_sum(:,:))
!  end do

!  if (all(abs(qvec) < 1d-10)) then    ! Gamma point correction
!     lattice_sum(1,:,0) = lattice_sum(1,:,0) + chi / sqrt(4*pi) &
!     &                                       - exponent / sqrt(4*pi)
!     write(use_unit,*)"singularity contribution", chi / sqrt(4*pi) -  exponent / sqrt(4*pi)
!     lattice_sum(1,:,1) =  lattice_sum(1,:,1) + chi(2) / sqrt(4*pi)
     ! The exponent/sqrt(4*pi) term comes in from Taylor expanding the
     ! Gaussian at zero to first order in k^2, which cancels the k^-2 from the
     ! Coulomb interaction.  While terms of this order are in principle
     ! neglected, we make one exception here.  Without this, the final result
     ! would (slightly) depend on the Ewald gamma.  -- JW
!  end if


  ! --- Transform to normalized YY products

  Vq_coeff = (0.d0, 0.d0)
  do L1 = 0, min(1,max_L1)
!  do L1 = 0, 1
     ! Normalization is done to unity multipole moment.
     norm_1 = double_factorial(2*L1-1) * sqrt(pi / 2.d0)
     do L2 = 0, min(1,max_L2)
!     do L2 = 0, 1
        norm_2 = double_factorial(2*L2-1) * sqrt(pi / 2.d0)
        phase = img_unit**(L2-L1)
        cfac = phase * (4*pi * recip_cell_volume / (norm_1 * norm_2))
!        if(L1+L2.le.1) then
!          kpower=1.d0
!        else
!          kpower=(sqrt(dot_product(qvec,qvec)))**(L1+L2-2)
!        endif
!        cfac=cfac*kpower
        
        do LL = abs(L1-L2), L1+L2, 2
           i_add_ksq = (L1+L2-LL) / 2
           call triple_Y_real(L1, L2, LL, max_L1, max_L2, triple_Yr, Ms)
           do M1 = -L1, L1
              i_LM1 = L1**2 + L1 + M1 + 1
              do M2 = -L2, L2
                 i_LM2 = L2**2 + L2 + M2 + 1
                 do i_MM = 1, 2
                    MM = Ms(M1, M2, i_MM)
                    if (MM >= -LL) then
                       i_LLMM = LL**2 + LL + MM + 1
                       Vq_coeff(i_LM1, i_LM2, :) = Vq_coeff(i_LM1, i_LM2, :) &
                       &       + triple_Yr(M1, M2, i_MM) &
                       &         * cfac * lattice_sum(i_LLMM, :)
!                      if(i_LM1 .eq. i_LM2) then
!                      write(use_unit,'(4I4,7f12.6)') i_LM1, i_LM2, i_LLMM, i_add_ksq, triple_Yr(M1, M2, i_MM), cfac, &
!                            lattice_sum(i_LLMM, 1, i_add_ksq), Vq(i_LM1, i_LM2, 1)
!                     endif
                    end if
                 end do
              end do
           end do
        end do
     end do
  end do
  deallocate(lattice_sum)
!  norm_1 = double_factorial(1) * sqrt(pi / 2.d0)
!  write(use_unit,*)recip_cell_volume, norm_1, 4*pi * recip_cell_volume / (norm_1 * norm_1)
! write(use_unit,*) "Vq", qvec(:)
! do i_LM1 = 1, max_n_LM, 1
!   write(use_unit,'(I4, 2f16.8)') i_LM1, Vq(i_LM1, i_LM1, 1)
! enddo

end subroutine get_gaussian_Vq_coeff
!******
