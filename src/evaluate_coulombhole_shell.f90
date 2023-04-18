!****s* FHI-aims/evaluate_coulombhole_shell
!  NAME
!   evaluate_coulombhole_shell
!  SYNOPSIS

      subroutine evaluate_coulombhole_shell &
           ( n_points, n_KS_states, &
             partition, &
             n_prod_compute, i_basbas, &
             KS_wave, wave, coulomb_matr, &
             coulomb_hole &
           )

!  PURPOSE
!  Subroutine evaluate_coulombhole_shell evaluates the static coulomb
!  hole part of the self energy:
!     sig_coh = < i | W(r,r) |i >, |i> being the KS state,
!  for a radial shell, and adds it to the overall
!  coulomb hole contribuation.
!
!  USES
      use dimensions
      use prodbas
      use constants

      implicit none

!  ARGUMENTS

      integer :: n_points
      integer :: n_KS_states

      integer n_prod_compute
      integer i_basbas(n_prod_compute)

      real*8 partition(n_points)
      real*8 wave(n_prod_compute,n_points)
      real*8 KS_wave(n_KS_states,n_points)
      real*8 coulomb_matr( n_basbas,n_basbas )

      real*8 coulomb_hole( n_KS_states )

!  INPUTS
!  o  n_points -- number of relevant points in the current radial integration shell
!  o  n_KS_states -- number of KS states for which the Coulomb hole contribution
!        needs to calculated
!  o  n_prod_compute -- number of relevant (non-zero) auxiliary basis funcitions for
!        the current integration shell
!  o  i_basbas -- specifies the non-zero auxiliary functions
!  o  partition -- the partition function used for the numerical integral
!  o  wave -- the value of the nonzero regular basis function in the current radial shell
!  o  KS_wave -- the value of the nonzero KS orbitals in the current radial shell
!  o  coulomb_matr -- the Coulomb interaction matrix  within the auxiliary basis
!  
!  INPUTS/OUTPUTS
!  o  coulomb_hole -- the Coulomb hole contribution to the self-energy
!
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

!  local variables

!     auxiliary matrices for Level 3 Blas matrix multiplications

! ??? Need we store partition * wave in an extra matrix? ??? need BLAS descriptions to go on ...
! ??? Does it matter that n_points (the summation index in the matrix multiplication) is
! ??? the outer one, not the inner one, in wave?
!     Yes, we need, but as long as we are inside the cache, no harm should follow...

      real*8 aux_coulomb_matr &
             (n_prod_compute,n_prod_compute)

      real*8 aux_coulomb_r &
             (n_prod_compute,n_points)

      real*8 scr_coulomb_r( n_points )

!     counters

      integer :: i_compute
      integer :: i_compute_1
      integer :: i_point
      integer :: i_state

!     begin work

!     Now Integrate:
!     We only integrate the upper triangle -- all matrices must be Hermitian.

!test
!      write(use_unit,*) "psi_times_psi(1): "
!      do i_point = 1, n_points, 1
!        if (i_basbas(1).eq.1) then
!          write(use_unit,*) i_point, psi_times_psi(1,i_point),
!     +                 psi_times_psi(n_basis,i_point)
!        end if
!      enddo
!test end

!   Filtering out the relevant matrix elements
      do i_compute = 1, n_prod_compute
       do i_compute_1 = 1, n_prod_compute
        aux_coulomb_matr(i_compute, i_compute_1) = &
          coulomb_matr(i_basbas(i_compute), i_basbas(i_compute_1))
       enddo
      enddo


      aux_coulomb_r = 0.0d0

!    evaluate W_v(r) =  \sum_u P_u(r) W_uv
      call dgemm('N', 'N', n_prod_compute, n_points, &
           n_prod_compute, 1.0d0, &
           aux_coulomb_matr, n_prod_compute, wave, &
           n_prod_compute, 0.0d0, aux_coulomb_r, &
           n_prod_compute )

!    evaluate W(r,r) =  \sum_v W_v(r) P_v(r) and
!    integrate Sigcoh_ii = <i|W(r,r)|i>
      do i_point = 1, n_points
        scr_coulomb_r(i_point) = 0.d0
        do i_compute =1, n_prod_compute
         scr_coulomb_r(i_point) = &
          scr_coulomb_r(i_point) + &
          aux_coulomb_r(i_compute, i_point) * &
          wave(i_compute, i_point)
        enddo

        do i_state = 1, n_KS_states
         coulomb_hole(i_state) = &
           coulomb_hole(i_state) + &
           scr_coulomb_r(i_point) * &
           KS_wave(i_state, i_point)* &
           KS_wave(i_state, i_point)* &
           partition(i_point)
        enddo
!         write(use_unit,*) i_point
!         write(use_unit,*) scr_coulomb_r(i_point)
!         write(use_unit,*) partition(i_point)
!         write(use_unit,*) KS_wave(1:2,i_point)
!         write(use_unit,*) coulomb_hole(1:2)
      enddo
!      do i_compute=1, n_KS_states
!        coulomb_matr(1,1) = coulomb_matr(1,1)+ KS_wave(i_compute,1)*
!     +                      KS_wave(i_compute,1)
!        coulomb_matr(1,2) = coulomb_matr(1,2)+KS_wave(i_compute,1)*
!     +                      KS_wave(i_compute,2)
!      enddo
!      write(use_unit,*)n_prod_compute
!      write(use_unit,*)"wave", wave(:,1)
!      write(use_unit,*)"W", scr_coulomb_r(:)
!      write(use_unit,*)"coulomb_hole", coulomb_hole(1:2)

      return
      end subroutine evaluate_coulombhole_shell
!---------------------------------------------------------------------
!******
