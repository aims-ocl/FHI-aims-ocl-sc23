!****s* FHI-aims/integrate_errorfunction
!  NAME
!   integrate_errorfunction
!  SYNOPSIS

      subroutine integrate_errorfunction(n_legendre, n_grid, n_max_grid, &
      &                                  r_grid, omega, expansion_matrix)

!  PURPOSE
!     subroutine to calculate spherical harmonics expansion for
!     erfc(\omega * (r1 - r2))/(r1 - r2)
!
!     based on j.A.Angyan, I.Gerber, M.Marsman, J.Phys.A:Math.Gen 39,8613-8630 (2006)
!
!  USES

      use mpi_tasks
      use constants
      use arch_specific
      use localorb_io,only:use_unit
      implicit none

!  ARGUMENTS

      integer, intent(IN) :: n_legendre
      integer, intent(IN) :: n_grid, n_max_grid
      real*8, intent(IN) :: r_grid(n_grid)
      real*8, intent(IN) :: omega
      real*8, intent(OUT) :: expansion_matrix(n_max_grid, n_max_grid, &
      &                                       n_legendre+1)
!  INPUTS
!  o n_legendre -- the highest order of the Legendre expansion (max_l)
!  o n_grid -- number of grid points
!  o n_max_grid -- array dimension
!  o r_grid -- actual grid points
!  o omega -- HSE paramter
!  OUTPUTS
!  o  expansion_matrix -- the expansion of error function in terms of Legendre
!        polynomials ?
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

      real*8 :: Xi_great, Xi_small
      real*8 :: Xi_prod
!
      real*8 :: upper_argument, lower_argument
      real*8 :: upper_exp, lower_exp
      real*8 :: upper_erfc,lower_erfc
!
!     EXTERNAL FUNCTIONS
!
      real*8 :: double_factorial
      real*8 :: factorial
      real*8 :: binomial_coefficient_neg
!
!     CUTOFF VARIABLES
!
      real*8  :: l_cut_cond
      integer :: l_cut
!
!     VARIABLES FOR CALCULATE AUX_F
!
      real*8 :: F_prefactor_array(n_legendre+1,n_legendre+1)
      real*8 :: F_prefactor_const, F_prefactor
      real*8 sum_up_low(2)
      real*8 :: aux_F_fn
      integer :: F_index
!
!     VARIABLES FOR CALCULATE AUX_H
!
      real*8 :: H_prefactor, H_prefactor_const
      real*8 :: H_upper_part, H_lower_part
      real*8  aux_H_fn(n_legendre+1)
      integer :: H_exponent
!
!     VARIABLE FOR GENERELL INTEGRATION
!
      real*8  greens_fn_array(n_legendre+1)
      real*8  greens_fn_array_aux(n_legendre+1)
      real*8 :: prefactor
      real*8 :: greens_fn,aux_greens_fn
!
!     VARIABLES FOR DAMPINGFUNCION
!
      real*8 :: damping_fn
      real*8 :: aux_damping_fn
      real*8 :: damp_prefactor
      integer :: i_damping
!
!     COUNTER
!
      integer :: i_order, j_order
      integer :: i_radial,j_radial
!
!     begin work
!
      do i_order = 0, n_legendre,1
!
         do j_order = 0,n_legendre,1
!
            F_prefactor_array(i_order+1,j_order+1) = &
                 factorial(i_order+j_order)/ &
                 (factorial(j_order)*factorial(i_order-j_order))
!
         enddo
!
      enddo
!
!
!     loop over the min-radial
!
      do i_radial = 1, n_grid
!
         Xi_small = omega * r_grid(i_radial)
!
!     loop over the max-radial
!
!            write(use_unit,*)'|  #  |i_radial',i_radial
!
!
         do j_radial = i_radial, n_grid
!
            Xi_great = omega * r_grid(j_radial)
!
!
!     FIXME: change order, use arch
!
!     Variables that needed in any case
!
            upper_argument = Xi_great + Xi_small
            lower_argument = Xi_great - Xi_small
!
            upper_erfc = arch_erfc(upper_argument)
            lower_erfc = arch_erfc(lower_argument)
!
            Xi_prod = Xi_small * Xi_great
!
            upper_exp = dexp(-(upper_argument*upper_argument))
            lower_exp = dexp(-(lower_argument*lower_argument))
!
            F_prefactor_const = dble(-1.d0/(4.d0*Xi_prod))
!
            sum_up_low(1) = upper_exp - lower_exp
            sum_up_low(2) = -(upper_exp + lower_exp)
!
!     First: get the cutoff_index
!
            if (Xi_prod .eq.1.0d0) then
!
               l_cut = n_legendre
!
            elseif (Xi_prod.gt.1.0d0)then
!
               l_cut = n_legendre
            else
               l_cut_cond = -3.0d0/(dlog10(Xi_prod))-0.5d0
!
               if((l_cut_cond).lt.0.0d0)then
!
                  l_cut = -1
!
               else
!
                  l_cut = min(int(l_cut_cond),n_legendre)
!
               endif
!
            endif
!
!     next: tabulate the Funktion H
!
!
            H_prefactor_const = 1.0d0/(Xi_prod)
!
            H_prefactor = 1.0d0/2.0d0
!
            do i_order = 0, l_cut,1
!
               H_prefactor = H_prefactor * H_prefactor_const
!
               H_exponent = 2*i_order+1
!
               H_upper_part =(Xi_great**(H_exponent) + &
                    Xi_small**(H_exponent)) * upper_erfc
!
               H_lower_part =(Xi_great**(H_exponent) - &
                    Xi_small**(H_exponent)) * lower_erfc
!
               aux_H_fn(i_order+1) = H_prefactor*(H_upper_part - &
                    H_lower_part)
!
            enddo
!
!     Function H done, now get Function F
!
            do i_order =0, l_cut,1
!
               F_prefactor = 1.0d0
!
               aux_F_fn = 0.d0
!
               do j_order = 0, i_order, 1
!
                  F_index = mod(i_order - j_order,2) + 1
!
                  F_prefactor = F_prefactor * F_prefactor_const
!
                  aux_F_fn = aux_F_fn &
                       +  F_prefactor * F_prefactor_array(i_order+1, &
                       j_order+1) * sum_up_low(F_index)
!
               enddo
!
               greens_fn_array(i_order+1) = aux_F_fn
!
               greens_fn_array_aux(i_order+1) = &
                  greens_fn_array(i_order+1)
            enddo
!
!     OPTIMIZED!!
!
            do i_order =1, l_cut,1
!
               F_prefactor = ( Xi_great**(2*i_order) + &
                    Xi_small ** (2*i_order) )/(Xi_prod) ** i_order
!
               do j_order = i_order,l_cut,1
!
                  greens_fn_array(j_order+1) = &
                       greens_fn_array(j_order+1) + &
                       greens_fn_array_aux(j_order-i_order+1) * &
                       F_prefactor
!
!
               enddo
            enddo
!
!
            do i_order = 0,l_cut,1
!
               prefactor = dble(2.0d0/(2.0d0 * i_order +1.0d0))
!
               prefactor = omega * prefactor
!
               greens_fn_array(i_order+1) = &
                    greens_fn_array(i_order+1) * 2.0d0 &
                    *pisqrt_inv &
                    + aux_H_fn(i_order+1)
!
               greens_fn_array(i_order+1) = &
                    greens_fn_array(i_order+1)* &
                    prefactor
!
					expansion_matrix(i_radial,j_radial,i_order+1) &
		         & = greens_fn_array(i_order+1)
!
		         expansion_matrix(j_radial,i_radial,i_order+1) &
		         & = greens_fn_array(i_order+1)

               
!
!
            enddo
!
!     END OF GENERERALL INTEGRATION, NOW BEGIN WITH THE DAMPING FUNCTION
!
            do i_order = (l_cut+1), n_legendre,1
!
               prefactor = dble(2.0d0/(2.0d0 * i_order + 1.0d0))
!
               prefactor = omega * prefactor
!
               greens_fn = 0.d0
!
               aux_greens_fn = 0.d0
!
!     FIXME: FIX THE DAMPING ORDER
!
!     first: Zeroth order of the damping function
!
               damping_fn = 0.0d0
!
               aux_damping_fn = 0.0d0
!
               do i_damping =1, i_order,1
!
                  aux_damping_fn = 2.0d0**i_damping * &
                       Xi_great**(2*i_damping) &
                       * double_factorial(2*i_order - 2*i_damping +1)
!
                  damping_fn = damping_fn+dble(1.0d0/ aux_damping_fn)
!
               enddo
!
               damping_fn = damping_fn * Xi_great**(2*i_order+1) * &
                    2.0d0**(i_order+1)*dexp(-Xi_great**2)*pisqrt_inv
!
               damping_fn = damping_fn + arch_erfc(Xi_great)
!
               aux_greens_fn = damping_fn * Xi_small** &
                    (i_order) / Xi_great**(i_order +1)
!
               greens_fn = greens_fn + aux_greens_fn
!
!     zeroth order done, now next
!
               do j_order=1,15,1
!
                  damping_fn =0.0d0
!
                  do i_damping = 1, j_order,1
!
                     aux_damping_fn = &
                          binomial_coefficient_neg(i_damping - &
                          j_order-1,i_damping-1)* &
                          2.0d0 ** (j_order -i_damping) &
                          * Xi_great**(2*(j_order - i_damping)) / &
                          double_factorial(2*i_order+2*j_order - &
                          2*i_damping +1)
!
                     damping_fn = damping_fn + aux_damping_fn
!
                  enddo
!
                  damp_prefactor = 0.d0
!
                  damp_prefactor = 2.0d0**(i_order + 1)* &
                       (2.0d0 * i_order + 1) &
                       /(2.d0*i_order+2.d0*j_order +1)
!
                  damp_prefactor = &
                       dble(damp_prefactor/factorial(j_order)) &
                       *pisqrt_inv
!
                  damping_fn = damping_fn *damp_prefactor * &
                       dexp(-Xi_great**2) * Xi_great**(2*i_order +1)
!
!     damping of j - order done now compare
!
                  aux_greens_fn = damping_fn * Xi_small** &
                       (i_order + 2 * j_order)/ Xi_great**(i_order+1)
!
                  greens_fn = greens_fn + aux_greens_fn
!
                  if (aux_greens_fn/greens_fn.lt.1.d-10) then
!
                     exit
!
                  endif
!
               enddo
!
               if (aux_greens_fn/greens_fn.gt.1.d-10) then
!
                  write(use_unit,*) "damping_fn failt"
!
                  stop
!
               endif
!
!         greens_fn_array(i_order+1) =  greens_fn * prefactor
					expansion_matrix(i_radial,j_radial,i_order+1) &
   	         & = greens_fn * prefactor
!
   	         expansion_matrix(j_radial,i_radial,i_order+1) &
   	         & = greens_fn * prefactor
!
!
            enddo
!
         enddo
!
      enddo
!
      return
!
      end subroutine integrate_errorfunction


