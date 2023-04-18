!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Generate Slater-type orbitals (STOs) from the input parameters
!!  sto_n, sto_l, and sto_zeta given by the user. Note that the
!!  central quantity here is u(r) and not u(r)/r.
!!
subroutine get_sto_basis_fns(i_species)

  use dimensions,      only: n_max_grid, use_basis_gradients
  use grids,           only: n_grid, r_grid
  use runtime_choices, only: flag_rel, out_basis, REL_x2c, REL_4c_dks
  use species_data,    only: n_sto, sto_n, sto_l, sto_k, sto_zeta, sto_wave, &
       & sto_wave_deriv, sto_kinetic, sto_wave_small, sto_kinetic_small
  use types,           only: dp

  implicit none

  interface
     real*8 function factorial(n)
       integer, intent(in) :: n
     end function factorial
  end interface

  integer, intent(in) :: i_species
  real(dp) :: wave(n_max_grid)
  real(dp) :: wave_1st(n_max_grid) ! First derivative
  real(dp) :: wave_2nd(n_max_grid) ! Second derivative
  real(dp) :: norm
  integer, parameter :: io = 27458
  character(16) :: file_base
  integer :: i_grid, i_sto

  do i_sto = 1, n_sto(i_species)
     associate(U => wave(:n_grid(i_species)), &
          & U1 => wave_1st(:n_grid(i_species)), &
          & U2 => wave_2nd(:n_grid(i_species)), &
          & R => r_grid(:n_grid(i_species),i_species), &
          & N => sto_n(i_species,i_sto), &
          & L => sto_l(i_species,i_sto), &
          & zeta => sto_zeta(i_species,i_sto))
       U = R**(N)*exp(-zeta*R)
       if (use_basis_gradients) U1 = (N - zeta*R)*R**(N-1)*exp(-zeta*R)
       U2 = (N*(N-1) - 2*zeta*N*R + zeta**2*R**2)*R**(N-2)*exp(-zeta*R)
       norm = sqrt(2*zeta/factorial(2*N))*(2*zeta)**(N)
       U = norm*U
       if (use_basis_gradients) U1 = norm*U1
       U2 = norm*U2
       sto_wave(:n_grid(i_species), i_species, i_sto) = U
       sto_kinetic(:n_grid(i_species), i_species, i_sto) = &
            & (-U2 + L*(L+1)/R**2*U)/2
       if (use_basis_gradients) &
            & sto_wave_deriv(:n_grid(i_species), i_species, i_sto) = U1
       ! Relativistic small component
       if(flag_rel == REL_x2c .or. flag_rel == REL_4c_dks)then
          call sigma_dot_p_large_wave( &
               n_grid(i_species), r_grid(1,i_species), &
               sto_l(i_species,i_sto), sto_k(i_species,i_sto), &
               sto_wave(1,i_species,i_sto), U1(1), &
               sto_wave_small(1,i_species,i_sto), &
               sto_kinetic_small(1,i_species,i_sto) )
       endif
       ! Print the orbital and the kinetic wave into files
       if (out_basis) then
          write(file_base, '(3(a, i0), i0, a)') &
               & 'S', i_species, '_', i_sto, '_', N, L, '_base'
          open(io, file=trim(file_base)//'.dat', form='formatted', &
               & action='write')
          do i_grid = 1, n_grid(i_species)
             write(io, *) r_grid(i_grid, i_species), U(i_grid)
          end do
          close(io)
          open(io, file=trim(file_base)//'_kin.dat', form='formatted', &
               & action='write')
          do i_grid = 1, n_grid(i_species)
             write(io, *) r_grid(i_grid, i_species), U2(i_grid)
          end do
          close(io)
       end if
     end associate
  end do

  if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) call check_outradius_sto(i_species)

end subroutine get_sto_basis_fns

 subroutine sigma_dot_p_large_wave(n_grid, r_grid, l, kappa, &
  wave, wave_deriv, wave_small, kinetic_small)
  use constants, only: light_speed, light_speed_sq
  implicit none
  integer,intent(in) :: n_grid, l, kappa
  real*8,intent(in)  :: r_grid(n_grid)
  real*8,intent(in)  :: wave(n_grid), wave_deriv(n_grid)
 ! wave_small: \frac{\sigma \dot p \chi^L}{2c} r
 ! kinetic_small: \frac{\sigma \dot p \chi^L}{c} r
  real*8,intent(out) :: wave_small(n_grid), kinetic_small(n_grid)

  integer :: i
  real*8 :: j, fac, c2

  if(kappa == l)then            ! spin down
     j = l - 0.5d0
  else if(kappa == (-l-1))then   ! spin up
     j = l + 0.5d0
  end if

  fac = l*(l+1.d0) - j*(j+1.d0) -0.25d0
  c2 = 0.5d0/light_speed

  do i=1, n_grid
     wave_small(i) = wave_deriv(i)*c2 + wave(i)*c2/r_grid(i)*fac
  enddo

  kinetic_small = 2.d0*wave_small*light_speed_sq
 end subroutine

 subroutine check_outradius_sto(i_species)
  use species_data, only: n_sto, sto_n, sto_l, sto_k, sto_zeta, sto_wave, &
       & sto_wave_deriv, sto_kinetic, sto_wave_small
  use grids,        only: n_grid, r_grid
  implicit none
  integer,intent(in) :: i_species
  integer :: i,j,k, i_sto
  real*8 :: thresh=1.d-8

  do i_sto=1, n_sto(i_species)
    do i=n_grid(i_species), 1, -1 
       if(dabs(sto_wave(i,i_species,i_sto)).gt.thresh)then
          sto_wave(i:n_grid(i_species),i_species,i_sto) = 0.d0
          write(6,"('sto_large=',i4,3x,'n=',i2,3x,'l=',i2,3x,'zeta=',f8.3,3x,'radius=',f13.7)")&
           i_sto,sto_n(i_species,i_sto),sto_l(i_species,i_sto), &
           sto_zeta(i_species,i_sto),r_grid(i,i_species)
          exit
       endif
    enddo
  enddo

  do i_sto=1, n_sto(i_species)
    do i=n_grid(i_species), 1, -1 
       if(dabs(sto_wave_small(i,i_species,i_sto)).gt.thresh)then
          sto_wave_small(i:n_grid(i_species),i_species,i_sto) = 0.d0
          write(6,"('sto_small=',i4,3x,'n=',i2,3x,'l=',i2,3x,'zeta=',f8.3,3x,'radius=',f13.7)")&
           i_sto,sto_n(i_species,i_sto),sto_l(i_species,i_sto), &
           sto_zeta(i_species,i_sto),r_grid(i,i_species)
          exit
       endif
    enddo
  enddo

  do i_sto=1, n_sto(i_species)
    do i=n_grid(i_species), 1, -1 
       if(dabs(sto_kinetic(i,i_species,i_sto)).gt.thresh)then
          sto_kinetic(i:n_grid(i_species),i_species,i_sto) = 0.d0
          write(6,"('sto_kinetic=',i4,3x,'n=',i2,3x,'l=',i2,3x,'zeta=',f8.3,3x,'radius=',f13.7)")&
           i_sto,sto_n(i_species,i_sto),sto_l(i_species,i_sto), &
           sto_zeta(i_species,i_sto),r_grid(i,i_species)
          exit
       endif
    enddo
  enddo

 end subroutine

