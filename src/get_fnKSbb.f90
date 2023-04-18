!****s* FHI-aims/get_fnKSbb
!  NAME
!    get_fnKSbb
!  SYNOPSIS

subroutine get_fnKSbb(n_state_block, n_state_dim, basbas2col, &
&                     coeff_fnfnbb_ten, KS_vec, coeff_fnKSbb)

  !  PURPOSE
  !
  !     Contract sparse 3fn coefficients in coeff_fnfnbb_ten with a block of
  !     KS eigenvectors (assumed to be multiplied by sqrt(occ)) and return
  !     fnKSbb.  Parallel distribution is on the third index (basbas) alone.
  !
  !  USES

  use dimensions
  use prodbas
  use sparse_tensor
  use mpi_tasks, only: aims_stop, myid
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_state_block, n_state_dim
  integer, intent(IN) :: basbas2col(n_basbas)
  type(sp_ten), intent(IN) :: coeff_fnfnbb_ten
  real*8, intent(IN) :: KS_vec(n_basis, n_state_dim)
  real*8, intent(OUT) :: coeff_fnKSbb(n_basis, n_state_dim, n_loc_prodbas)


  !  INPUTS
  !    o n_state_block -- Number of states in this block
  !    o n_state_dim -- Array dimension of KS_vec and coeff_fnKSbb.
  !    o basbas2col -- i_basbas -> i_loc_prodbas
  !    o coeff_fnfnbb_ten -- Sparse tensor of 3fn coefficients.
  !                          saved as 'U' with diagonal (fnfn-symmetry).
  !    o KS_vec -- Part of KS_eigenvectors, possibly multiplied by sqrt(occ).
  !  OUTPUTS
  !    o coeff_fnKSbb -- 2fn1KS array.
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
  !    Release version, FHI-aims (2010).
  !  SOURCE

  integer :: i_nzb, top(3), shp(3), off, str_2, str_prod
  integer :: i_basbas, i_loc_prodbas, off_prod, off_prod2
  integer :: i_basis_1, i_basis_2, i_state, i_val
  character(*), parameter :: func = 'get_fnKSbb'

  if (any(coeff_fnfnbb_ten%top(1,:) > coeff_fnfnbb_ten%top(2,:))) then
     call aims_stop('Tensor not upper triangular', func)
  else if (allocated(coeff_fnfnbb_ten%str)) then
     call aims_stop('Non-default stride not implemented', func)
  end if

  coeff_fnKSbb = 0.d0
  do i_nzb = 1, coeff_fnfnbb_ten%n_nzb
     top = coeff_fnfnbb_ten%top(:, i_nzb)
     shp = coeff_fnfnbb_ten%shp(:, i_nzb)
     off = coeff_fnfnbb_ten%off(i_nzb)
     str_2 = shp(1)
     str_prod = str_2 * shp(2)

     ! Could this looping structure go into its own sub ...
     do i_basbas = top(3), top(3)+shp(3)-1
        i_loc_prodbas = basbas2col(i_basbas)
        if (i_loc_prodbas <= 0) call aims_stop('Non-local index', func)
        if (map_prodbas(i_loc_prodbas, myid+1) /= i_basbas) then
           call aims_stop('Wrong local index', func)
        end if
        off_prod = off + (i_basbas - top(3)) * str_prod
        do i_basis_2 = top(2), top(2)+shp(2)-1
           off_prod2 = off_prod + (i_basis_2 - top(2)) * str_2
           do i_state = 1, n_state_block

              ! coeff_fnfnbb is saved as 'U' with diagonal.

              ! Off-diagonal: Make sure that i_basis_2 > i_basis_1
              i_val = off_prod2
              do i_basis_1 = top(1), min(i_basis_2-1, top(1)+shp(1)-1)
                 i_val = i_val + 1

                 ! Direct: fnKS[bb] += fnfn^T[bb] * fnKS
                 coeff_fnKSbb(i_basis_2, i_state, i_loc_prodbas) = &
                 & coeff_fnKSbb(i_basis_2, i_state, i_loc_prodbas) + &
                 &   coeff_fnfnbb_ten%val(i_val) * &
                 &   KS_vec(i_basis_1, i_state)

                 ! Transposed: fnKS[bb] += fnfn[bb] * fnKS
                 coeff_fnKSbb(i_basis_1, i_state, i_loc_prodbas) = &
                 & coeff_fnKSbb(i_basis_1, i_state, i_loc_prodbas) + &
                 &   coeff_fnfnbb_ten%val(i_val) * &
                 &   KS_vec(i_basis_2, i_state)
              end do

              ! Diagonal
              if (i_basis_2>=top(1) .and. i_basis_2<top(1)+shp(1)) then
                 i_basis_1 = i_basis_2
                 i_val = off_prod2 + (i_basis_1 - top(1)) + 1
                 coeff_fnKSbb(i_basis_2, i_state, i_loc_prodbas) = &
                 & coeff_fnKSbb(i_basis_2, i_state, i_loc_prodbas) + &
                 &   KS_vec(i_basis_1, i_state) * &
                 &   coeff_fnfnbb_ten%val(i_val)
              end if
           end do
        end do
     end do
  end do

end subroutine get_fnKSbb
!******
