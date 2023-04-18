!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Compute the dipolar or the direct spin-spin coupling tensors.
!!
!!                   1     3 R_AB R_AB^T
!!  H = alpha^2 ( ------ - ------------- ),
!!                R_AB^3       R_AB^5
!!
!!  where R_AB is the distance vector between nuclei A and B.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
subroutine compute_dipolar_couplings(dipolar_tensor)

  use constants, only: light_speed_sq
  use geometry,  only: coords
  use MR_global, only: active_nuclei
  use tools,     only: norm2
  use types,     only: dp

  implicit none

  real(dp), intent(out) :: &
       & dipolar_tensor(3,3,size(active_nuclei),size(active_nuclei))
  real(dp) :: R_AB(3)
  integer :: i_dir, j_dir, i_atom, j_atom
  do i_atom = 1, size(active_nuclei)
     do j_atom = i_atom+1, size(active_nuclei)
        R_AB = coords(:,active_nuclei(i_atom)) - coords(:,active_nuclei(j_atom))
        do i_dir = 1, 3
           dipolar_tensor(i_dir,i_dir,i_atom,j_atom) = &
                & norm2(R_AB)**(-3) - 3*R_AB(i_dir)*R_AB(i_dir)/norm2(R_AB)**5
           do j_dir = 1, 3
              if (i_dir /= j_dir) dipolar_tensor(i_dir,j_dir,i_atom,j_atom) = &
                   & -3*R_AB(i_dir)*R_AB(j_dir)/norm2(R_AB)**5
           end do
        end do
     end do
  end do
  dipolar_tensor = dipolar_tensor/light_speed_sq
end subroutine compute_dipolar_couplings
