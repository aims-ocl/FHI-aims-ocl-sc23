!****s* FHI-aims/evaluate_first_order_rho_polar_reduce_memory
!  NAME
!   evaluate_first_order_rho_polar_reduce_memory
!  SYNOPSIS

subroutine evaluate_first_order_rho_polar_reduce_memory( & 
           n_points, n_compute_c, i_basis_index, & 
           wave, &
           first_order_density_matrix_compute, &
           first_order_rho)

!  PURPOSE
!  Evaluates first_order_rho = first_order_DM*wave*wave
!  shanghui  2020.02.27

!  USES

  use dimensions
  use runtime_choices
  use basis ! basis_atom
  use pbc_lists
  implicit none

!  ARGUMENTS

  integer, intent(in) :: n_points
  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8 ,  intent(in) :: first_order_density_matrix_compute(n_compute_c,n_compute_c)

  real*8, dimension(n_points), intent(inout) :: first_order_rho


! INPUTS
! o n_points -- number of grid points
! o first_oredr_density_matrix -- relevant numbers of density matrix

!
!  OUTPUT
! o first_order_rho : first_order_rho= rho_gradien+ first_order_DM*wave*wave
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



  !     counters
  integer :: i_point
  integer :: i_state 
  integer :: i_basis,j_basis,i_compute,j_compute
  real*8, external :: ddot
  real*8 :: tmp_rho(n_compute_c, n_points)

  !  begin work
  if (n_compute_c.eq.0) then
     first_order_rho(1:n_points) = 0.0
     return
  else

  tmp_rho = 0.d0
  call dgemm("N","N",n_compute_c,n_points,n_compute_c,&
             1.d0,first_order_density_matrix_compute,n_compute_c,&
             wave(1:n_compute_c,:),n_compute_c,0.d0,tmp_rho,n_compute_c)

   ! The following amounts to considering only the diagonal elements of the result of a matrix multiplication

   !-------(v3) this is shanghui's version for degmm + tmp_rho * wave, which is even faster than ddot
     do i_point = 1, n_points,1
          first_order_rho(i_point) = &
           dot_product(tmp_rho(1:n_compute_c,i_point), wave(1:n_compute_c,i_point))
     end do

  endif
 
! old method: nested loops (slow)

!  first_order_rho(1:3,1:n_points) = 0.0d0
!
! !-------------second: first_order_rho---------------------------------
!  do i_coord = 1, 3, 1 
!
!       do i_point=1, n_points,1
!         
!          !-------shanghui add for safe---------------
!          first_order_rho(i_coord,  i_point)=0.0d0
!
!         do i_compute=1,n_compute_c
!         do j_compute=1,n_compute_c
!         
!          i_basis=i_basis_index(i_compute)
!          j_basis=i_basis_index(j_compute)          
!
!          first_order_rho(i_coord,  i_point)= &
!          first_order_rho(i_coord,  i_point)+ & 
!          first_order_density_matrix(i_coord, i_basis, j_basis) * & 
!          wave(i_compute,i_point)*wave(j_compute,i_point)           
!         enddo
!         enddo
!
!
!
!      end do  ! i_point
!  end do   ! i_coord


end subroutine evaluate_first_order_rho_polar_reduce_memory
!---------------------------------------------------------------------
!******	 
