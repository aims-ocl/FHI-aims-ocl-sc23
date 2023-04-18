!****s* FHI-aims/evaluate_first_order_rho_p1
!  NAME
!   evaluate_first_order_rho_p1
!  SYNOPSIS

subroutine evaluate_first_order_rho_p1( & 
           n_points,n_compute_c, i_basis_index, & 
           wave,gradient_basis_wave, &
           density_matrix_compute,first_order_density_matrix_compute, & 
           first_order_rho)

!  PURPOSE
!  Evaluates first_order_rho= rho_gradien+ first_order_DM*wave*wave
!  shanghui  2012.05.02

!  USES

  use dimensions
  use runtime_choices
  use basis ! basis_atom
  use pbc_lists
  use mpi_tasks, only: check_allocation
  implicit none

!  ARGUMENTS

  integer :: n_points

  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8, dimension(n_compute_c,n_compute_c),intent(in) :: density_matrix_compute
  real*8, dimension(3,n_centers_in_sc_DFPT,n_compute_c,n_compute_c),intent(in) :: & 
                   first_order_density_matrix_compute

  real*8, dimension(3,n_centers_in_sc_DFPT,n_points), intent(out) :: first_order_rho



! INPUTS
! o n_points -- number of grid points
! o first_oredr_density_matrix_sparse -- relevant numbers of density matrix
! o wave -- basis functions
! o gradient_basis_wave -- gradient of basis functions
! o KS_ev_compute -- KS eigenvector

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
  integer :: i_coord
  integer :: i_point
  integer :: i_state,info 
  integer :: i_basis,j_basis,i_compute,j_compute,i_basis_uc,j_basis_uc,i_cell,j_cell
  integer :: i_center, j_center
  integer :: i_place


  real*8, dimension(:,:,:), allocatable :: rho_gradient_atom 
  real*8, allocatable ::  work(:,:)
  

 !     external functions
  real*8, external ::  ddot

   !----------------allocate local -------------------------------
   allocate(rho_gradient_atom(3, n_centers_in_sc_DFPT, n_points),stat=info) 
   call check_allocation(info, 'rho_gradient_atom')        
   allocate(work(n_compute_c, n_points))


  rho_gradient_atom(1:3,1:n_centers_in_sc_DFPT,1:n_points)  = 0.0d0
  first_order_rho(1:3,1:n_centers_in_sc_DFPT,1:n_points) = 0.0d0


  !--------------------first: rho_gradient_atom----------------------

   do i_compute = 1,n_compute_c
      i_basis   = i_basis_index(i_compute)
      i_center  = center_to_center_in_sc_DFPT(Cbasis_to_center(i_basis))
   do j_compute = 1,n_compute_c
      j_basis   = i_basis_index(j_compute)
      j_center  = center_to_center_in_sc_DFPT(Cbasis_to_center(j_basis))
 
         do i_coord = 1, 3, 1
         do i_point=1, n_points,1

          rho_gradient_atom(i_coord, i_center, i_point)= &
          rho_gradient_atom(i_coord, i_center, i_point)+ & 
          density_matrix_compute(i_compute, j_compute) * & 
          (-gradient_basis_wave(i_compute,i_coord,i_point)*wave(j_compute,i_point)) 

          rho_gradient_atom(i_coord, j_center, i_point)= &
          rho_gradient_atom(i_coord, j_center, i_point)+ & 
          density_matrix_compute(i_compute, j_compute) * & 
          (-gradient_basis_wave(j_compute,i_coord,i_point)*wave(i_compute,i_point)) 

         enddo
         enddo

   enddo
   enddo

  


  !-------------second: first_order_rho---------------------------------
   do i_coord = 1, 3, 1 
   do i_center = 1, n_centers_in_sc_DFPT, 1

       ! \sum_{u} wave(u,r) [\sum_{v} DM^{(1)}_{uv} wave(v,r) ]----------------
      work = 0.0d0
      call dgemm("N","N", n_compute_c, n_points, n_compute_c,&  !transa, transb, l, n, m,
           1.0d0,first_order_density_matrix_compute(i_coord, i_center, 1:n_compute_c, 1:n_compute_c), & 
           n_compute_c, wave(1:n_compute_c,1:n_points), n_compute_c,& !alpha, a, lda, b, ldb 
           0.0d0, work,n_compute_c)  ! beta, c, ldc

      do i_point = 1, n_points,1
          first_order_rho(i_coord, i_center, i_point) = rho_gradient_atom(i_coord, i_center, i_point) +  &
          dot_product(wave(1:n_compute_c,i_point), work(1:n_compute_c,i_point))
      end do



!---------------begin old version--------------
!       do i_point=1, n_points,1
!        do i_compute=1,n_compute_c
!        do j_compute=1,n_compute_c
!        
!         i_basis=i_basis_index(i_compute)
!         j_basis=i_basis_index(j_compute)          
!
!        first_order_DM_rho(i_coord, i_center, i_point)= &
!        first_order_DM_rho(i_coord, i_center, i_point)+ & 
!        first_order_density_matrix_compute(i_coord, i_center, i_compute, j_compute) * & 
!        wave(i_compute,i_point)*wave(j_compute,i_point)           
!
!        enddo
!        enddo
!      first_order_rho(i_coord, i_center, i_point) =                     &
!      rho_gradient_atom(i_coord, i_center, i_point)  +                    &
!      first_order_DM_rho(i_coord, i_center, i_point)
!      end do  ! i_point
!---------------end old version--------------
 
   end do ! i_center
   end do ! i_coord


  if(allocated(rho_gradient_atom)) deallocate (rho_gradient_atom ) 
  if(allocated(work))  deallocate(work)

end subroutine evaluate_first_order_rho_p1
!---------------------------------------------------------------------
!******	 
