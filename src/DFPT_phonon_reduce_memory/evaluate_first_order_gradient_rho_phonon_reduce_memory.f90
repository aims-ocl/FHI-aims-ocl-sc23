!****s* FHI-aims/evaluate_first_order_gradient_rho_phonon_reduce_memory
!  NAME
!   evaluate_first_order_gradient_rho_phonon_reduce_memory
!  SYNOPSIS

subroutine evaluate_first_order_gradient_rho_phonon_reduce_memory( & 
           n_points,n_compute_c, i_basis_index, index_hessian, & 
           wave,gradient_basis_wave, hessian_basis_wave, &
           density_matrix_compute,first_order_density_matrix_compute, & 
           i_q_point, j_atom, j_coord, &
           first_order_gradient_rho)

!  PURPOSE
!  Evaluates first_order_gradient_rho = first_order_DM*gradient(wave*wave) + 
!                                       DM * first_order(gradient(wave*wave))
!  shanghui  2017 @ Berlin

!  USES

  use dimensions
  use runtime_choices
  use basis ! basis_atom
  use pbc_lists
  implicit none

!  ARGUMENTS

  integer :: n_points

  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  integer , intent(in) :: index_hessian(3,3)
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham, 6, n_points),intent(in) :: hessian_basis_wave

  real*8     :: density_matrix_compute(n_compute_c,n_compute_c)
  complex*16 :: first_order_density_matrix_compute(n_compute_c,n_compute_c)
  integer , intent(in) :: i_q_point
  integer , intent(in) :: j_atom
  integer , intent(in) :: j_coord
 
  complex*16, dimension(3, n_points), intent(inout) :: first_order_gradient_rho


! INPUTS
! o n_points -- number of grid points
! o first_oredr_density_matrix -- relevant numbers of density matrix
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
  integer :: i_state 
  integer :: i_basis,j_basis,i_compute,j_compute
  real*8  :: Gr(3)
  complex*16  ::  exp_iqr(3)
  real*8, external :: ddot

!---------------------------(2) This is the fast verion---------------------------
!----here we only first consider i_q_point=0, further we will add more q point---- 
 
  real*8, allocatable ::  work(:,:) 
  real*8, allocatable ::  first_order_density_matrix_compute_real(:,:) 
  real*8, allocatable ::  first_order_wave(:,:) 
  real*8, allocatable ::  first_order_gradient_wave(:,:,:) 
 

  allocate(work(n_compute_c, n_points))
  allocate(first_order_density_matrix_compute_real(n_compute_c,n_compute_c))
  allocate(first_order_wave(n_compute_c, n_points))
  allocate(first_order_gradient_wave(n_compute_c,3, n_points))

 
 
  if (n_compute_c.eq.0) then
     first_order_gradient_rho(1:3,1:n_points) = (0.0d0,0.0d0)
     return
  else
 
  first_order_gradient_rho(1:3, 1:n_points)    = (0.0d0, 0.0d0)
  do i_compute = 1,n_compute_c 
    do j_compute = 1,n_compute_c 
       first_order_density_matrix_compute_real(i_compute,j_compute)  = & 
            dble(first_order_density_matrix_compute(i_compute,j_compute)) 
    enddo 
  enddo
 
 !-------------------------(2.0) wave^{(1)}(u,r) and gradient_wave^{(1)}(u,r)--------------------
   first_order_wave(1:n_compute_c,1:n_points) = 0.0d0 
   ! first_order_wave = delta(i_compute,j_atom) d wave(i_compute)/dR(j_atom)
   do i_compute=1,n_compute_c
      i_basis=i_basis_index(i_compute)
      if(Cbasis_to_atom(i_basis).eq.j_atom) then 
        first_order_wave(i_compute,1:n_points) = -gradient_basis_wave(i_compute,j_coord,1:n_points)
      endif
   enddo

   first_order_gradient_wave(1:n_compute_c,1:3,1:n_points) = 0.0d0
   ! first_order_gradient_wave = delta(i_compute,j_atom) d gradient_wave(i_compute)/dR(j_atom)
   do i_coord = 1, 3
   do i_compute=1,n_compute_c
      i_basis=i_basis_index(i_compute)
      if(Cbasis_to_atom(i_basis).eq.j_atom) then 
        first_order_gradient_wave(i_compute,i_coord,1:n_points) = &
        -hessian_basis_wave(i_compute,index_hessian(i_coord,j_coord),1:n_points)
      endif
   enddo
   enddo 
 

 
 !------------- (2.1) begin for first_order_density_matrix_compute terms------------------------
 !                    \sum_{u} gradient_wave(u,r) [\sum_{v} DM^{(1)}_{uv} wave(v,r) ]
 !                   +\sum_{u} wave(u,r) [\sum_{v} DM^{(1)}_{uv} gradient_wave(v,r) ]         
 !--------------------end  for first_order_density_matrix_compute terms----------------------------
  do i_coord = 1,3
    work = 0.0d0
    call dgemm("N","N", n_compute_c, n_points, n_compute_c,&  !transa, transb, l, n, m,
         1.0d0,first_order_density_matrix_compute_real, n_compute_c, wave(1:n_compute_c,1:n_points), n_compute_c,& !alpha, a, lda, b, ldb 
         0.0d0, work,n_compute_c)  ! beta, c, ldc
    do i_point = 1, n_points,1
        first_order_gradient_rho(i_coord, i_point) = &
          dcmplx( dot_product(gradient_basis_wave(1:n_compute_c,i_coord,i_point), work(1:n_compute_c,i_point)), 0.0d0 )  
    end do
 
    work = 0.0d0
    call dgemm("N","N", n_compute_c, n_points, n_compute_c,&  !transa, transb, l, n, m,
         1.0d0,first_order_density_matrix_compute_real, n_compute_c, & 
         gradient_basis_wave(1:n_compute_c,i_coord,1:n_points), n_compute_c,& !alpha, a, lda, b, ldb 
         0.0d0, work,n_compute_c)  ! beta, c, ldc
    do i_point = 1, n_points,1
        first_order_gradient_rho(i_coord, i_point) = &
        first_order_gradient_rho(i_coord, i_point) + &
          dcmplx( dot_product(wave(1:n_compute_c,i_point), work(1:n_compute_c,i_point)), 0.0d0 )  
    end do
  enddo !i_coord 

 !------------- (2.2) begin for density_matrix_compute terms------------------------
 !                    \sum_{u} gradient_wave^{(1)}(u,r) [\sum_{v} DM_{uv} wave(v,r) ]    
 !                   +\sum_{u} gradient_wave(u,r)       [\sum_{v} DM_{uv} wave^{(1)}(v,r) ]
 !                   +\sum_{u} wave^{(1)}(u,r) [\sum_{v} DM_{uv} gradient_wave(v,r) ]
 !                   +\sum_{u} wave(u,r)       [\sum_{v} DM_{uv} gradient_wave^{(1)}(v,r) ] 
 !--------------------end  for density_matrix_compute terms----------------------------
  do i_coord = 1 ,3 

    work = 0.0d0
    call dgemm("N","N", n_compute_c, n_points, n_compute_c,&  !transa, transb, l, n, m,
         1.0d0,density_matrix_compute, n_compute_c, wave(1:n_compute_c,1:n_points), n_compute_c,& !alpha, a, lda, b, ldb 
         0.0d0,work,n_compute_c)  ! beta, c, ldc
    do i_point = 1, n_points,1
        first_order_gradient_rho(i_coord,i_point) = & 
        first_order_gradient_rho(i_coord,i_point) + &
        dcmplx( dot_product(first_order_gradient_wave(1:n_compute_c,i_coord,i_point), work(1:n_compute_c,i_point)), 0.0d0 )  
    end do

    work = 0.0d0
    call dgemm("N","N", n_compute_c, n_points, n_compute_c,&  !transa, transb, l, n, m,
         1.0d0,density_matrix_compute, n_compute_c, first_order_wave, n_compute_c,& !alpha, a, lda, b, ldb 
         0.0d0,work,n_compute_c)  ! beta, c, ldc
 
    do i_point = 1, n_points,1
        first_order_gradient_rho(i_coord, i_point) = & 
        first_order_gradient_rho(i_coord, i_point) + &
        dcmplx( dot_product(gradient_basis_wave(1:n_compute_c,i_coord,i_point), work(1:n_compute_c,i_point)), 0.0d0 )  
    end do
 
    work = 0.0d0
    call dgemm("N","N", n_compute_c, n_points, n_compute_c,&  !transa, transb, l, n, m,
         1.0d0,density_matrix_compute, n_compute_c, & 
         gradient_basis_wave(1:n_compute_c,i_coord,1:n_points), n_compute_c,& !alpha, a, lda, b, ldb 
         0.0d0, work,n_compute_c)  ! beta, c, ldc
    do i_point = 1, n_points,1
        first_order_gradient_rho(i_coord, i_point) = &
        first_order_gradient_rho(i_coord, i_point) + &
        dcmplx( dot_product(first_order_wave(1:n_compute_c,i_point), work(1:n_compute_c,i_point)), 0.0d0 ) 
    end do
 
    work = 0.0d0
    call dgemm("N","N", n_compute_c, n_points, n_compute_c,&  !transa, transb, l, n, m,
         1.0d0,density_matrix_compute, n_compute_c, & 
         first_order_gradient_wave(1:n_compute_c,i_coord,1:n_points), n_compute_c,& !alpha, a, lda, b, ldb 
         0.0d0, work,n_compute_c)  ! beta, c, ldc
    do i_point = 1, n_points,1
        first_order_gradient_rho(i_coord, i_point) = &
        first_order_gradient_rho(i_coord, i_point) + &
          dcmplx( dot_product(wave(1:n_compute_c,i_point), work(1:n_compute_c,i_point)), 0.0d0 )    
    end do
 
  enddo ! i_coord 


 
  endif !n_compute_c 
 
 
  deallocate(work)
  deallocate(first_order_density_matrix_compute_real)
  deallocate(first_order_wave)
  deallocate(first_order_gradient_wave)
 

end subroutine evaluate_first_order_gradient_rho_phonon_reduce_memory
!---------------------------------------------------------------------
!******	 
