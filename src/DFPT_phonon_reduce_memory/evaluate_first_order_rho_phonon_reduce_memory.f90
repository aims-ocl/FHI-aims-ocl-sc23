!****s* FHI-aims/evaluate_first_order_rho_phonon_reduce_memory
!  NAME
!   evaluate_first_order_rho_phonon_reduce_memory
!  SYNOPSIS

subroutine evaluate_first_order_rho_phonon_reduce_memory( & 
           n_points,n_compute_c, i_basis_index, & 
           wave,gradient_basis_wave, &
           density_matrix_compute,first_order_density_matrix_compute,coords_npoints, & 
           i_q_point, j_atom, j_coord, &
           first_order_rho)

!  PURPOSE
!  Evaluates first_order_rho= rho_gradien+ first_order_DM*wave*wave
!  shanghui  2012.05.02


!  for phonon_reduce_memory 
!  shanghui 2015.07.30
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
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave

  real*8     :: density_matrix_compute(n_compute_c,n_compute_c)
  complex*16 :: first_order_density_matrix_compute(n_compute_c,n_compute_c)
  real*8, dimension(3,n_points), intent(in) :: coords_npoints
  integer , intent(in) :: i_q_point
  integer , intent(in) :: j_atom
  integer , intent(in) :: j_coord
 
  complex*16, dimension(n_points) :: first_order_rho


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

  
!------------------(1) This is the slow verion---------------------------
!  complex*16, allocatable ::  rho_gradient_atom(:) 
!  complex*16, allocatable ::  first_order_DM_rho(:) 
!  !-------begin allocatation-----------------
!  allocate(rho_gradient_atom(n_points))
!  allocate(first_order_DM_rho(n_points))
!  first_order_DM_rho(1:n_points) = (0.0d0, 0.0d0)
!  rho_gradient_atom(1:n_points)  = (0.0d0, 0.0d0)
!  first_order_rho(1:n_points)    = (0.0d0, 0.0d0)
! 
! 
! !--------------------first: rho_gradient_atom----------------------
!        do i_compute=1,n_compute_c
!        do j_compute=1,n_compute_c
!        
!           do i_point=1, n_points,1
! 
!           i_basis=i_basis_index(i_compute)
!           j_basis=i_basis_index(j_compute)          
! 
!           Gr(1) = ddot(3,recip_lattice_vector(1:3,1),1, coords_npoints(1:3,i_point),1)
!           Gr(2) = ddot(3,recip_lattice_vector(1:3,2),1, coords_npoints(1:3,i_point),1)
!           Gr(3) = ddot(3,recip_lattice_vector(1:3,3),1, coords_npoints(1:3,i_point),1)
!   
!           exp_iqr(1) = exp(-(0,1)*Gr(1)*k_point_list(i_q_point,1))
!           exp_iqr(2) = exp(-(0,1)*Gr(2)*k_point_list(i_q_point,2))
!           exp_iqr(3) = exp(-(0,1)*Gr(3)*k_point_list(i_q_point,3))
! 
!           if(j_atom.eq.Cbasis_to_atom(i_basis)) then
!           rho_gradient_atom(i_point)= &
!           rho_gradient_atom(i_point)+ & 
!           density_matrix_compute(i_compute, j_compute) * & 
!           (-gradient_basis_wave(i_compute,j_coord,i_point)*wave(j_compute,i_point))* & 
!            k_phase(center_to_cell(Cbasis_to_center( i_basis )),i_q_point)* & 
!            exp_iqr(1)*exp_iqr(2)*exp_iqr(3)   
!           endif
! 
!           if(j_atom.eq.Cbasis_to_atom(j_basis)) then
!           rho_gradient_atom(i_point)= &
!           rho_gradient_atom(i_point)+ & 
!           density_matrix_compute(i_compute, j_compute) * &
!           (-gradient_basis_wave(j_compute,j_coord,i_point)*wave(i_compute,i_point))* &
!            k_phase(center_to_cell(Cbasis_to_center( j_basis )),i_q_point)* &
!            exp_iqr(1)*exp_iqr(2)*exp_iqr(3)
!           endif
! 
!          enddo ! i_point
!  
!        enddo
!        enddo
! 
! !-------------second: first_order_rho---------------------------------
!        do i_compute=1,n_compute_c
!        do j_compute=1,n_compute_c
!        
!         do i_point=1, n_points,1
!         first_order_DM_rho( i_point)= &
!         first_order_DM_rho( i_point)+ & 
!         first_order_density_matrix_compute(i_compute, j_compute) * & 
!         wave(i_compute,i_point)*wave(j_compute,i_point)           
!         enddo !  i_point 
! 
!        enddo
!        enddo
! 
! 
!    
!        do i_point = 1,n_points
!           first_order_rho(i_point)    = &
!           rho_gradient_atom(i_point)  + & 
!           first_order_DM_rho(i_point)
!        enddo 
! 
! 
!  !-------begin deallocatation-----------------
!  deallocate(rho_gradient_atom)
!  deallocate(first_order_DM_rho)

!---------------------------(2) This is the fast verion---------------------------
!----here we only first consider i_q_point=0, further we will add more q point---- 
 
  real*8, allocatable ::  work(:,:) 
  real*8, allocatable ::  first_order_density_matrix_compute_real(:,:) 
  real*8, allocatable ::  first_order_wave(:,:) 
 

  allocate(work(n_compute_c, n_points))
  allocate(first_order_density_matrix_compute_real(n_compute_c,n_compute_c))
  allocate(first_order_wave(n_compute_c, n_points))

 
 
  if (n_compute_c.eq.0) then
     first_order_rho(1:n_points) = (0.0d0,0.0d0)
     return
  else
 
  first_order_rho(1:n_points)    = (0.0d0, 0.0d0)
  do i_compute = 1,n_compute_c 
    do j_compute = 1,n_compute_c 
       first_order_density_matrix_compute_real(i_compute,j_compute)  = & 
            dble(first_order_density_matrix_compute(i_compute,j_compute)) 
    enddo 
  enddo
 
 !-------------------------(2.0) wave^{(1)}(u,r) --------------------
   first_order_wave(1:n_compute_c,1:n_points) = 0.0d0
   do i_compute=1,n_compute_c
      i_basis=i_basis_index(i_compute)
      if(Cbasis_to_atom(i_basis).eq.j_atom) then 
        first_order_wave(i_compute,1:n_points) = -gradient_basis_wave(i_compute,j_coord,1:n_points)
      endif
   enddo
 

 
 !-------------------------(2.1) \sum_{u} wave(u,r) [\sum_{v} DM^{(1)}_{uv} wave(v,r) ]----------------
    work = 0.0d0
    call dgemm("N","N", n_compute_c, n_points, n_compute_c,&  !transa, transb, l, n, m,
         1.0d0,first_order_density_matrix_compute_real, n_compute_c, wave(1:n_compute_c,1:n_points), n_compute_c,& !alpha, a, lda, b, ldb 
         0.0d0, work,n_compute_c)  ! beta, c, ldc
 
    do i_point = 1, n_points,1
        first_order_rho(i_point) = &
          dcmplx( dot_product(wave(1:n_compute_c,i_point), work(1:n_compute_c,i_point)), 0.0d0 )  
    end do
 
 
 !-------------------------(2.2) \sum_{u} wave(u,r) [\sum_{v} DM_{uv} wave^{(1)}(v,r) ]----------------
    work = 0.0d0
    call dgemm("N","N", n_compute_c, n_points, n_compute_c,&  !transa, transb, l, n, m,
         1.0d0,density_matrix_compute, n_compute_c, first_order_wave, n_compute_c,& !alpha, a, lda, b, ldb 
         0.0d0,work,n_compute_c)  ! beta, c, ldc
 
    do i_point = 1, n_points,1
        first_order_rho(i_point) = first_order_rho(i_point) + &
          dcmplx( dot_product(wave(1:n_compute_c,i_point), work(1:n_compute_c,i_point)), 0.0d0 )  
    end do
 
 
 !-------------------------(2.3) \sum_{u} wave^{(1)}(u,r) [\sum_{v} DM_{uv} wave(v,r) ]----------------
    work = 0.0d0
    call dgemm("N","N", n_compute_c, n_points, n_compute_c,&  !transa, transb, l, n, m,
         1.0d0,density_matrix_compute, n_compute_c, wave(1:n_compute_c,1:n_points), n_compute_c,& !alpha, a, lda, b, ldb 
         0.0d0,work,n_compute_c)  ! beta, c, ldc
 
    do i_point = 1, n_points,1
        first_order_rho(i_point) = first_order_rho(i_point) + &
          dcmplx( dot_product(first_order_wave(1:n_compute_c,i_point), work(1:n_compute_c,i_point)), 0.0d0 )  
    end do
 
 
  endif !n_compute_c 
 
 
  deallocate(work)
  deallocate(first_order_density_matrix_compute_real)
  deallocate(first_order_wave)
 
 

end subroutine evaluate_first_order_rho_phonon_reduce_memory
!---------------------------------------------------------------------
!******	 
