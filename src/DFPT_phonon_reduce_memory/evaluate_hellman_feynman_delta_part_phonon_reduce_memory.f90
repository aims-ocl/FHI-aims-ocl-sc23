!****s* FHI-aims/evaluate_hellman_feynman_delta_part_phonon_reduce_memory
!  NAME
!    evaluate_hellman_feynman_delta_part_phonon_reduce_memory
!  SYNOPSIS

subroutine evaluate_hellman_feynman_delta_part_phonon_reduce_memory(  & 
           n_points, partition_tab, &
           rho, first_order_rho, & 
           coords_npoints,  &
           i_q_point, j_atom, j_coord,  & 
           hellman_feynman_dynamical_matrix_delta_part )



!  PURPOSE

!  hellman_feynman_delta_part for phonon_reduce_memory 
!  using Ewald method.
!  shanghui, 2015.07.30
!  USES

  use constants, only: sqrt_pi
  use runtime_choices, only: Ewald_radius
  use dimensions
  use pbc_lists
  use arch_specific, only : Arch_erfc
  use geometry, only: recip_lattice_vector, species
  use localorb_io, only: use_unit

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) ::  partition_tab

  real*8, dimension(n_points), intent(in) :: rho
  complex*16, dimension(n_points), intent(in) :: first_order_rho
  real*8, dimension(3,n_points), intent(in) :: coords_npoints

  integer , intent(in) :: i_q_point
  integer , intent(in) :: j_atom
  integer , intent(in) :: j_coord

  complex*16, dimension(3,n_atoms,3,n_atoms), intent(inout) :: hellman_feynman_dynamical_matrix_delta_part

!  INPUTS
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!
!  OUTPUT
!  hellman_feynman_dynamical_matrix_delta_part

  integer :: i_point, i_coord, i_atom,i_cell, n, i_center,current_center
  real*8  :: Gr(3), F_table_1,r, nx, ny, nz, k(3), k2 
  complex*16 :: exp_iqr(3) 
  real*8, external :: ddot
  real*8 ::  Ewald_radius_my 

  complex*16 :: Ewald_term_real(3,n_atoms), Ewald_term_reciprocal

  Ewald_radius_my = Ewald_radius

  do i_point = 1,n_points 

       Gr(1) = ddot(3,recip_lattice_vector(1:3,1),1, coords_npoints(1:3,i_point),1)
       Gr(2) = ddot(3,recip_lattice_vector(1:3,2),1, coords_npoints(1:3,i_point),1)
       Gr(3) = ddot(3,recip_lattice_vector(1:3,3),1, coords_npoints(1:3,i_point),1)

       exp_iqr(1) = exp((0,1)*Gr(1)*k_point_list(i_q_point,1))
       exp_iqr(2) = exp((0,1)*Gr(2)*k_point_list(i_q_point,2))
       exp_iqr(3) = exp((0,1)*Gr(3)*k_point_list(i_q_point,3))


     !-----------real part------------------------------
       Ewald_term_real(1:3,1:n_atoms) = (0.0d0,0.0d0) 
       do i_coord = 1,3

       do i_center = 1,  n_centers_in_hamiltonian
           current_center = centers_in_hamiltonian(i_center)
           i_cell = center_to_cell(current_center)
           i_atom = center_to_atom(current_center)

           r=dsqrt( (coords_npoints(1,i_point)-coords_center(1,current_center))**2 + & 
                    (coords_npoints(2,i_point)-coords_center(2,current_center))**2 + & 
                    (coords_npoints(3,i_point)-coords_center(3,current_center))**2) 
           if(r.lt.1.0d-6) then  
           write(use_unit,*) 'shanghui in evaluate_hellman_feynman_delta_part_phonon_reduce_memory: r is too small:', r
           endif  

           !if(r.gt.1.0d-5) then  !i_center.ne.i_atom) then 
           F_table_1 =  - (( (2.0d0 * exp(- (r**2 /Ewald_radius_my**2 )))/(Ewald_radius_my * sqrt_pi ) & 
                           + Arch_erfc(r /Ewald_radius_my )/r ) /r**2 )

           Ewald_term_real(i_coord,i_atom) =  Ewald_term_real(i_coord, i_atom) + & 
           species(i_atom) * F_table_1 * & 
           (coords_center(i_coord,current_center)-coords_npoints(i_coord,i_point))  * & 
            dconjg(k_phase(i_cell,i_q_point)) * exp_iqr(1)*exp_iqr(2)*exp_iqr(3)
           !else 
           !write(use_unit,*) 'shanghui in evaluate_hellman_feynman_delta_part_phonon_reduce_memory: r is too small:', r
           !F_table_1 = 0.0d0
           !endif

!           if(i_coord.eq.1.and.i_atom.eq.1.and.dabs(F_table_1).gt.10) then 
!           write(use_unit,*) 'myid,icenter:',myid,i_center
!           write(use_unit,*) 'F_tabe, Ewald:',F_table_1 ,Ewald_term_real(i_coord,i_atom)     
!           write(use_unit,*) 'r,ZI,R-r',r, species(i_atom), coords_center(i_coord,current_center)-coords_npoints(i_coord,i_point) 
!           write(use_unit,*) '1:',(2.0d0 * exp(- (r**2 /Ewald_radius_my**2 )))/(Ewald_radius_my * sqrt_pi )
!           write(use_unit,*) '2:',Arch_erfc(r /Ewald_radius_my )/r 
! 
!           call aims_stop
!           endif

        enddo ! i_center 
        enddo ! i_coord
  
        !-----------reciprocal part------------------------------
!        do i_coord = 1, 3
!        do i_atom = 1,n_atoms
!        
!          Ewald_term_reciprocal = (0.0d0, 0.0d0)
!          do n = 1, n_k_points_hartree
!          ! K = G-q
!          nx = k_points_hartree(n, 1) - k_point_list(i_q_point,1)
!          ny = k_points_hartree(n, 2) - k_point_list(i_q_point,2)
!          nz = k_points_hartree(n, 3) - k_point_list(i_q_point,3)
! 
!          k(1) = nx * recip_lattice_vector(1,1) + ny * recip_lattice_vector(1,2) +  nz * recip_lattice_vector(1,3)
!          k(2) = nx * recip_lattice_vector(2,1) + ny * recip_lattice_vector(2,2) +  nz * recip_lattice_vector(2,3)
!          k(3) = nx * recip_lattice_vector(3,1) + ny * recip_lattice_vector(3,2) +  nz * recip_lattice_vector(3,3)
! 
!          k2 = k(1)**2 + k(2)**2 + k(3)**2
! 
!          Ewald_term_reciprocal = Ewald_term_reciprocal + & 
!             species(i_atom) * (-1) * (0,1) * k(i_coord) *  & 
!             4*pi / cell_volume /k2 * exp(-Ewald_radius_my**2 * k2 / 4.0d0) * &
!             exp ( (0,1)*(k(i_coord)* (coords_npoints(i_coord,i_point) - coords_center(i_coord, i_atom))   ) )
!          enddo  ! i_k_point_hartree    
! 
!        hellman_feynman_dynamical_matrix_delta_part(i_coord,i_atom,j_coord,j_atom) = & 
!        hellman_feynman_dynamical_matrix_delta_part(i_coord,i_atom,j_coord,j_atom) - &    
!        partition_tab(i_point)* first_order_rho(i_point) * &  
!        !Ewald_term_reciprocal
!        (Ewald_term_real(i_coord,i_atom) + Ewald_term_reciprocal) 
! 
!        enddo ! i_atom 
!        enddo ! i_coord 

  enddo ! i_point

end subroutine evaluate_hellman_feynman_delta_part_phonon_reduce_memory
!******
