!****s* FHI-aims/evaluate_first_order_H_p1
!  NAME
!    evaluate_first_order_H_p1
!  SYNOPSIS

subroutine evaluate_first_order_H_p1( & 
           n_points, partition_tab, & 
           n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave, & 
           H_times_psi, &  
           first_order_rho,first_order_potential,dVxc_drho, & 
           first_order_H_sparse & 
           ,test_T,test_V_hartree,test_V_xc)

!  PURPOSE
!<1>  calculate the first-order Hamiltion matrix elements for phonon_gamma
!     four terms:

!   (1) <X0| free_V_hartree(1)  |X0>  }
!   (2) <X0| delta_V_hartree(1) |X0>  } ===> Hellman-Feynman term in this subroutine 
!   (2) <X0| dVxc/drho * rho(1) |X0>  }
 
!   (3) <X0| Hks(0)             |X1>  } ===> Pulay term in this subroutine
!   (4) <X1| Hks(0)             |X0>  }
!
!                                     ------ shanghui  2013.12.30

!<2>  change to sparse format for H.  ------ shanghui  2014.12.06

!  USES

  use dimensions
  use species_data
  use geometry
  use runtime_choices
  use basis 
  use pbc_lists

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) :: partition_tab

  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: H_times_psi

  real*8, dimension(3, n_centers_in_sc_DFPT, n_points), intent(in) :: first_order_rho
  real*8, dimension(3, n_centers_in_sc_DFPT, n_points), intent(in) :: first_order_potential
  real*8, dimension(n_points), intent(in) :: dVxc_drho


  real*8, dimension(3, n_centers_in_sc_DFPT, n_hamiltonian_matrix_size), intent(inout) ::  & 
                                                                        first_order_H_sparse

  logical, intent(in) :: test_T 
  logical, intent(in) :: test_V_hartree 
  logical, intent(in) :: test_V_xc 
!  INPUTS
!  o  n_points -- number of grid points in this grid batch
!  o  partition_tab -- values of partition function
!  o  n_compute_c -- basis set computed at this subroutine
!  o  i_basis_index -- i_compute to i_basis
!  o  wave -- wave function in unit cell
!  o  gradient_basis_wave -- gradient wave function in unit cell
!  o  H_times_psi           -- H*psi  
!  o  first_order_rho       -- rho_tot(1)
!  o  first_order_potential -- V_free_hartree(1)+delta_V_hartree(1)
!  o  dVxc_drho             -- d Vxc/ drho
!
!  OUTPUT
!  first_order_H --- PM_index


  integer :: i_coord, i_point, i_place
  integer :: i_compute,j_compute
  integer :: i_basis,j_basis
  integer :: i_basis_uc, j_basis_uc
  integer :: i_atom,j_atom
  integer :: i_cell_in_hamiltonian,j_cell_in_hamiltonian
  integer :: i_cell_in_sc_DFPT,j_cell_in_sc_DFPT

  integer :: i_center_trans, j_center_trans

  integer :: k_center_in_sc_DFPT, k_cell_in_sc_DFPT, k_cell_trans, k_atom, k_center_trans 

  real*8  :: point_term
 


  do i_compute=1,n_compute_c,1
     i_basis    = i_basis_index(i_compute)
     i_basis_uc = Cbasis_to_basis(i_basis)
     i_atom     = Cbasis_to_atom(i_basis)
     i_cell_in_hamiltonian     = center_to_cell(Cbasis_to_center(i_basis))
     i_cell_in_sc_DFPT         = center_in_sc_DFPT_to_cell_in_sc_DFPT( &
                                 center_to_center_in_sc_DFPT(Cbasis_to_center(i_basis)) )


  do j_compute=1,n_compute_c,1
     j_basis    = i_basis_index(j_compute)
     j_basis_uc = Cbasis_to_basis(j_basis)


     if(j_basis_uc <= i_basis_uc) then
       j_atom = Cbasis_to_atom(j_basis)
       j_cell_in_hamiltonian = center_to_cell(Cbasis_to_center(j_basis))
       j_cell_in_sc_DFPT     = center_in_sc_DFPT_to_cell_in_sc_DFPT( &
                               center_to_center_in_sc_DFPT(Cbasis_to_center(j_basis)) )


     do i_place = &
        index_hamiltonian(1,position_in_hamiltonian(i_cell_in_hamiltonian,j_cell_in_hamiltonian), i_basis_uc),  & 
        index_hamiltonian(2,position_in_hamiltonian(i_cell_in_hamiltonian,j_cell_in_hamiltonian), i_basis_uc)  

     if( column_index_hamiltonian( i_place) == j_basis_uc)then

         !----------get the i_place related basis, and basis related center_trans-----
          i_center_trans = cell_and_atom_to_center_sc_DFPT( &
                           cell_diff_sc_DFPT(i_cell_in_sc_DFPT,j_cell_in_sc_DFPT),i_atom)

          j_center_trans = j_atom
         !---------end get the i_place realated basis--------------------------------



      do i_point = 1, n_points, 1
         point_term = partition_tab(i_point)    & 
                    * wave(i_compute,i_point) * wave(j_compute,i_point)   

!-------------------(1) Hellman-Feynman term---------------------
         do k_center_in_sc_DFPT = 1, n_centers_in_sc_DFPT
            k_cell_in_sc_DFPT = center_in_sc_DFPT_to_cell_in_sc_DFPT(k_center_in_sc_DFPT)
            k_atom = center_in_sc_DFPT_to_atom(k_center_in_sc_DFPT)

            k_cell_trans = cell_diff_sc_DFPT(k_cell_in_sc_DFPT,j_cell_in_sc_DFPT)
            k_center_trans = cell_and_atom_to_center_sc_DFPT(k_cell_trans, k_atom)
 
            do i_coord = 1,3

!           if(test_T) then 
!            first_order_H_sparse(i_coord,k_center_trans, i_place) = &
!            first_order_H_sparse(i_coord,k_center_trans, i_place) + & 
!                   0.0d0
!           else if(test_V_hartree) then 
!            first_order_H_sparse(i_coord,k_center_trans, i_place) = &
!            first_order_H_sparse(i_coord,k_center_trans, i_place) + & 
!                   point_term*first_order_potential(i_coord, k_center_in_sc_DFPT, i_point) 
!           else if(test_V_xc) then
!            first_order_H_sparse(i_coord,k_center_trans, i_place) = &
!            first_order_H_sparse(i_coord,k_center_trans, i_place) + & 
!                   point_term*dVxc_drho(i_point)*first_order_rho(i_coord,k_center_in_sc_DFPT, i_point) !(3) 
!           else ! all
            first_order_H_sparse(i_coord,k_center_trans, i_place) = &
            first_order_H_sparse(i_coord,k_center_trans, i_place) + & 
                   point_term*first_order_potential(i_coord, k_center_in_sc_DFPT, i_point)+ & 
                   point_term*dVxc_drho(i_point)*first_order_rho(i_coord,k_center_in_sc_DFPT, i_point) !(3) 
!           endif

            enddo
         enddo
 
!-------------------(2) Pulay term-------------------------------
            do i_coord = 1,3
           first_order_H_sparse(i_coord, i_center_trans,i_place) = & 
           first_order_H_sparse(i_coord, i_center_trans,i_place) - & 
                     partition_tab(i_point)                              &  
                     *gradient_basis_wave(i_compute,i_coord,i_point)     &  
                     *H_times_psi(j_compute,i_point)   
 
 
           first_order_H_sparse(i_coord, j_center_trans,i_place) = & 
           first_order_H_sparse(i_coord, j_center_trans,i_place) - & 
                     partition_tab(i_point)                              &  
                     *gradient_basis_wave(j_compute,i_coord,i_point)     &  
                     *H_times_psi(i_compute,i_point)   
            enddo
      enddo ! i_point
         
      endif ! column

      enddo ! i_place 
      endif ! j_basis_uc <i_basis_uc


  enddo  ! j_compute
  enddo  ! i_compute 



end subroutine evaluate_first_order_H_p1
!******
