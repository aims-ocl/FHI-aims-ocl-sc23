subroutine get_momentum_matrix &
     (KS_eigen, KS_vec, KS_vec_complex, occ_numbers, &
      partition_tab, l_shell_max, i_coord, momentum_matrix)

!  PURPOSE
!
!  Wrapper function for calculating momentum_matrix 
!  1.function 'calculate_mommat_p0', real space integral of atom centered basis 
!    function i with the gradient of function j for in each cell (cell distance)
!  2.function 'construct_overlap' 'fourier transform' of real space integrals by
!    summing up the phases resulting from cell distances
!    funtion is called once for the upper triangle and once for the lower
!    triangle of the sparse matrix from first step
!  3.function 'calc_moment_p0' basis transformation from atom centered basis to 
!    KS-orbitals

!  USES
  use calculate_mommat_base
  use dimensions
  use runtime_choices
  use localorb_io
  use mpi_utilities

  implicit none

  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) :: KS_eigen
  real*8, dimension(n_basis, n_states, n_spin, n_k_points_task), INTENT(IN)::  KS_vec
  complex*16, dimension(n_basis, n_states, n_spin, n_k_points_task), INTENT(IN)::  &
                                                                 KS_vec_complex
  real*8, dimension(n_states, n_spin,n_k_points), INTENT(IN) :: occ_numbers
  real*8, target, dimension(n_full_points), INTENT(IN) :: partition_tab
  integer, dimension(n_species), INTENT(IN) :: l_shell_max 
  integer, INTENT(IN) :: i_coord

  complex*16, INTENT(OUT) ::  & 
    momentum_matrix(n_states,n_states,n_k_points_task)

!  INPUTS
!   o KS_eigen
!   o KS_vec/KS_vec_complex
!   o occ_numbers
!   o partition_tab
!   o l_shell_max
!  OUTPUT
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



  !  local variables
  character*128 :: info_str
  integer :: info

  !  counters
  integer :: i_k_point, i_k_task, i_place, i_state, j_state,i_basis_1,i_basis_2
  integer :: num_base
  complex*16 ::  momentum_matrix_temp(((n_states+1)*(n_states)/2)*((n_spin*(n_spin+1))/2))
  complex*16 ::  momentum_matrix_complex(n_basis,n_basis) 

  !  begin work
    write(info_str,'(6X,A,1X,I4)') "get_momentum_matrix begin"
    call localorb_info ( info_str )
 
      call allocate_mommat()
      call allocate_mommat_k()

      call calculate_mommat_p0( partition_tab, l_shell_max, &
                                 mommat_full_oned_up, mommat_full_oned_low,i_coord)

      i_k_task = 0
      do i_k_point = 1,n_k_points, 1
         if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
          i_k_task = i_k_task + 1
          call construct_overlap( mommat_full_oned_up, mommat_full_w_up,&
                                  mommat_full_w_complex_up, i_k_point,&
                                  work_ovl_mom )
          call construct_overlap( mommat_full_oned_low, mommat_full_w_low,  & 
                                  mommat_full_w_complex_low, i_k_point,&
                                  work_ovl_mom )
          num_base = 0
          do i_basis_2 = 1, n_basis
             do i_basis_1 = 1, i_basis_2
                num_base=num_base+1
                momentum_matrix_complex(i_basis_1, i_basis_2)=mommat_full_w_complex_up(num_base)
                momentum_matrix_complex(i_basis_2, i_basis_1)=conjg(mommat_full_w_complex_low(num_base))
             enddo
          enddo
         
          write(use_unit,*) 'ik:',i_k_point
          do i_basis_1 = 1, n_basis 
             write(use_unit,*) momentum_matrix_complex(i_basis_1,1:n_basis) 
          enddo     
          write(use_unit,*)  '  ' 

          call calc_moment_p0(momentum_matrix_temp,mommat_full_w_up, mommat_full_w_low, & 
                              mommat_full_w_complex_up, mommat_full_w_complex_low, & 
                              KS_vec(:,:,:,i_k_task), KS_vec_complex(:,:,:,i_k_task), & 
                              i_k_point, i_coord, 1, n_states)


          i_place = 0 
          do i_state = 1, n_states 
          do j_state = i_state, n_states
             i_place = i_place + 1
             momentum_matrix(i_state,j_state,i_k_task) = & 
             momentum_matrix_temp(i_place)

             momentum_matrix(j_state,i_state,i_k_task) = & 
             conjg(momentum_matrix_temp(i_place))
          enddo 
          enddo 
           write(use_unit,*) 'ik: momentum_matrix_MO:',i_k_point,momentum_matrix(1,2,i_k_task) 
          endif
      enddo

      call clean_mommat()
      call clean_mommat_final()

    write(info_str,'(6X,A,1X,I4)') "get_momentum_matrix end"
    call localorb_info ( info_str )
 
    
    write(use_unit,*) momentum_matrix(1,1:2,1) 

end subroutine get_momentum_matrix

