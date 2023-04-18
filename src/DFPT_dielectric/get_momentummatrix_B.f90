!****s* FHI-aims/get_angular_grid
!  NAME
!    get_angular_grid
!  SYNOPSIS

subroutine get_momentummatrix_B &
     (KS_eigen, KS_vec, KS_vec_complex, occ_numbers, chemical_potential, partition_tab,&
       l_shell_max, j_coord, momentum_matrix)

!  PURPOSE
!
!  Wrapper function for calculating and outputting of Momentummatrix elements
!
!  USES
  use calculate_mommat_base
  use dimensions
  use runtime_choices
  use localorb_io
  use mpi_utilities
!  ARGUMENTS

  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) :: KS_eigen
  complex*16, dimension(n_basis, n_states, n_spin, n_k_points), INTENT(IN)::  &
                                                                 KS_vec_complex
  real*8, dimension(n_basis, n_states, n_spin, n_k_points), INTENT(IN)::  KS_vec
  real*8, dimension(n_states, n_spin,n_k_points), INTENT(IN) :: occ_numbers
  real*8, INTENT(IN) :: chemical_potential
  real*8, target, dimension(n_full_points), INTENT(IN) :: partition_tab
  integer, dimension(n_species), INTENT(IN) :: l_shell_max 

  integer, INTENT(IN) :: j_coord

  complex*16, INTENT(OUT) ::  & 
    momentum_matrix(n_states,n_states,n_k_points_task)


!  INPUTS
!   o KS_eigen
!   o KS_vec/KS_vec_complex
!   o occ_numbers
!   o chemical_potential
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
  integer :: n_state_min_in
  integer :: n_state_max_in
  integer :: i_state, j_state

  !  counters

  integer :: i_k
  integer :: new_k_point
  !  begin work

    write(info_str,'(6X,A,1X,I4)') "Momentum Matrix post processing starts"
    call localorb_info ( info_str )
    call get_state_minmax_k(KS_eigen, n_state_min_in, n_state_max_in)

    write(use_unit,*) 'n_state_min_in, n_state_max_in:', n_state_min_in, n_state_max_in
    write(use_unit,*) 'j_coord:', j_coord

   call allocate_mommat()
   call allocate_mommat_k()
   !if (k_point_mom==0)then
   !   call allocate_moment_one(n_state_min_in, n_state_max_in)
   !else
      call allocate_moment_one(n_state_min_in, n_state_max_in)
      call allocate_moment_two(n_state_min_in, n_state_max_in)
      call allocate_moment_three(n_state_min_in, n_state_max_in)
   !endif

   call calculate_mommat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up,&
                              mommat_full_oned_low,j_coord)
   i_k = 0
   do new_k_point = 1,n_k_points
      if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points) then
        ! k-point stored on current mpi task
        i_k = i_k+1 !new_k_point
        if (new_k_point==k_point_mom.or.k_point_mom==0)then			
	  call construct_overlap( mommat_full_oned_up, mommat_full_w_up, &
                                  mommat_full_w_complex_up, new_k_point, &
                                  work_ovl_mom )
	  call construct_overlap( mommat_full_oned_low, mommat_full_w_low, &
                                  mommat_full_w_complex_low, new_k_point, &
                                  work_ovl_mom )
	  call calc_moment_p0(moment_one,mommat_full_w_up, mommat_full_w_low, &
                           mommat_full_w_complex_up, mommat_full_w_complex_low,&
                           KS_vec(:,:,:,i_k) ,KS_vec_complex(:,:,:,i_k) , &
                           new_k_point, j_coord, n_state_min_in, n_state_max_in)





!        num = 0 
!        do i_state = 1, n_states 
!        do j_state = i_state, n_states
!           num = num + 1
!           momentum_matrix(i_state,j_state,i_k) = & 
!           moment_one(num)
!
!           !momentum_matrix(j_state,i_state,i_k) = & 
!           !conjg(moment_onenum)
!        enddo 
!        enddo 


	endif
      endif
   enddo

!  call calculate_mommat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up, &
!                             mommat_full_oned_low,2)
!  i_k = 0
!  do new_k_point = 1,n_k_points
!    if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points) then
!       ! k-point stored on current mpi task
!       i_k = i_k+1 !new_k_point
!       if (new_k_point==k_point_mom.or.k_point_mom==0)then
!         call construct_overlap( mommat_full_oned_up, mommat_full_w_up, &
!                                 mommat_full_w_complex_up, new_k_point, &
!                                 work_ovl_mom )
!         call construct_overlap( mommat_full_oned_low, mommat_full_w_low, &
!                                 mommat_full_w_complex_low, new_k_point, &
!                                 work_ovl_mom )
!         call calc_moment_p0(moment_two,mommat_full_w_up, mommat_full_w_low, &
!                          mommat_full_w_complex_up, mommat_full_w_complex_low,&
!                           KS_vec(:,:,:,i_k) ,KS_vec_complex(:,:,:,i_k) , &
!                           new_k_point, 2, n_state_min_in, n_state_max_in)
!       endif
!     endif
!   enddo
!   write(use_unit,*) '2'

!  call calculate_mommat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up,&
!                             mommat_full_oned_low,3)
!  i_k = 0
!  do new_k_point = 1,n_k_points
!    if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points) then
!      ! k-point stored on current mpi task
!      i_k = i_k+1 !new_k_point
!      if (new_k_point==k_point_mom.or.k_point_mom==0)then
!        call construct_overlap( mommat_full_oned_up, mommat_full_w_up, &
!                                mommat_full_w_complex_up, new_k_point, &
!                                work_ovl_mom )
!        call construct_overlap( mommat_full_oned_low, mommat_full_w_low, &
!                                mommat_full_w_complex_low, new_k_point, &
!                                work_ovl_mom )
!        call calc_moment_p0(moment_three,mommat_full_w_up, mommat_full_w_low,&
!                         mommat_full_w_complex_up, mommat_full_w_complex_low,&
!                         KS_vec(:,:,:,i_k) ,KS_vec_complex(:,:,:,i_k) , &
!                         new_k_point, 3, n_state_min_in, n_state_max_in)
!      endif
!    endif
!  enddo
!  write(use_unit,*) '3'
!
   
    call clean_mommat()
    call clean_mommat_final()
   !if (allocated(mommat_full_oned_up))  write(use_unit,*) 'haha' !deallocate(mommat_full_oned_up)
   ! deallocate(mommat_full_oned_up) 
   !if (allocated(mommat_full_oned_low)) write(use_unit,*) 'haha1' !deallocate(mommat_full_oned_low)
   ! deallocate(mommat_full_oned_low)

   !write(use_unit,*) 'test'
  !  write(info_str,'(6X,A,1X,I4)') "Momentum Matrix post processing finished"
 
end subroutine get_momentummatrix_B
!******	
