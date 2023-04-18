!****s* FHI-aims/get_angular_grid
!  NAME
!    get_angular_grid
!  SYNOPSIS

subroutine get_momentummatrix &
     (KS_eigen, KS_vec, KS_vec_complex, occ_numbers, chemical_potential, partition_tab,&
       l_shell_max)

!  PURPOSE
!
!  Wrapper function for calculating and outputting of Momentummatrix elements
!
!  USES
  use calculate_mommat_base, only:moment_one,moment_two,moment_three,&
                                  get_state_minmax_k,allocate_mommat,allocate_mommat_k,&
                                  allocate_moment_one,allocate_moment_two,allocate_moment_three,&
                                  calculate_mommat_p0,calc_moment_p0,clean_mommat,clean_mommat_final,&
                                  mommat_full_oned_up,mommat_full_oned_low,mommat_full_w_up,&
                                  mommat_full_w_low,mommat_full_w_complex_up,&
                                  mommat_full_w_complex_low,out_moment,work_ovl_mom
  use dimensions,only:n_states,n_spin,n_k_points,n_basis,n_k_points_task,n_species,n_full_points
  use runtime_choices,only:k_point_mom
  use localorb_io,only:localorb_info
  use mpi_tasks,only:n_tasks,myid,mpi_comm_global
  use hdf5_tools,only:HID_T,open_hdf5,open_hdf5_kart,out_bands_kart,out_k_points_kart,&
                      outmetafile,out_moment_kart,close_hdf5
  implicit none
!  ARGUMENTS

  real*8 , dimension(n_states, n_spin, n_k_points), INTENT(IN) :: KS_eigen
  complex*16, dimension(n_basis, n_states, n_spin, n_k_points_task), INTENT(IN)::  &
                                                                 KS_vec_complex
  real*8, dimension(n_basis, n_states, n_spin, n_k_points_task), INTENT(IN)::  KS_vec
  real*8, dimension(n_states, n_spin,n_k_points), INTENT(IN) :: occ_numbers
  real*8, INTENT(IN) :: chemical_potential
  real*8, target, dimension(n_full_points), INTENT(IN) :: partition_tab
  integer, dimension(n_species), INTENT(IN) :: l_shell_max 



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

  INTEGER(HID_T) :: file_id                  ! File identifier
  INTEGER(HID_T) :: plist_id                 ! Property list identifier
  character*128 :: info_str
  integer :: info
  integer :: n_state_min_in
  integer :: n_state_max_in
  !  counters

  integer :: i_k
  integer :: new_k_point
  !  begin work

    write(info_str,'(6X,A,1X,I4)') "Momentum Matrix post processing starts"
    call localorb_info ( info_str )
    call get_state_minmax_k(KS_eigen, n_state_min_in, n_state_max_in)
    if(k_point_mom==0)then	
	call open_hdf5('mommat.h5', file_id, plist_id)
	call open_hdf5_kart("Momentummatrix", file_id, plist_id,n_state_min_in,&
                            n_state_max_in)
	call mpi_barrier(mpi_comm_global,info)
	call out_bands_kart(KS_eigen,file_id,plist_id,n_state_min_in,&
                            n_state_max_in)
	call out_k_points_kart(file_id,plist_id)
	call outmetafile(KS_eigen,occ_numbers,chemical_potential,file_id,&
                         plist_id,n_state_min_in,n_state_max_in)
   endif
   call allocate_mommat()
   call allocate_mommat_k()
   !if (k_point_mom==0)then
   !   call allocate_moment_one(n_state_min_in, n_state_max_in)
   !else
      call allocate_moment_one(n_state_min_in, n_state_max_in)
      call allocate_moment_two(n_state_min_in, n_state_max_in)
      call allocate_moment_three(n_state_min_in, n_state_max_in)
   !endif
   moment_one = (0d0, 0d0)
   moment_two = (0d0, 0d0)
   moment_three = (0d0, 0d0)
   call calculate_mommat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up,&
                              mommat_full_oned_low,1)
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
                           new_k_point, 1, n_state_min_in, n_state_max_in)
	  if (k_point_mom==0)then	
             call out_moment_kart(moment_one,KS_eigen,new_k_point,&
                                  "Momentummatrix", file_id, plist_id,1,&
                                  n_state_min_in,n_state_max_in)
          endif
	endif
      endif
   enddo
   call calculate_mommat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up, &
                              mommat_full_oned_low,2)
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
	  call calc_moment_p0(moment_two,mommat_full_w_up, mommat_full_w_low, &
                           mommat_full_w_complex_up, mommat_full_w_complex_low,&
                            KS_vec(:,:,:,i_k) ,KS_vec_complex(:,:,:,i_k) , &
                            new_k_point, 2, n_state_min_in, n_state_max_in)
	  if (k_point_mom==0)then
            call out_moment_kart(moment_two,KS_eigen,new_k_point,&
                                 "Momentummatrix", file_id, plist_id,2, &
                                  n_state_min_in,n_state_max_in)
          endif
	endif
      endif
    enddo
    call calculate_mommat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up,&
                               mommat_full_oned_low,3)
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

	  call calc_moment_p0(moment_three,mommat_full_w_up, mommat_full_w_low,&
                           mommat_full_w_complex_up, mommat_full_w_complex_low,&
                           KS_vec(1:n_basis,1:n_states,1:n_spin,i_k) ,&
                           KS_vec_complex(1:n_basis,1:n_states,1:n_spin,i_k) , &
                           new_k_point, 3, n_state_min_in, n_state_max_in)

	  if (k_point_mom==0)then 
            call out_moment_kart(moment_three,KS_eigen,new_k_point,&
                                 "Momentummatrix", file_id, plist_id,3, &
                                 n_state_min_in,n_state_max_in)
          endif
	endif
      endif
    enddo
    if (k_point_mom/=0.and.myid==modulo(k_point_mom, n_tasks))then 
         call out_moment(moment_one,moment_two,moment_three, KS_eigen, &
                       k_point_mom,n_state_min_in,n_state_max_in)
    endif
    call clean_mommat()
    call clean_mommat_final()
    if(k_point_mom==0)then
      call mpi_barrier(mpi_comm_global,info)
      call close_hdf5(file_id,plist_id)
    endif
    write(info_str,'(6X,A,1X,I4)') "Momentum Matrix post processing finished"
 
end subroutine get_momentummatrix
!******	
