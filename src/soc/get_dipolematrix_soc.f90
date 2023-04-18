!****s* FHI-aims/get_dipolematrix_soc
!  NAME
!    get_dipolematrix_soc
!  SYNOPSIS
subroutine get_dipolematrix_soc &
     (KS_eigen, occ_numbers, KS_vec_complex, chemical_potential, &
      partition_tab, l_shell_max)
!  PURPOSE
!  Wrapper function for calculating and outputting of Dipolematrix elements
!  USES
  use calculate_mommat_base, only : moment_one, moment_two, moment_three, &
                                    get_state_minmax_k, allocate_mommat, &
                                    allocate_mommat_k, allocate_moment_one_p1_soc, &
                                    allocate_moment_two_p1_soc, &
                                    allocate_moment_three_p1_soc, calculate_dipmat_p0, &
                                    calc_moment_p0_soc, clean_mommat, clean_mommat_final, &
                                    out_moment, mommat_full_oned_up, &
                                    mommat_full_oned_low, mommat_full_w_complex_low, &
                                    mommat_full_w_complex_up, mommat_full_w_low, &
                                    mommat_full_w_up, work_ovl_mom
  use dimensions,            only : n_k_points, n_k_points_task, n_species, &
                                    n_full_points, flag_out_boys, n_basis
  use dimensions_soc,        only : n_basis_soc, n_basis_soc_coll, n_basis_soc_ncoll, &
                                    n_saved_states_soc
  use runtime_choices,       only : k_point_mom
  use mpi_tasks,             only : mpi_comm_global, myid, n_tasks, aims_stop
  use localorb_io,           only : localorb_info
  use hdf5_tools,            only : HID_T, close_HDF5, out_moment_kart, &
                                    open_hdf5, open_hdf5_kart, out_bands_kart, &
                                    out_occ_kart, out_k_points_kart, outmetafile
  use soc_utilities,         only : convert_sr_to_soc_environment, revert_soc_to_sr_environment
  use boys,                  only : boys_centers
  implicit none
!  ARGUMENTS
  real*8 ,    dimension(n_saved_states_soc, 1, n_k_points),                   INTENT(IN) :: KS_eigen
  real*8,     dimension(n_saved_states_soc, 1 ,n_k_points),                   INTENT(IN) :: occ_numbers
  complex*16, dimension(n_basis_soc, n_saved_states_soc, 1, n_k_points_task), INTENT(IN) :: KS_vec_complex
  real*8,                                                                     INTENT(IN) :: chemical_potential
  real*8,     dimension(n_full_points), target,                               INTENT(IN) :: partition_tab
  integer,    dimension(n_species),                                           INTENT(IN) :: l_shell_max 
!  INPUTS
!    o KS_eigen
!    o occ_numbers
!    o KS_vec_complex
!    o chemical_potential
!    o partition_tab
!    o l_shell_max
!  OUTPUT
!    None
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  TODO
!    Unfork this from the SR version.
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

  character(*), parameter :: func = 'get_dipolematrix_soc'
  !  begin work

  write(info_str,'(6X,A,1X,I4)') "Dipole Matrix post processing starts"
  call localorb_info ( info_str )

  if (n_basis_soc .ne. 2*n_basis) then
    write(info_str,'(1X,A,A)') '* Incorrect number of SOC-perturbed basis elements specified.  Exiting.'
    call aims_stop(info_str, func)
  end if
  if (mod(n_basis_soc_coll,2) .eq. 1) then                                       
    call aims_stop('The number of collinear basis functions in SOC is odd, exiting.', func)
  end if                                                                                                                                                          
  if (n_basis_soc_ncoll .gt. 0) then                                            
    call aims_stop('This subroutine does not support usage of non-collinear basis functions&
                   & in SOC, exiting.', func)                                        
  end if                                                                         

  call convert_sr_to_soc_environment ()
  call get_state_minmax_k(KS_eigen, n_state_min_in, n_state_max_in)

  if(k_point_mom==0)then
      call open_hdf5('dipmat.h5', file_id, plist_id)
      call open_hdf5_kart("_Dipolematrix_", file_id, plist_id,n_state_min_in,&
                          n_state_max_in)
      call mpi_barrier(mpi_comm_global,info)
      call out_bands_kart(KS_eigen,file_id,plist_id,n_state_min_in,&
                          n_state_max_in)
      call out_occ_kart(occ_numbers,file_id,plist_id,n_state_min_in,&
                          n_state_max_in)
      call out_k_points_kart(file_id,plist_id)
      call outmetafile(KS_eigen,occ_numbers,chemical_potential,file_id,&
                       plist_id,n_state_min_in,n_state_max_in)
  endif


  write(info_str,'(2X,A,I7,I7)') 'n_state_min_in, n_state_max_in', n_state_min_in, n_state_max_in
  call localorb_info ( info_str )
  call allocate_mommat()
  call allocate_mommat_k()

  call allocate_moment_one_p1_SOC(n_state_min_in, n_state_max_in)
  call allocate_moment_two_p1_SOC(n_state_min_in, n_state_max_in)
  call allocate_moment_three_p1_SOC(n_state_min_in, n_state_max_in)

  call calculate_dipmat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up,&
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
        call calc_moment_p0_SOC(moment_one, &
                    mommat_full_w_complex_up, &
                    mommat_full_w_complex_low,&
                    KS_vec_complex(:,:,:,i_k), new_k_point, 1, &
                    n_state_min_in, n_state_max_in)
        if (k_point_mom==0) then
          call out_moment_kart(moment_one,KS_eigen,new_k_point,&
                               "_Dipolematrix_", file_id, plist_id,1,&
                               n_state_min_in,n_state_max_in)
        endif
      endif
    endif
  enddo
  call calculate_dipmat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up, &
                             mommat_full_oned_low,2)

  i_k = 0
  do new_k_point = 1,n_k_points
    if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points) then
      ! k-point stored on current mpi task
      i_k = i_k+1 !new_k_point
      if (new_k_point==k_point_mom.or.k_point_mom==0) then
        call construct_overlap( mommat_full_oned_up, mommat_full_w_up, &
                                 mommat_full_w_complex_up, new_k_point, &
                                 work_ovl_mom )
        call construct_overlap( mommat_full_oned_low, mommat_full_w_low, &
                                 mommat_full_w_complex_low, new_k_point, &
                                 work_ovl_mom )
        call calc_moment_p0_SOC(moment_two, &
                    mommat_full_w_complex_up, &
                    mommat_full_w_complex_low,&
                    KS_vec_complex(:,:,:,i_k), new_k_point, 2, &
                    n_state_min_in, n_state_max_in)
        if (k_point_mom==0) then
          call out_moment_kart(moment_two,KS_eigen,new_k_point,&
                                "_Dipolematrix_", file_id, plist_id,2, &
                                 n_state_min_in,n_state_max_in)
        endif
      endif
    endif
  enddo
  call calculate_dipmat_p0 ( partition_tab, l_shell_max, mommat_full_oned_up,&
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
        call calc_moment_p0_SOC(moment_three, &
                    mommat_full_w_complex_up, &
                    mommat_full_w_complex_low,&
                    KS_vec_complex(:,:,:,i_k), new_k_point, 3, &
                    n_state_min_in, n_state_max_in)
        if (k_point_mom==0)then 
          call out_moment_kart(moment_three,KS_eigen,new_k_point,&
                                "_Dipolematrix_", file_id, plist_id,3, &
                                n_state_min_in,n_state_max_in)
        endif
      endif
    endif
  enddo

  call revert_soc_to_sr_environment ()

  ! Here all necessary ingredients for the calculation of the localized
  ! Boys orbitals have been calculated. Now it is just optimizing the
  ! functional from the external routine and outputting it. Boys orbitals
  ! are for non periodic systems only (k_point_mom=1). Implementation of
  ! periodic Wannier orbitals shouldn't be so hard... Also, all occupied
  ! orbitals should be taken into account! How to always consider all of them?
  if (flag_out_boys) then
     if (k_point_mom/=0.and.myid==modulo(k_point_mom, n_tasks))then
        write(info_str,'(6X,A,1X,I4)') 'Boys centers and rotation matrix calculation'
        call localorb_info ( info_str )
        call boys_centers(moment_one, moment_two, moment_three, n_state_max_in) 
     endif
  else ! just output the momentum matrix
     if (k_point_mom/=0.and.myid==modulo(k_point_mom, n_tasks))then 
         call out_moment(moment_one,moment_two,moment_three, KS_eigen, &
                     k_point_mom,n_state_min_in,n_state_max_in)
     endif
  endif

  call clean_mommat()
  call clean_mommat_final()
  if(k_point_mom==0)then
    call mpi_barrier(mpi_comm_global,info)
    call close_hdf5(file_id,plist_id)
  endif
  write(info_str,'(6X,A,1X,I4)') "Dipole Matrix post processing finished"
  call localorb_info ( info_str )
 
end subroutine get_dipolematrix_soc
!******	
