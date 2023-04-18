!****h* FHI-aims/ks_wrapper
!  NAME
!    ks_wrapper
!  AUTHOR
!    Victor Yu, Duke University, FHI-aims team
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is
!    subject to the terms and conditions of the respective license
!    agreement.
!  SYNOPSIS
module ks_wrapper

   implicit none

   private

   public :: solve_KS_elsi_serial
   public :: solve_KS_elsi_parallel
   public :: solve_KS_elsi_dm_parallel

   interface solve_KS_elsi_serial
      module procedure solve_KS_elsi_serial_real
      module procedure solve_KS_elsi_serial_cmplx
   end interface

   interface solve_KS_elsi_parallel
      module procedure solve_KS_elsi_parallel_real
      module procedure solve_KS_elsi_parallel_cmplx
   end interface

contains
!-------------------------------------------------------------------------------
!****s* FHI-aims/!****s* FHI-aims/solve_KS_elsi_serial_real
!  NAME
!    solve_KS_elsi_serial_real
!  SYNOPSIS
subroutine solve_KS_elsi_serial_real(ham_s,ovlp_s,eval,evec_s,i_spin,i_kpt)
!  AUTHOR
!    Victor Yu, Duke University, FHI-aims team
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is
!    subject to the terms and conditions of the respective license
!    agreement.
!  USES
   use aims_memory_tracking, only: aims_allocate,aims_deallocate
   use dimensions, only: n_states_k,n_states,n_basis,n_core_states
   use elsi_wrapper, only: eh_scf,rwh_w,aims_elsi_ev,aims_elsi_write_mat,&
       aims_elsi_set_illcond_check,aims_elsi_get_n_illcond,&
       aims_elsi_set_output
   use runtime_choices, only: out_h_elsi,out_s_elsi,flag_KS_k_points,&
       frozen_core_scf
   use separate_core_states, only: permute_ham_lapack,permute_ovlp_lapack,&
       permute_evec_lapack,fetch_core_eval_lapack

   implicit none

   real*8, intent(in) :: ham_s(n_basis*(n_basis+1)/2)
   real*8, intent(in) :: ovlp_s(n_basis*(n_basis+1)/2)
   real*8, intent(out) :: eval(n_states)
   real*8, intent(out) :: evec_s(n_basis,n_states)
   integer, intent(in) :: i_spin
   integer, intent(in) :: i_kpt

   real*8, dimension(:), allocatable :: eval_tmp
   real*8, dimension(:,:), allocatable :: ham_work
   real*8, dimension(:,:), allocatable :: ovlp_work
   real*8, dimension(:,:), allocatable :: evec_work

   integer :: i_row
   integer :: i_col
   integer :: i_count
   integer :: n_kpt

   character(50) :: mat_file

   call aims_allocate(ovlp_work,n_basis,n_basis,"ovlp_work")
   call aims_allocate(ham_work,n_basis,n_basis,"ham_work")

   ovlp_work = 0.d0
   ham_work = 0.d0

   ! Unpack packed matrices
   i_count = 0

   do i_col = 1,n_basis
      do i_row = 1,i_col
         i_count = i_count+1

         ovlp_work(i_row,i_col) = ovlp_s(i_count)
         ovlp_work(i_col,i_row) = ovlp_s(i_count)
         ham_work(i_row,i_col) = ham_s(i_count)
         ham_work(i_col,i_row) = ham_s(i_count)
      end do
   end do

   call aims_allocate(evec_work,n_basis,n_basis,"evec_work")
   call aims_allocate(eval_tmp,n_basis,"eval_tmp")

   ! Skip singularity check when possible
   if(flag_KS_k_points(i_kpt) == 0) then
      call aims_elsi_set_illcond_check(eh_scf,0)
   else
      call aims_elsi_set_illcond_check(eh_scf,1)
   end if

   if(frozen_core_scf) then
      call permute_ovlp_lapack(ovlp_work,evec_work)
      call permute_ham_lapack(ham_work,ovlp_work,evec_work)
   end if

   ! Output matrices
   if(out_h_elsi) then
      write(mat_file,"(A,I2.2,A,I6.6,A)") "H_spin_",i_spin,"_kpt_",i_kpt,".csc"
      call aims_elsi_write_mat(rwh_w,trim(mat_file),ham_work)
   end if

   if(out_s_elsi .and. i_spin == 1) then
      write(mat_file,"(A,I6.6,A)") "S_spin_01_kpt_",i_kpt,".csc"
      call aims_elsi_write_mat(rwh_w,trim(mat_file),ovlp_work)
   end if

   if(frozen_core_scf) then
      call fetch_core_eval_lapack(ham_work,ovlp_work,eval_tmp(1:n_core_states))
      call aims_elsi_ev(eh_scf,&
           ham_work(n_core_states+1:n_basis,n_core_states+1:n_basis),&
           ovlp_work(n_core_states+1:n_basis,n_core_states+1:n_basis),&
           eval_tmp(n_core_states+1:n_basis),&
           evec_work(n_core_states+1:n_basis,n_core_states+1:n_basis))
   else
      call aims_elsi_ev(eh_scf,ham_work,ovlp_work,eval_tmp,evec_work)
   end if

   if(frozen_core_scf) then
      call permute_evec_lapack(ovlp_work,evec_work)
   end if

   eval = eval_tmp(1:n_states)
   evec_s = evec_work(:,1:n_states)

   call aims_deallocate(ovlp_work,"ovlp_work")
   call aims_deallocate(ham_work,"ham_work")
   call aims_deallocate(evec_work,"evec_work")
   call aims_deallocate(eval_tmp,"eval_tmp")

   call aims_elsi_get_n_illcond(eh_scf,n_states_k(i_kpt))
   n_states_k(i_kpt) = n_basis-n_states_k(i_kpt)

   ! Overlap singular for this k-point?
   call aims_elsi_get_n_illcond(eh_scf,flag_KS_k_points(i_kpt))

   ! Avoid abundant output
   call aims_elsi_set_output(eh_scf,0)

end subroutine
!******
!-------------------------------------------------------------------------------
!****s* FHI-aims/!****s* FHI-aims/solve_KS_elsi_serial_cmplx
!  NAME
!    solve_KS_elsi_serial_cmplx
!  SYNOPSIS
subroutine solve_KS_elsi_serial_cmplx(ham_s,ovlp_s,eval,evec_s,i_spin,i_kpt)
!  AUTHOR
!    Victor Yu, Duke University, FHI-aims team
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is
!    subject to the terms and conditions of the respective license
!    agreement.
!  USES
   use aims_memory_tracking, only: aims_allocate,aims_deallocate
   use dimensions, only: n_states_k,n_states,n_basis,n_core_states
   use elsi_wrapper, only: eh_scf,rwh_w,aims_elsi_ev,aims_elsi_write_mat,&
       aims_elsi_set_illcond_check,aims_elsi_get_n_illcond,&
       aims_elsi_set_output
   use runtime_choices, only: out_h_elsi,out_s_elsi,flag_KS_k_points,&
       frozen_core_scf
   use separate_core_states, only: permute_ham_lapack,permute_ovlp_lapack,&
       permute_evec_lapack,fetch_core_eval_lapack

   implicit none

   complex*16, intent(in) :: ham_s(n_basis*(n_basis+1)/2)
   complex*16, intent(in) :: ovlp_s(n_basis*(n_basis+1)/2)
   real*8, intent(out) :: eval(n_states)
   complex*16, intent(out) :: evec_s(n_basis,n_states)
   integer, intent(in) :: i_spin
   integer, intent(in) :: i_kpt

   real*8, dimension(:), allocatable :: eval_tmp
   complex*16, dimension(:,:), allocatable :: ham_work
   complex*16, dimension(:,:), allocatable :: ovlp_work
   complex*16, dimension(:,:), allocatable :: evec_work

   integer :: i_row
   integer :: i_col
   integer :: i_count
   integer :: n_kpt

   character(50) :: mat_file

   call aims_allocate(ovlp_work,n_basis,n_basis,"ovlp_work")
   call aims_allocate(ham_work,n_basis,n_basis,"ham_work")

   ovlp_work = (0.d0,0.d0)
   ham_work = (0.d0,0.d0)

   ! Unpack packed matrices
   i_count = 0

   do i_col = 1,n_basis
      do i_row = 1,i_col
         i_count = i_count+1

         ovlp_work(i_row,i_col) = ovlp_s(i_count)
         ovlp_work(i_col,i_row) = conjg(ovlp_s(i_count))
         ham_work(i_row,i_col) = ham_s(i_count)
         ham_work(i_col,i_row) = conjg(ham_s(i_count))
      end do
   end do

   call aims_allocate(evec_work,n_basis,n_basis,"evec_work")
   call aims_allocate(eval_tmp,n_basis,"eval_tmp")

   ! Skip singularity check when possible
   if(flag_KS_k_points(i_kpt) == 0) then
      call aims_elsi_set_illcond_check(eh_scf,0)
   else
      call aims_elsi_set_illcond_check(eh_scf,1)
   end if

   if(frozen_core_scf) then
      call permute_ovlp_lapack(ovlp_work,evec_work)
      call permute_ham_lapack(ham_work,ovlp_work,evec_work)
   end if

   ! Output matrices
   if(out_h_elsi) then
      write(mat_file,"(A,I2.2,A,I6.6,A)") "H_spin_",i_spin,"_kpt_",i_kpt,".csc"
      call aims_elsi_write_mat(rwh_w,trim(mat_file),ham_work)
   end if

   if(out_s_elsi .and. i_spin == 1) then
      write(mat_file,"(A,I6.6,A)") "S_spin_01_kpt_",i_kpt,".csc"
      call aims_elsi_write_mat(rwh_w,trim(mat_file),ovlp_work)
   end if

   if(frozen_core_scf) then
      call fetch_core_eval_lapack(ham_work,ovlp_work,eval_tmp(1:n_core_states))
      call aims_elsi_ev(eh_scf,&
           ham_work(n_core_states+1:n_basis,n_core_states+1:n_basis),&
           ovlp_work(n_core_states+1:n_basis,n_core_states+1:n_basis),&
           eval_tmp(n_core_states+1:n_basis),&
           evec_work(n_core_states+1:n_basis,n_core_states+1:n_basis))
   else
      call aims_elsi_ev(eh_scf,ham_work,ovlp_work,eval_tmp,evec_work)
   end if

   if(frozen_core_scf) then
      call permute_evec_lapack(ovlp_work,evec_work)
   end if

   eval = eval_tmp(1:n_states)
   evec_s = evec_work(:,1:n_states)

   call aims_deallocate(ovlp_work,"ovlp_work")
   call aims_deallocate(ham_work,"ham_work")
   call aims_deallocate(evec_work,"evec_work")
   call aims_deallocate(eval_tmp,"eval_tmp")

   call aims_elsi_get_n_illcond(eh_scf,n_states_k(i_kpt))
   n_states_k(i_kpt) = n_basis-n_states_k(i_kpt)

   ! Overlap singular for this k-point?
   call aims_elsi_get_n_illcond(eh_scf,flag_KS_k_points(i_kpt))

   ! Avoid abundant output
   call aims_elsi_set_output(eh_scf,0)

end subroutine
!******
!-------------------------------------------------------------------------------
!****s* FHI-aims/!****s* FHI-aims/solve_KS_elsi_parallel_real
!  NAME
!    solve_KS_elsi_parallel_real
!  SYNOPSIS
subroutine solve_KS_elsi_parallel_real(eval,evec_s,i_spin)
!  AUTHOR
!    Victor Yu, Duke University, FHI-aims team
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is
!    subject to the terms and conditions of the respective license
!    agreement.
!  USES

   use aims_memory_tracking, only: aims_allocate,aims_deallocate
   use dimensions, only: n_states_k,n_states,n_basis,n_core_states,n_periodic
   use elsi_wrapper, only: eh_scf,rwh_w,aims_elsi_ev,aims_elsi_write_mat,&
       aims_elsi_get_n_illcond
   use runtime_choices, only: packed_matrix_format,PM_index,elsi_solver,&
       frozen_core_scf,collect_eigenvectors,out_h_elsi,out_s_elsi
   use scalapack_wrapper, only: ham,ovlp,eigenvec,my_k_point,full_ovlp_ready,&
       my_scalapack_id,set_full_matrix_real,collect_eigenvectors_scalapack
   use separate_core_states, only: permute_ham_scalapack,&
       permute_ovlp_scalapack,permute_evec_scalapack,fetch_core_eval_scalapack,&
       ham_real_vv,ovlp_real_vv,evec_real_vv

   implicit none

   real*8, intent(inout) :: eval(n_states)
   real*8, intent(out) :: evec_s(n_basis,n_states)
   integer, intent(in) :: i_spin

   real*8, dimension(:), allocatable :: eval_tmp

   character(50) :: mat_file

   call set_full_matrix_real(ham(:,:,i_spin))

   if(.not. full_ovlp_ready) then
      call set_full_matrix_real(ovlp)

      full_ovlp_ready = .true.

      if(frozen_core_scf) then
         call permute_ovlp_scalapack()
      end if
   end if

   if(frozen_core_scf) then
      call permute_ham_scalapack(i_spin)
   end if

   call aims_allocate(eval_tmp,n_basis,"eval_tmp")

   if(elsi_solver == 5) then
      eval_tmp(1:n_states) = eval
   end if

   ! Output matrices
   if(out_h_elsi) then
      write(mat_file,"(A,I2.2,A,I6.6,A)") "H_spin_",i_spin,"_kpt_",my_k_point,&
         ".csc"
      call aims_elsi_write_mat(rwh_w,trim(mat_file),ham(:,:,i_spin))
   end if

   if(out_s_elsi .and. i_spin == 1) then
      write(mat_file,"(A,I6.6,A)") "S_spin_01_kpt_",my_k_point,".csc"
      call aims_elsi_write_mat(rwh_w,trim(mat_file),ovlp)
   end if

   if(frozen_core_scf) then
      call fetch_core_eval_scalapack(i_spin,eval_tmp(1:n_core_states))
      call aims_elsi_ev(eh_scf,ham_real_vv,ovlp_real_vv,&
           eval_tmp(n_core_states+1:n_basis),evec_real_vv)
   else
      call aims_elsi_ev(eh_scf,ham(:,:,i_spin),ovlp,eval_tmp,&
           eigenvec(:,:,i_spin))
   end if

   eval = eval_tmp(1:n_states)

   call aims_deallocate(eval_tmp,"eval_tmp")

   if(frozen_core_scf) then
      call permute_evec_scalapack(i_spin)
   end if

   if(collect_eigenvectors) then
      call collect_eigenvectors_scalapack(evec_s,i_spin)
   end if

   call aims_elsi_get_n_illcond(eh_scf,n_states_k(my_k_point))
   n_states_k(my_k_point) = n_basis-n_states_k(my_k_point)

   if(n_periodic > 0 .or. packed_matrix_format == PM_index) then
      if(my_scalapack_id /= 0) then
         eval = 0.d0
         n_states_k = 0
      end if
   end if

end subroutine
!******
!-------------------------------------------------------------------------------
!****s* FHI-aims/!****s* FHI-aims/solve_KS_elsi_parallel_cmplx
!  NAME
!    solve_KS_elsi_parallel_cmplx
!  SYNOPSIS
subroutine solve_KS_elsi_parallel_cmplx(eval,evec_s,i_spin)
!  AUTHOR
!    Victor Yu, Duke University, FHI-aims team
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is
!    subject to the terms and conditions of the respective license
!    agreement.
!  USES

   use aims_memory_tracking, only: aims_allocate,aims_deallocate
   use dimensions, only: n_states_k,n_states,n_basis,n_core_states,n_periodic
   use elsi_wrapper, only: eh_scf,rwh_w,aims_elsi_ev,aims_elsi_write_mat,&
       aims_elsi_get_n_illcond
   use runtime_choices, only: packed_matrix_format,PM_index,elsi_solver,&
       frozen_core_scf,collect_eigenvectors,out_h_elsi,out_s_elsi
   use scalapack_wrapper, only: ham_complex,ovlp_complex,eigenvec_complex,&
       my_k_point,full_ovlp_ready,my_scalapack_id,set_full_matrix_complex,&
       collect_eigenvectors_scalapack_complex
   use separate_core_states, only: permute_ham_scalapack,&
       permute_ovlp_scalapack,permute_evec_scalapack,fetch_core_eval_scalapack,&
       ham_cmplx_vv,ovlp_cmplx_vv,evec_cmplx_vv

   implicit none

   real*8, intent(inout) :: eval(n_states)
   complex*16, intent(out) :: evec_s(n_basis,n_states)
   integer, intent(in) :: i_spin

   real*8, dimension(:), allocatable :: eval_tmp

   character(50) :: mat_file

   call set_full_matrix_complex(ham_complex(:,:,i_spin))

   if(.not. full_ovlp_ready) then
      call set_full_matrix_complex(ovlp_complex)

      full_ovlp_ready = .true.

      if(frozen_core_scf) then
         call permute_ovlp_scalapack()
      end if
   end if

   if(frozen_core_scf) then
      call permute_ham_scalapack(i_spin)
   end if

   call aims_allocate(eval_tmp,n_basis,"eval_tmp")

   if(elsi_solver == 5) then
      eval_tmp(1:n_states) = eval
   end if

   ! Output matrices
   if(out_h_elsi) then
      write(mat_file,"(A,I2.2,A,I6.6,A)") "H_spin_",i_spin,"_kpt_",my_k_point,&
         ".csc"
      call aims_elsi_write_mat(rwh_w,trim(mat_file),ham_complex(:,:,i_spin))
   end if

   if(out_s_elsi .and. i_spin == 1) then
      write(mat_file,"(A,I6.6,A)") "S_spin_01_kpt_",my_k_point,".csc"
      call aims_elsi_write_mat(rwh_w,trim(mat_file),ovlp_complex)
   end if

   if(frozen_core_scf) then
      call fetch_core_eval_scalapack(i_spin,eval_tmp(1:n_core_states))
      call aims_elsi_ev(eh_scf,ham_cmplx_vv,ovlp_cmplx_vv,&
           eval_tmp(n_core_states+1:n_basis),evec_cmplx_vv)
   else
      call aims_elsi_ev(eh_scf,ham_complex(:,:,i_spin),ovlp_complex,eval_tmp,&
           eigenvec_complex(:,:,i_spin))
   end if

   eval = eval_tmp(1:n_states)

   call aims_deallocate(eval_tmp,"eval_tmp")

   if(frozen_core_scf) then
      call permute_evec_scalapack(i_spin)
   end if

   if(collect_eigenvectors) then
      call collect_eigenvectors_scalapack_complex(evec_s,i_spin)
   end if

   call aims_elsi_get_n_illcond(eh_scf,n_states_k(my_k_point))
   n_states_k(my_k_point) = n_basis-n_states_k(my_k_point)

   if(n_periodic > 0 .or. packed_matrix_format == PM_index) then
      if(my_scalapack_id /= 0) then
         eval = 0.d0
         n_states_k = 0
      end if
   end if

end subroutine
!******
!-------------------------------------------------------------------------------
!****s* FHI-aims/!****s* FHI-aims/solve_KS_elsi_dm_parallel
!  NAME
!    solve_KS_elsi_dm_parallel
!  SYNOPSIS
subroutine solve_KS_elsi_dm_parallel(ebs,i_spin)
!  AUTHOR
!    Victor Yu, Duke University, FHI-aims team
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is
!    subject to the terms and conditions of the respective license
!    agreement.
!  USES

   use elsi_wrapper, only: eh_scf,rwh_w,pexsi_mu_min,pexsi_mu_max,aims_elsi_dm,&
       aims_elsi_write_mat,aims_elsi_set_pexsi_mu_max,&
       aims_elsi_set_pexsi_mu_min,aims_elsi_get_pexsi_mu_max,&
       aims_elsi_get_pexsi_mu_min
   use pbc_lists, only: k_weights
   use physics, only: dv_hartree_min,dv_hartree_max
   use runtime_choices, only: real_eigenvectors,first_elsi_call,elsi_solver,&
       out_h_elsi,out_s_elsi
   use scalapack_wrapper, only: ham,ovlp,eigenvec,ham_complex,ovlp_complex,&
       eigenvec_complex,my_k_point,full_ovlp_ready,set_full_matrix_real,&
       set_full_matrix_complex

   implicit none

   real*8, intent(out) :: ebs
   integer, intent(in) :: i_spin

   character(50) :: mat_file

   if(real_eigenvectors) then
      call set_full_matrix_real(ham(:,:,i_spin))
   else
      call set_full_matrix_complex(ham_complex(:,:,i_spin))
   end if

   if(.not. full_ovlp_ready) then
      if(real_eigenvectors) then
         call set_full_matrix_real(ovlp)
      else
         call set_full_matrix_complex(ovlp_complex)
      end if

      full_ovlp_ready = .true.
   end if

   if(elsi_solver == 3) then
      if(first_elsi_call) then
         pexsi_mu_min = -2.d0
         pexsi_mu_max = 1.d0
      else
         pexsi_mu_min = max(pexsi_mu_min+dv_hartree_min,-2.d0)
         pexsi_mu_max = min(pexsi_mu_max+dv_hartree_max,1.d0)
      end if

      call aims_elsi_set_pexsi_mu_min(eh_scf,pexsi_mu_min)
      call aims_elsi_set_pexsi_mu_max(eh_scf,pexsi_mu_max)
   end if

   ! Output matrices
   if(out_h_elsi) then
      write(mat_file,"(A,I2.2,A,I6.6,A)") "H_spin_",i_spin,"_kpt_",my_k_point,&
         ".csc"

      if(real_eigenvectors) then
         call aims_elsi_write_mat(rwh_w,trim(mat_file),ham(:,:,i_spin))
      else
         call aims_elsi_write_mat(rwh_w,trim(mat_file),ham_complex(:,:,i_spin))
      end if
   end if

   if(out_s_elsi .and. i_spin == 1) then
      write(mat_file,"(A,I6.6,A)") "S_spin_01_kpt_",my_k_point,".csc"

      if(real_eigenvectors) then
         call aims_elsi_write_mat(rwh_w,trim(mat_file),ovlp)
      else
         call aims_elsi_write_mat(rwh_w,trim(mat_file),ovlp_complex)
      end if
   end if

   if(real_eigenvectors) then
      ! Use eigenvec to store density matrix
      call aims_elsi_dm(eh_scf,ham(:,:,i_spin),ovlp,eigenvec(:,:,i_spin),ebs)

      ! Apply k weights
      ham(:,:,i_spin) = k_weights(my_k_point)*eigenvec(:,:,i_spin)
   else
      ! Use eigenvec_complex to store density matrix
      call aims_elsi_dm(eh_scf,ham_complex(:,:,i_spin),ovlp_complex,&
           eigenvec_complex(:,:,i_spin),ebs)

      ! Apply k weights
      ham_complex(:,:,i_spin) = k_weights(my_k_point)&
         *eigenvec_complex(:,:,i_spin)
   end if

   if(elsi_solver == 3) then
      call aims_elsi_get_pexsi_mu_min(eh_scf,pexsi_mu_min)
      call aims_elsi_get_pexsi_mu_max(eh_scf,pexsi_mu_max)
   end if

end subroutine
!******
end module
!******
