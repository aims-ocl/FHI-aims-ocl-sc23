!****s* FHI-aims/!****s* FHI-aims/solve_KS_eigen
!  NAME
!    solve_KS_eigen
!  SYNOPSIS
subroutine solve_KS_eigen(overlap_matrix,hamiltonian,KS_eigenvalue,&
   KS_eigenvector,KS_eigenvector_complex)
!  PURPOSE
!    When use_elsi is .true., this subroutine computes the Kohn-Sham orbitals or
!    the density matrix by ELSI, using one of the solvers supported in ELSI:
!    ELPA, LAPACK, libOMM, NTPoly, PEXSI, SLEPc-SIPs.
!    When use_elsi is .false., this subroutine computes the Kohn-Sham orbitals
!    using ELPA (2013), ScaLAPACK, or LAPACK.
!  AUTHOR
!    Victor Yu, Duke University, FHI-aims team
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is
!    subject to the terms and conditions of the respective license
!    agreement.
!  USES
   use aims_memory_tracking, only: aims_allocate,aims_deallocate
   use applicable_citations, only: cite_reference
   use dimensions, only: use_constraint,n_states_k,n_periodic,use_lc_wpbeh,&
       use_periodic_hf,n_centers_basis_I,use_hartree_fock,use_periodic_hf,&
       n_hamiltonian_matrix_size,n_states,n_k_points,n_basis,n_k_points_task,&
       n_spin,n_core_states
   use elsi_wrapper, only: eh_scf,aims_elsi_set_output,&
       aims_elsi_set_sips_ev_max,aims_elsi_set_sips_ev_min,&
       aims_elsi_get_ovlp_ev_max,aims_elsi_get_ovlp_ev_min
   use hartree_fock_p0, only: hf_exchange_matr_real,hf_exchange_matr_real_SR,&
       hf_exchange_matr_complex,hf_exchange_matr_complex_SR
   use ks_wrapper, only: solve_KS_elsi_serial,solve_KS_elsi_parallel,&
       solve_KS_elsi_dm_parallel
   use localorb_io, only: OL_NORM,localorb_info,use_unit
   use mpi_tasks, only: myid,n_tasks,aims_stop
   use physics, only: ev_sum
   use rel_x2c_mod, only: dim_matrix_rel
   use runtime_choices, only: use_elsi,use_elsi_dm,use_elsi_ev,elsi_solver,&
       elsi_extrap_dm,elsi_write_dm,first_elsi_call,out_matrices_elsi,&
       out_h_elsi,out_s_elsi,packed_matrix_format,PM_index,use_load_balancing,&
       use_scalapack,real_eigenvectors,use_density_matrix,flag_rel,REL_x2c,&
       REL_4c_dks,hybrid_coeff,lc_dielectric_constant,elsi_out_level,&
       flag_KS_k_points
   use scalapack_wrapper, only: ham,ham_complex,my_k_point,l_row,l_col,&
       setup_hamiltonian_scalapack,solve_evp_scalapack,&
       solve_evp_scalapack_complex
   use synchronize_mpi, only: sync_eigenvalues,sync_vector_integer,sync_vector
   use synchronize_mpi_basic, only: get_min_double

   implicit none

   real*8, intent(in) :: hamiltonian(n_hamiltonian_matrix_size,n_spin)
   real*8, intent(in) :: overlap_matrix(n_hamiltonian_matrix_size)
   real*8, intent(inout) :: KS_eigenvalue(n_states,n_spin,n_k_points)
   real*8, intent(inout) :: KS_eigenvector(n_basis,n_states,n_spin,&
      n_k_points_task)
   complex*16, intent(inout) :: KS_eigenvector_complex(n_basis,n_states,n_spin,&
      n_k_points_task)

   real*8, dimension(:), allocatable :: eval_tmp
   real*8, dimension(:,:,:), allocatable :: work_ham ! For PM_none
   real*8, dimension(:,:), allocatable :: work_ovlp ! For PM_none
   real*8, dimension(:,:), allocatable :: ham_work
   real*8, dimension(:), allocatable :: ovlp_work
   real*8, dimension(:,:), allocatable :: ham_work_s
   real*8, dimension(:,:), allocatable :: ovlp_work_s
   real*8, dimension(:,:), allocatable :: evec_work_s
   real*8, dimension(:), allocatable :: ovlp_ev_min_k
   real*8, dimension(:), allocatable :: ovlp_ev_max_k

   complex*16, dimension(:,:), allocatable :: ham_work_c
   complex*16, dimension(:), allocatable :: ovlp_work_c
   complex*16, dimension(:,:), allocatable :: ham_work_cs
   complex*16, dimension(:,:), allocatable :: ovlp_work_cs
   complex*16, dimension(:,:), allocatable :: evec_work_cs

   integer :: i
   integer :: j
   integer :: i_spin
   integer :: i_kpt
   integer :: i_count
   integer :: n_kpt
   integer :: ierr
   integer :: n_states_new

   real*8 :: ev_shift
   real*8 :: ham_diag_min

   logical :: sing_output
   logical :: lapack_output

   character :: aux
   character(50) :: mat_file
   character(200) :: msg

   if(.not. use_elsi) then
      KS_eigenvalue = 0.d0
      sing_output = .false.
   else
      call cite_reference("ELSI")

      ! Use SLEPc-SIPs?
      if(elsi_solver /= 5) then
         KS_eigenvalue = 0.d0
      else
         ham_diag_min = 0.d0

         do i = 1,n_basis
            if(l_col(i) > 0 .and. l_row(i) > 0) then
               ham_diag_min = min(ham_diag_min,ham(l_row(i),l_col(i),1))
            end if
         end do

         call get_min_double(ev_shift,ham_diag_min)

         call aims_elsi_set_sips_ev_min(eh_scf,ev_shift-0.2d0)
         call aims_elsi_set_sips_ev_max(eh_scf,2.d0)

         ev_shift = ev_shift-KS_eigenvalue(1,1,1)
         KS_eigenvalue = KS_eigenvalue+ev_shift
      end if

      ! Ill-conditioning output
      if(first_elsi_call .and. use_elsi_ev .and. elsi_solver == 1) then
         sing_output = .true.
      else
         sing_output = .false.
      end if

      write(aux,"(I1)") int(log(real(n_tasks,8))/log(10.d0),4)+1

      if(sing_output) then
         call aims_allocate(ovlp_ev_min_k,n_k_points,"ovlp_ev_min_k")
         call aims_allocate(ovlp_ev_max_k,n_k_points,"ovlp_ev_max_k")

         ovlp_ev_min_k = 0.d0
         ovlp_ev_max_k = 0.d0

         if(use_elsi_ev .and. elsi_solver == 1 .and. n_periodic > 0) then
            if(use_scalapack) then
               write(msg,"(2X,A)") "Singularity check in k-point 1 (analysis"//&
                  " for other k-points may follow below):"
               call localorb_info(msg,use_unit,"(A)",OL_norm)
            else if(n_k_points >= n_tasks) then
               write(msg,"(2X,A,I"//aux//",A)") "Singularity check in"//&
                  " k-point ",n_tasks,", task 0 (analysis for other"//&
                  " k-points/tasks may follow below):"
               call localorb_info(msg,use_unit,"(A)",OL_norm)
            end if
         end if
      end if

      ! Corner case
      if(.not. first_elsi_call .and. .not. use_scalapack .and. use_elsi_ev&
         .and. elsi_solver == 1 .and. n_periodic > 0&
         .and. n_k_points >= n_tasks) then
         if(flag_KS_k_points(n_tasks) /= 0) then
            write(msg,"(2X,A,I"//aux//",A)") "Singularity check in k-point ",&
               n_tasks,", task 0:"
            call localorb_info(msg,use_unit,"(A)",OL_norm)
         end if
      end if

      ! ELSI matrices output
      if(first_elsi_call) then
         if((elsi_extrap_dm .and. elsi_write_dm) .or. out_matrices_elsi) then
            out_s_elsi = .true.
         end if

         if(out_matrices_elsi) then
            out_h_elsi = .true.
         end if
      else
         out_h_elsi = .false.
         out_s_elsi = .false.
      end if
   end if ! Use ELSI?

   n_states_k = 0

   ! Main work
   if(n_periodic > 0 .or. packed_matrix_format == PM_index&
      .or. flag_rel == REL_x2c .or. flag_rel == REL_4c_dks) then
      if(real_eigenvectors) then
         if(use_scalapack) then
            if(use_periodic_hf .and. allocated(hf_exchange_matr_real)) then
               ! FIXME This is ugly and should not be in this subroutine
               ! Add exchange matrix to Hamiltonian
               if(use_lc_wpbeh) then
                  if(hybrid_coeff /= 0.d0) then
                     ham(:,:,:) = ham(:,:,:)-((1/lc_dielectric_constant)&
                        *hf_exchange_matr_real(:,:,1,:)+hybrid_coeff&
                        *hf_exchange_matr_real_SR(:,:,1,:))
                  else
                     ham(:,:,:) = ham(:,:,:)-(1/lc_dielectric_constant)&
                        *hf_exchange_matr_real(:,:,1,:)
                  end if
               else
                  ham(:,:,:) = ham(:,:,:)-hybrid_coeff&
                     *hf_exchange_matr_real(:,:,1,:)
               end if
            end if

            do i_spin = 1,n_spin
               if(.not. use_elsi) then
                  call solve_evp_scalapack(KS_eigenvalue(:,i_spin,my_k_point),&
                       KS_eigenvector(:,:,i_spin,1),i_spin)
               else if(use_elsi_ev) then
                  call solve_KS_elsi_parallel(&
                       KS_eigenvalue(:,i_spin,my_k_point),&
                       KS_eigenvector(:,:,i_spin,1),i_spin)

                  if(sing_output) then
                     call aims_elsi_get_ovlp_ev_min(eh_scf,&
                          ovlp_ev_min_k(my_k_point))
                     call aims_elsi_get_ovlp_ev_max(eh_scf,&
                          ovlp_ev_max_k(my_k_point))
                  end if
               else if(use_elsi_dm) then
                  call solve_KS_elsi_dm_parallel(ev_sum,i_spin)
               else
                  write(msg,"(2X,A)") '*** ERROR: Invalid "elsi_method".'
                  call localorb_info(msg,use_unit,"(A)",OL_norm)
                  call aims_stop
               end if
            end do
         else ! Use ScaLAPACK?
            n_kpt = 0

            do i_kpt = 1,n_k_points
               ! TODO: The use_load_balancing implementation here is blocking
               ! across MPI tasks and, worse yet, ugly. Needs to be cleaned up.
               if((myid == mod(i_kpt,n_tasks) .and. myid <= n_k_points)&
                  .or. use_load_balancing) then
                  n_kpt = n_kpt+1

                  if(.not. allocated(ham_work)) then
                     call aims_allocate(ham_work,n_basis*(n_basis+1)/2,n_spin,&
                          "ham_work")
                     call aims_allocate(ovlp_work,n_basis*(n_basis+1)/2,&
                          "ovlp_work")
                     call aims_allocate(ham_work_c,1,n_spin,"ham_work_c")
                     call aims_allocate(ovlp_work_c,1,"ovlp_work_c")
                  end if

                  call construct_hamiltonian_and_ovl(hamiltonian,&
                       overlap_matrix,ham_work,ovlp_work,ham_work_c,&
                       ovlp_work_c,i_kpt)

                  if(.not. (myid == mod(i_kpt,n_tasks)&
                     .and. myid <= n_k_points)) then
                     n_kpt = n_kpt-1
                     cycle
                  end if

                  ! FIXME This is ugly and should not be in this subroutine
                  ! Add the exchange matrix here if needed
                  if(use_periodic_hf&
                     .and. allocated(hf_exchange_matr_real)) then
                     if(use_lc_wpbeh) then
                        if(hybrid_coeff /= 0.d0) then
                           call get_hf_hamiltonian_real_p0(&
                                (1/lc_dielectric_constant)&
                                *hf_exchange_matr_real+hybrid_coeff&
                                *hf_exchange_matr_real_SR,ham_work,n_kpt)
                        else
                           call get_hf_hamiltonian_real_p0(&
                                (1/lc_dielectric_constant)&
                                *hf_exchange_matr_real,ham_work,n_kpt)
                        end if
                     else
                        call get_hf_hamiltonian_real_p0(hf_exchange_matr_real,&
                             ham_work,n_kpt)
                     end if
                  end if

                  ! Solve KS-equations for both spins
                  do i_spin = 1,n_spin
                     if(.not. use_elsi) then
                        if(n_kpt == 1 .and. .not. use_constraint) then
                           lapack_output = .true.
                        else
                           lapack_output = .false.
                        end if

                        call improve_real_eigenfunctions(ovlp_work,&
                             ham_work(:,i_spin),lapack_output,&
                             KS_eigenvalue(:,i_spin,i_kpt),&
                             KS_eigenvector(:,:,i_spin,n_kpt),i_kpt)
                     else
                        call solve_KS_elsi_serial(ham_work(:,i_spin),ovlp_work,&
                             KS_eigenvalue(:,i_spin,i_kpt),&
                             KS_eigenvector(:,:,i_spin,n_kpt),i_spin,i_kpt)

                        if(sing_output) then
                           call aims_elsi_get_ovlp_ev_min(eh_scf,&
                                ovlp_ev_min_k(i_kpt))
                           call aims_elsi_get_ovlp_ev_max(eh_scf,&
                                ovlp_ev_max_k(i_kpt))
                        end if
                     end if
                  end do
               end if
            end do

            if(allocated(ham_work)) then
               call aims_deallocate(ham_work,"ham_work")
            end if
            if(allocated(ovlp_work)) then
               call aims_deallocate(ovlp_work,"ovlp_work")
            end if
            if(allocated(ham_work_c)) then
               call aims_deallocate(ham_work_c,"ham_work_c")
            end if
            if(allocated(ovlp_work_c)) then
               call aims_deallocate(ovlp_work_c,"ovlp_work_c")
            end if
         end if ! Use ScaLAPACK?
      else ! Real eigenvectors?
         if(use_scalapack) then
            if(use_periodic_hf .and. allocated(hf_exchange_matr_complex)) then
               ! FIXME This is ugly and should not be in this subroutine
               ! Add exchange matrix to Hamiltonian
               if(use_lc_wpbeh) then
                  if(hybrid_coeff /= 0.d0) then
                     ham_complex(:,:,:) = ham_complex(:,:,:)&
                        -((1/lc_dielectric_constant)&
                        *hf_exchange_matr_complex(:,:,1,:)+hybrid_coeff&
                        *hf_exchange_matr_complex_SR(:,:,1,:))
                  else
                     ham_complex(:,:,:) = ham_complex(:,:,:)&
                        -(1/lc_dielectric_constant)&
                        *hf_exchange_matr_complex(:,:,1,:)
                  end if
               else
                  ham_complex(:,:,:) = ham_complex(:,:,:)-hybrid_coeff&
                     *hf_exchange_matr_complex(:,:,1,:)
               end if
            end if

            do i_spin = 1,n_spin
               if(.not. use_elsi) then
                  call solve_evp_scalapack_complex(&
                       KS_eigenvalue(:,i_spin,my_k_point),&
                       KS_eigenvector_complex(:,:,i_spin,1),i_spin)
               else if(use_elsi_ev) then
                  call solve_KS_elsi_parallel(&
                       KS_eigenvalue(:,i_spin,my_k_point),&
                       KS_eigenvector_complex(:,:,i_spin,1),i_spin)

                  if(sing_output) then
                     call aims_elsi_get_ovlp_ev_min(eh_scf,&
                          ovlp_ev_min_k(my_k_point))
                     call aims_elsi_get_ovlp_ev_max(eh_scf,&
                          ovlp_ev_max_k(my_k_point))
                  end if
               else if(use_elsi_dm) then
                  call solve_KS_elsi_dm_parallel(ev_sum,i_spin)
               else
                  write(msg,"(2X,A)") '*** ERROR: Invalid "elsi_method".'
                  call localorb_info(msg,use_unit,"(A)",OL_norm)
                  call aims_stop
               end if
            end do
         else ! Use ScaLAPACK?
            if(packed_matrix_format /= PM_index) then
               call aims_allocate(work_ham,n_centers_basis_I,n_centers_basis_I,&
                    n_spin,"work_ham")
               call aims_allocate(work_ovlp,n_centers_basis_I,&
                    n_centers_basis_I,"work_ovlp")
            else
               ! Dummy allocation
               call aims_allocate(work_ham,1,1,n_spin,"work_ham")
               call aims_allocate(work_ovlp,1,1,"work_ovlp")
            end if

            n_kpt = 0

            do i_kpt = 1,n_k_points
               if((myid == mod(i_kpt,n_tasks) .and. myid <= n_k_points)&
                  .or. use_load_balancing) then
                  n_kpt = n_kpt+1

                  if(flag_rel == REL_x2c .or. flag_rel == REL_4c_dks) then
                     if(n_kpt == 1) then
                        call aims_allocate(ham_work_c,&
                             dim_matrix_rel*(2*dim_matrix_rel+1),n_spin,&
                             "ham_work_c")
                        call aims_allocate(ovlp_work_c,&
                             dim_matrix_rel*(2*dim_matrix_rel+1),"ovlp_work_c")
                     end if

                     call x2c_elsi_lapack_wrapper(n_kpt,i_kpt,n_spin,&
                          ham_work_c,ovlp_work_c,flag_KS_k_points(i_kpt),&
                          n_states_k(i_kpt),KS_eigenvalue,&
                          KS_eigenvector_complex)
                  else
                     if(.not. allocated(ham_work_c)) then
                        call aims_allocate(ham_work_c,n_basis*(n_basis+1)/2,&
                             n_spin,"ham_work_c")
                        call aims_allocate(ovlp_work_c,n_basis*(n_basis+1)/2,&
                             "ovlp_work_c")

                        ! Dummy allocation for real
                        call aims_allocate(ham_work,1,n_spin,"ham_work")
                        call aims_allocate(ovlp_work,1,"ovlp_work")
                     end if

                     call construct_hamiltonian_and_ovl(hamiltonian,&
                          overlap_matrix,ham_work,ovlp_work,ham_work_c,&
                          ovlp_work_c,i_kpt,work_ham,work_ovlp)

                     if(.not. (myid == mod(i_kpt,n_tasks)&
                        .and. myid <= n_k_points)) then
                        n_kpt = n_kpt-1
                        cycle
                     end if

                     ! FIXME This is ugly and should not be in this subroutine
                     ! Add the exchange matrix here if needed
                     if(use_periodic_hf&
                        .and. allocated(hf_exchange_matr_complex)) then
                        if(use_lc_wpbeh) then
                           if(hybrid_coeff /= 0.d0) then
                              call get_hf_hamiltonian_complex_p0(&
                                   (1/lc_dielectric_constant)&
                                   *hf_exchange_matr_complex+hybrid_coeff&
                                   *hf_exchange_matr_complex_SR,ham_work_c,&
                                   n_kpt)
                           else
                              call get_hf_hamiltonian_complex_p0(&
                                   (1/lc_dielectric_constant)&
                                   *hf_exchange_matr_complex,ham_work_c,n_kpt)
                           end if
                        else
                           call get_hf_hamiltonian_complex_p0(&
                                hf_exchange_matr_complex,ham_work_c,n_kpt)
                        end if
                     end if

                     do i_spin = 1,n_spin
                        if(.not. use_elsi) then
                           call improve_complex_eigenfunctions(ovlp_work_c,&
                                ham_work_c(:,i_spin),&
                                KS_eigenvalue(:,i_spin,i_kpt),&
                                KS_eigenvector_complex(:,:,i_spin,n_kpt),i_kpt)
                        else
                           call solve_KS_elsi_serial(ham_work_c(:,i_spin),&
                                ovlp_work_c,KS_eigenvalue(:,i_spin,i_kpt),&
                                KS_eigenvector_complex(:,:,i_spin,n_kpt),&
                                i_spin,i_kpt)

                           if(sing_output) then
                              call aims_elsi_get_ovlp_ev_min(eh_scf,&
                                   ovlp_ev_min_k(i_kpt))
                              call aims_elsi_get_ovlp_ev_max(eh_scf,&
                                   ovlp_ev_max_k(i_kpt))
                           end if
                        end if
                     end do
                  end if ! Fully-relativistic?
               end if ! myid?
            end do ! i_kpt
         end if ! Use ScaLAPACK?

         if(allocated(work_ham)) then
            call aims_deallocate(work_ham,"work_ham")
         end if
         if(allocated(work_ovlp)) then
            call aims_deallocate(work_ovlp,"work_ovlp")
         end if
         if(allocated(ham_work_c)) then
            call aims_deallocate(ham_work_c,"ham_work_c")
         end if
         if(allocated(ovlp_work_c)) then
            call aims_deallocate(ovlp_work_c,"ovlp_work_c")
         end if
         if(allocated(ham_work)) then
            call aims_deallocate(ham_work,"ham_work")
         end if
         if(allocated(ovlp_work)) then
            call aims_deallocate(ovlp_work,"ovlp_work")
         end if
      end if ! Real eigenvectors?
   else ! n_periodic > 0, not periodic, not packed matrices
      ! Dummy setting for output only!
      i_kpt = 1

      if(use_scalapack) then
         ! Normal case: all states treated at once
         if(use_hartree_fock) then
            if(use_periodic_hf) then
               if(allocated(hf_exchange_matr_real)) then
                  ! FIXME This is ugly and should not be in this subroutine
                  ! Add exchange matrix to Hamiltonian
                  if(use_lc_wpbeh) then
                     if(hybrid_coeff /= 0.d0) then
                        ham(:,:,:) = ham(:,:,:)-((1/lc_dielectric_constant)&
                           *hf_exchange_matr_real(:,:,1,:)+hybrid_coeff&
                           *hf_exchange_matr_real_SR(:,:,1,:))
                     else
                        ham(:,:,:) = ham(:,:,:)-(1/lc_dielectric_constant)&
                           *hf_exchange_matr_real(:,:,1,:)
                     end if
                  else
                     ham(:,:,:) = ham(:,:,:)-hybrid_coeff&
                        *hf_exchange_matr_real(:,:,1,:)
                  end if
               end if
            else ! use_periodic_hf
               call setup_hamiltonian_scalapack(hamiltonian)
            end if ! use_periodic_hf
         end if ! use_hartree_fock

         do i_spin = 1,n_spin
            if(.not. use_elsi) then
               call solve_evp_scalapack(KS_eigenvalue(:,i_spin,1),&
                    KS_eigenvector(:,:,i_spin,1),i_spin)
            else if(use_elsi_ev) then
               call solve_KS_elsi_parallel(KS_eigenvalue(:,i_spin,1),&
                    KS_eigenvector(:,:,i_spin,1),i_spin)

               if(sing_output) then
                  call aims_elsi_get_ovlp_ev_min(eh_scf,ovlp_ev_min_k(1))
                  call aims_elsi_get_ovlp_ev_max(eh_scf,ovlp_ev_max_k(1))
               end if
            else if(use_elsi_dm) then
               write(msg,"(2X,A)") '*** ERROR: "elsi_method dm" not supported.'
               call localorb_info(msg,use_unit,"(A)",OL_norm)
               call aims_stop
            else
               write(msg,"(2X,A)") '*** ERROR: Invalid "elsi_method".'
               call localorb_info(msg,use_unit,"(A)",OL_norm)
               call aims_stop
            end if
         end do
      else ! Use ScaLAPACK?
         call aims_allocate(ham_work,n_basis*(n_basis+1)/2,n_spin,"ham_work")
         call aims_allocate(ovlp_work,n_basis*(n_basis+1)/2,"ovlp_work")

         ham_work = hamiltonian
         ovlp_work = overlap_matrix

         ! FIXME This is ugly and should not be in this subroutine
         ! Add the exchange matrix here if needed
         if(use_periodic_hf .and. allocated(hf_exchange_matr_real))then
            if(use_lc_wpbeh) then
               if(hybrid_coeff /= 0.d0) then
                  call get_hf_hamiltonian_real_p0((1/lc_dielectric_constant)&
                       *hf_exchange_matr_real+hybrid_coeff&
                       *hf_exchange_matr_real_SR,ham_work,1)
               else
                  call get_hf_hamiltonian_real_p0((1/lc_dielectric_constant)&
                       *hf_exchange_matr_real,ham_work,1)
               end if
            else
               call get_hf_hamiltonian_real_p0(hf_exchange_matr_real,ham_work,1)
            end if
         end if

         do i_spin = 1,n_spin
            if(.not. use_elsi) then
               lapack_output = .true.

               call improve_real_eigenfunctions(ovlp_work,ham_work(:,i_spin),&
                    lapack_output,KS_eigenvalue(:,i_spin,i_kpt),&
                    KS_eigenvector(:,:,i_spin,i_kpt),i_kpt)
            else
               call solve_KS_elsi_serial(ham_work(:,i_spin),ovlp_work,&
                    KS_eigenvalue(:,i_spin,i_kpt),&
                    KS_eigenvector(:,:,i_spin,i_kpt),i_spin,i_kpt)

               if(sing_output) then
                  call aims_elsi_get_ovlp_ev_min(eh_scf,ovlp_ev_min_k(1))
                  call aims_elsi_get_ovlp_ev_max(eh_scf,ovlp_ev_max_k(1))
               end if
            end if
         end do

         call aims_deallocate(ham_work,"ham_work")
         call aims_deallocate(ovlp_work,"ovlp_work")
      end if ! Use ScaLAPACK?
   end if ! n_periodic > 0

   if(.not. use_elsi_dm) then
      if(n_periodic > 0 .or. packed_matrix_format == PM_index) then
         call sync_eigenvalues(KS_eigenvalue)
         call sync_vector_integer(n_states_k,n_k_points)

         ! Corner case: orbital update + non-periodic geometry + LAPACK
         ! Broadcast eigenvectors
         if(.not. use_density_matrix .and. n_periodic == 0&
            .and. .not. use_scalapack) then
            if(n_tasks > 1 .and. myid /= 1) then
               KS_eigenvector = 0.d0 ! Only k-point 1 contributes
            end if

            call sync_vector(KS_eigenvector,n_basis*n_states*n_spin)
         end if
      end if
   end if

   if(use_elsi) then
      ! ELSI output
      if(myid == 0) then
         call aims_elsi_set_output(eh_scf,elsi_out_level)
      end if

      first_elsi_call = .false.
      out_h_elsi = .false.
      out_s_elsi = .false.

      ! Singular basis output
      if(sing_output) then
         call print_illconditioning(ovlp_ev_min_k,ovlp_ev_max_k,n_states_k)

         sing_output = .false.
      end if
   end if

   ! Set n_states to be the minimal n_nonsingular across k-points
   if(.not. use_elsi_dm .and. elsi_solver /= 5) then
      n_states_new = minval(n_states_k)
      n_states = min(n_states,n_states_new)
      n_states_k = n_states
   end if

   if(allocated(ovlp_ev_min_k)) then
      call aims_deallocate(ovlp_ev_min_k,"ovlp_ev_min_k")
   end if
   if(allocated(ovlp_ev_max_k)) then
      call aims_deallocate(ovlp_ev_max_k,"ovlp_ev_max_k")
   end if

end subroutine
!******
