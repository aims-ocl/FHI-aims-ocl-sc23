!****s* FHI-aims/cpscf_solver_polar_reduce_memory
!  NAME
!    cpscf_solver_polar_reduce_memory
!  SYNOPSIS

    subroutine cpscf_solver_polar_reduce_memory &
    (converged)

!  PURPOSE
!  This is the DFPT_polar calcialtion using the sparse matrix
!  (packed_matrix_format=PM_index). 
!  Comapred to the dense matrix version, there are four points that changed : 
!
!  (1) evaluate_first_order_H_polar_reduce_memory.f90 
!      H1_compute(n_compute,n_compute)        ---> H1_sparse(n_hamiltonian_matrix_size) 
!
!  (2) prune_density_matrix_sparse_polar_reduce_memory.f9 
!       DM1_sparse(n_hamiltonian_matrix_size) ---> DM1_compute(n_compute, n_compute)
!
!  (3) construct_first_order_ham_polar_reduce_memory_scalapack 
!       H1_sparse ---> H1_scalapack_dense
!
!  (4) get_first_order_dm_polar_reduce_memory_scalapack 
!      DM1_scalapack_dense ---> DM1_sparse
! 
! The difference between dense version and the sparse version are listed here: 
!-------------------------------------------------------------------------------
!               DFPT_polarizablity            DFPT_polar_reduce_memeory
!-------------------------------------------------------------------------------
! H1 and DM1  (3, n_basis, n_basis, n_spin)   (n_hamiltonian_matrix_size, n_spin)
! lapack            yes                           no
! scalapack         yes                           yes
! HF                yes                           no
! ------------------------------------------------------------------------------


!  USES

      use dimensions
      use timing
      use physics
      use species_data
      use localorb_io
      use geometry

      use synchronize_mpi_basic, only: sync_vector
      use debugmanager, only: module_is_debugged, debugprint
      use runtime_choices
      use DFPT_pulay_mixing,    only: pulay_mix, cleanup_pulay_mixing
      use mpi_tasks, only: aims_stop, myid
  
      use aims_memory_tracking, only: aims_allocate, aims_deallocate 
      use pbc_lists, only: index_hamiltonian, column_index_hamiltonian 
      ! wyj
      use scalapack_wrapper
      use load_balancing

      use mpi_utilities

      use opencl_util, only: opencl_util_init, load_balance_finished, use_opencl
      use hartree_potential_storage, only: use_rho_multipole_shmem
      implicit none

!  ARGUMENTS

      logical, intent(OUT) :: converged

!  INPUTS
!    none
!  OUTPUT
!    o  converged -- did the cpscf cycle converged or not.
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
!
      ! imported variables


      ! local variables
      character*1000 :: info_str
      logical :: below_it_limit

      character*8  :: cdate
      character*10 :: ctime
      character(*), parameter :: deffmt = '2X'

      ! counters
      integer :: j_coord
      integer :: i_spin
      integer :: i_basis, j_basis

      character(*), parameter :: func = 'cpscf_solver'
      logical :: use_elsi_dm_cpscf
      !logical, parameter ::  use_elsi_dm_cpscf = .true.   ! use NTpoly DM1 matrix solver, scale as  O(N)
      !logical, parameter ::  use_elsi_dm_cpscf = .false. ! use normal DM1 matrix solver, scale as  O(N3) 

      !--------for DFPT debug----------
      !real*8 first_order_H_dense(n_basis, n_basis, n_spin) 
      !real*8 first_order_density_matrix_dense(n_basis, n_basis, n_spin) 
      !integer i_place


      real*8  polarizability(3,3), polar_mean

      real*8  change_of_first_order_DM
      real*8  time_start,time_end

      real*8 :: time1, time_c


  !-------------------begin define------------------------------------
     !----------------(1) grid---------------------------------------
     real*8, allocatable :: first_order_rho(:,:)
     real*8, allocatable :: first_order_total_rho(:)
     real*8, allocatable :: first_order_total_potential(:)
     !----------------(2) matrix---------------------------------------
     real*8, allocatable :: first_order_density_matrix(:,:)
     real*8, allocatable :: old_first_order_density_matrix(:,:)
     real*8, allocatable :: first_order_H(:,:)
  !-------------------end define------------------------------------- 

    !-------wuyangjun begin local index  define variables-----------

    logical first_iteration
    !integer :: d_i1
    ! n_matrix_local_size or n_hamiltonian_matrix_size
    integer :: n_matrix_size
    !-------wuyangjun end local index define variables----------

    !-------wuyangjun begin initialize local index variables----------
    ! wyj: maybe should add in runtime_choice.f90 as a global variable
    first_iteration = .true.

    !-------wuyangjun end initialize local index variables-----

    !call init_prof(0)
    !call start_prof()

  !-------------------begin allocate------------------------------------
     !----------------(1) grid---------------------------------------
     call aims_allocate(first_order_rho, n_spin,     n_full_points, "+first_order_rho")      
     call aims_allocate(first_order_total_rho,       n_full_points, "+first_order_total_rho") 
     call aims_allocate(first_order_total_potential, n_full_points, "+first_order_total_potential")
     !----------------(2) matrix---------------------------------------
     if (.not. (use_local_index .and. use_load_balancing)) then
         call aims_allocate(first_order_density_matrix,     n_hamiltonian_matrix_size, n_spin, "+first_order_density_matrix")
         call aims_allocate(old_first_order_density_matrix, n_hamiltonian_matrix_size, n_spin, "+old_first_order_density_matrix")
         call aims_allocate(first_order_H,                  n_hamiltonian_matrix_size, n_spin, "+first_order_H")
         n_matrix_size = n_hamiltonian_matrix_size
         use_batch_permutation = 0
         get_batch_weights = .false.
     endif
  !-------------------end allocate------------------------------------



     if(myid.eq.0) then
     write(use_unit,*) "-------------------------------------------------------------------------------"
     write(use_unit,*) "|           ENTERING DFPT_POLARIZABILITY (NON-PERIODIC CALCULATION)           |"
     write(use_unit,*) "|                                                                             |"
     write(use_unit,*) "|  Details on the implementation can be found in the following reference:     |"
     write(use_unit,*) "|                                                                             |"
     write(use_unit,*) "|    Honghui Shang, Nathaniel Raimbault, Patrick Rinke,                       |"
     write(use_unit,*) "|     Matthias Scheffler,  Mariana Rossi and Christian Carbogno,              |"
     write(use_unit,*) "|    'All-Electron, Real-Space Perturbation Theory for Homogeneous            |"
     write(use_unit,*) "|    Electric Fields: Theory, Implementation, and Application within dft'     |"
     write(use_unit,*) "|                                                                             |"
     write(use_unit,*) "|  Please cite New Journal of Physics, 20(7):073040, 2018                     |"
     write(use_unit,*) "-------------------------------------------------------------------------------"
     endif

     if (use_opencl) then
        call opencl_init()
     endif

     if (use_local_index .and. use_load_balancing .and. use_gga) &
         call aims_stop("local_index doesn't supports gga")
     
     !-------check if DFPT_polar_reduce_memory can be used for cluster system--------
     if ( .not.( use_scalapack .and. real_eigenvectors ) ) then
        call aims_stop('You are making DFPT_polar_reduce_memeory calcualtion, & 
                        you must use_scalapack. so we stop here.', func)
     end if

     !-------check if we should use ntpoly, we keep the same as DFT part
     if(elsi_solver.eq.6) then 
       use_elsi_dm_cpscf = .true.  ! use NTpoly DM1, O(N)
     else
       use_elsi_dm_cpscf = .false. ! use normal DM1, O(N3)
     endif 

     !-------check if this is a degenerate system---------
     if(n_spin.eq.1.and.mod(n_electrons,2.0d0).ne.0) then 
       write(info_str,'(A)') &
               "                                ^              "
       call localorb_info(info_str, use_unit,'(A)', OL_norm  )
       write(info_str,'(A)') &
               "The system looks like this:   --|--  ----, where the ground state is degenerate (degenerace=2)"
       call localorb_info(info_str, use_unit,'(A)', OL_norm  )
       write(info_str,'(A)') &
               "In order to deal with this system correctly, you need to use  'spin   collinear' and 'default_initial_moment 1'"
       call localorb_info(info_str, use_unit,'(A)', OL_norm  )
        
       call aims_stop('You are making a non-degenerate DFPT calcualtion for a degenerate system, so we stop here.', func) 
     endif

    ! wyj
    if (use_local_index .and. use_load_balancing) then 
        if (n_local_matrix_size .ne. batch_perm(n_bp_integ)%n_local_matrix_size) &
            call aims_stop("cpscf n_local_matrix_size error")
        !n_local_matrix_size = batch_perm(n_bp_integ)%n_local_matrix_size
        call aims_allocate(first_order_density_matrix, n_local_matrix_size, n_spin, "+first_order_density_matrix_sparse_local_index")
        call aims_allocate(old_first_order_density_matrix, n_local_matrix_size, n_spin, "+old_first_order_density_matrix_sparse_local_index")
        call aims_allocate(first_order_H, n_local_matrix_size, n_spin, "+first_order_H_local_index")
        n_matrix_size = n_local_matrix_size
    endif

    !------------------begin polar calculation for each coord---------------------
    do j_coord =1 ,3

      write(info_str,'(A)') ''
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
      write(info_str,'(A)') ''
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
      write(info_str,'(A)') "==========================================================================="
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
      write(info_str,'(2X,A,1X,I4)') 'CPSCF working for j_coord =',j_coord
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      write(info_str,'(A)') ''
      call localorb_info(info_str, use_unit,'(A)', OL_norm )
      write(info_str,'(A)') "=========================================================================="
      call localorb_info(info_str, use_unit,'(A)', OL_norm )
      write(info_str,'(A)') ''
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )


      ! begin work
      converged = .false.
      number_of_loops = 0
      below_it_limit = (number_of_loops.lt.sc_iter_limit)
     

      first_order_H=0.0d0
      first_order_density_matrix=0.0d0
      old_first_order_density_matrix=0.0d0
      first_order_rho=0.0d0
      first_order_total_rho=0.0d0
      first_order_total_potential=0.0d0

       !print *, 'wyj_H1', myid, 'mat_size=', n_local_matrix_size, 'full_points=', n_full_points, &
       !    batch_perm(1)%n_local_matrix_size, batch_perm(1)%n_full_points, &
       !    batch_perm(2)%n_local_matrix_size, batch_perm(2)%n_full_points, &
       !    batch_perm(3)%n_local_matrix_size, batch_perm(3)%n_full_points

    if (use_local_index .and. use_load_balancing) then
        if (first_iteration) then
            get_batch_weights = .true.
        else
            get_batch_weights = .false.
        endif
        ! wyj:TODO debug
        !get_batch_weights = .false.
        use_batch_permutation = n_bp_integ
        !n_local_matrix_size = batch_perm(n_bp_integ)%n_local_matrix_size
    endif

    !print *, 'n_basis=', n_basis, n_centers_basis_I, n_centers_basis_T, n_spin 
    !Starts calculation of U1, which at this point only contains -r from H1
    !wyj: TODO debug
    !if (myid == 0) print *, 'hartree_potential(1:5)=', hartree_potential(1:5)
    !hartree_potential = 1.0d0 ! not be used
    !if (myid == 0) print *, 'rho(1:5)=', rho(1,1:5)
    !rho = 1.0d0 ! not be used

    call integrate_first_order_H_polar_reduce_memory &
        (hartree_potential,first_order_total_potential, rho, rho_gradient,&
        first_order_rho, &
        partition_tab, l_shell_max, &
        j_coord,  & 
        first_order_density_matrix, &
        first_order_H, n_matrix_size &
        )

    !if (myid == 0) print *, myid, 'debug_H=' 
    !! wyj: [PASS]
    !do i_spin = 1, n_spin
    !if (use_local_index .and. use_load_balancing) then
    !    call print_sparse_to_dense_local_index_cpscf(use_batch_permutation,&
    !        first_order_H(:,i_spin), 1)
    !else
    !    if (myid == 0) call print_sparse_to_dense_global_index_cpscf(first_order_H(:,i_spin))
    !endif
    !enddo

       !print *, 'wyj_H2', myid, 'mat_size=', n_local_matrix_size, 'full_points=', n_full_points, &
       !    batch_perm(1)%n_local_matrix_size, batch_perm(1)%n_full_points, &
       !    batch_perm(2)%n_local_matrix_size, batch_perm(2)%n_full_points, &
       !    batch_perm(3)%n_local_matrix_size, batch_perm(3)%n_full_points
    if(use_local_index .and. use_load_balancing ) then
        call init_comm_full_local_matrix(&
            batch_perm(n_bp_integ)%n_basis_local, &
            batch_perm(n_bp_integ)%i_basis_local)

       ! write(use_unit,*) 'myid=', myid, '[H2] n_local_matrix_size',n_local_matrix_size, batch_perm(1)%n_local_matrix_size
        ! wyj: TODO BUG?
        do i_spin = 1, n_spin
        call get_set_full_local_matrix_scalapack_cpscf(first_order_H(:,i_spin), 1, i_spin)
        enddo
        !wyj: scalapack is equal to first_order_H
        !print *, 'print_ham_cpscf'
        !call print_ham_cpscf(first_order_ham_polar_reduce_memory_scalapack)

        if (first_iteration) then
            call reset_batch_permutation(n_bp_integ)
            call compute_balanced_batch_distribution(n_bp_integ)
        endif

        get_batch_weights = .false.
        use_batch_permutation = 0
    endif

       !print *, 'wyj_H3', myid, 'mat_size=', n_local_matrix_size, 'full_points=', n_full_points, &
       !    batch_perm(1)%n_local_matrix_size, batch_perm(1)%n_full_points, &
       !    batch_perm(2)%n_local_matrix_size, batch_perm(2)%n_full_points, &
       !    batch_perm(3)%n_local_matrix_size, batch_perm(3)%n_full_points
    

    !-------shanghui begin debug output------
    if (module_is_debugged("DFPT")) then
     !---------------begin sparse to dense matrix for cluster case--------------------
     ! do i_basis=1,n_basis
     ! do j_basis=1,i_basis ! sparse matrix only store the j_basis < i_basis case
     !    do i_place = index_hamiltonian(1,1, i_basis), index_hamiltonian(2,1,i_basis), 1
     !       if(column_index_hamiltonian(i_place) == j_basis) then
     !          first_order_H_dense(i_basis,j_basis,:) =  & 
     !          first_order_H(i_place, :)
     !          first_order_H_dense(j_basis,i_basis,:) =  & 
     !          first_order_H(i_place, :)
     !       endif
     !    enddo
     ! enddo
     ! enddo
     !---------------end dense to sparse matrix for cluster case--------------------
     !if(myid.eq.0) then
     !do i_spin = 1, n_spin
     !write(use_unit,*) '************shanghui begain first_order_H(X)****************'
     !  do i_basis=1,n_basis
     !  write(use_unit,'(7f15.9)') (first_order_H_dense(i_basis,j_basis,i_spin),j_basis=1,n_basis )
     !  enddo
     !write(use_unit,*) '************shanghui end first_order_H****************'
     !write(use_unit,*) ''
     !enddo
     !endif
    endif
    !-------shanghui end debug output--------



    !print *, myid, 'wyj_before_DM=' 
    if(use_local_index .and. use_load_balancing ) then
        call init_comm_full_local_matrix(&
            batch_perm(n_bp_integ)%n_basis_local, &
            batch_perm(n_bp_integ)%i_basis_local)
        ! wyj
        if (n_local_matrix_size .ne. batch_perm(n_bp_integ)%n_local_matrix_size) call aims_stop("cpscf n_local_matrix_size != batch_perm")

        if (allocated(first_order_density_matrix)) deallocate(first_order_density_matrix)
        if (allocated(old_first_order_density_matrix)) deallocate(old_first_order_density_matrix)
        if (allocated(first_order_H)) deallocate(first_order_H)
        call aims_allocate(first_order_density_matrix, n_local_matrix_size, n_spin, "+first_order_density_matrix_sparse_local_index")
        call aims_allocate(old_first_order_density_matrix, n_local_matrix_size, n_spin, "+old_first_order_density_matrix_sparse_local_index")
        call aims_allocate(first_order_H, n_local_matrix_size, n_spin, "+first_order_H_local_index")
        first_order_density_matrix = 0.d0
        old_first_order_density_matrix = 0.d0
        first_order_H = 0.d0
        n_matrix_size = n_local_matrix_size
    endif

      if(use_elsi_dm_cpscf) then
        call evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory(  &
             first_order_H,  &
             first_order_density_matrix, n_matrix_size)
      else ! not use_elsi_dm_cpscf
        call evaluate_first_order_DM_polar_reduce_memory(  &
             first_order_H,  & 
             KS_eigenvalue, occ_numbers,   &
             first_order_density_matrix, n_matrix_size)
      endif 

    !print *, myid, 'wyj_after_DM=' 
      !! wyj: [PASS]
      !if (myid == 0) print *, myid, 'density_debug'
      !do i_spin = 1, n_spin
      !if (use_local_index .and. use_load_balancing) then
      !    !call print_sparse_to_dense_local_index_cpscf(use_batch_permutation,&
      !    call print_sparse_to_dense_local_index_cpscf(n_bp_integ,&
      !        first_order_density_matrix(:,i_spin), 2)
      !else
      !    if (myid == 0) call print_sparse_to_dense_global_index_cpscf(first_order_density_matrix(:,i_spin))
      !endif
      !enddo

      ! wyj: debug in here

      !if (use_local_index .and. use_load_balancing) then 
      !    if(allocated(first_order_H)) deallocate(first_order_H)
      !endif
      

! ------------------------ self-consistency loop -------------->>
  SCF_LOOP: do while ( (.not.converged) .and.  &
  &                    below_it_limit )
        number_of_loops = number_of_loops + 1

        write(info_str,'(A)') ''
        call localorb_info(info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "-------------------------------------------------------------"
        call localorb_info(info_str, use_unit,'(A)', OL_norm  )
        
          write(info_str,'(10X,A,1X,I4)') "Begin CP-self-consistency iteration #", number_of_loops
        call localorb_info(info_str, use_unit,'(A)', OL_norm )
        
        write(info_str,'(A)') ''
        call localorb_info(info_str, use_unit,'(A)', OL_norm )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, ", Time     :  ", ctime
        call localorb_info(info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "------------------------------------------------------------"
        call localorb_info(info_str, use_unit,'(A)', OL_norm )

       call get_timestamps ( time_cpscf_loop, clock_time_cpscf_loop )

!---------------(1) Begin  update first_order_rho-----------------------------------
          
    !-------shanghui begin debug output------
      if (module_is_debugged("DFPT")) then
       !---------------begin sparse to dense matrix for cluster case--------------------
       !do i_basis=1,n_basis
       !do j_basis=1,i_basis ! sparse matrix only store the j_basis < i_basis case 
       !   do i_place = index_hamiltonian(1,1, i_basis), index_hamiltonian(2,1,i_basis), 1
       !      if(column_index_hamiltonian(i_place) == j_basis) then
       !         first_order_density_matrix_dense(i_basis,j_basis,:) =  & 
       !         first_order_density_matrix(i_place, :)
       !         first_order_density_matrix_dense(j_basis,i_basis,:) =  & 
       !         first_order_density_matrix(i_place, :)
       !      endif
       !   enddo

       !enddo
       !enddo
       !---------------end dense to sparse matrix for cluster case--------------------
       !if(myid.eq.0) then
       !do i_spin = 1,n_spin
       !write(use_unit,*) '************shanghui begain first_order_DM(X)****************'
       !  do i_basis=1,n_basis
       !  write(use_unit,'(7f15.9)') (first_order_density_matrix_dense(i_basis,j_basis,i_spin),j_basis=1,n_basis )
       !  enddo
       !write(use_unit,*) '************shanghui end first_order_DM****************'
       !write(use_unit,*) ''
       !enddo
       !endif
      endif
    !-------shanghui end debug output--------





       
       !---(1.1): we need to perform density matrix mixing here
         change_of_first_order_DM =0.0d0         

        !do i_basis = 1, n_hamiltonian_matrix_size - 1
        ! wyj: TODO why -1?
        do i_basis = 1, n_matrix_size - 1
        do i_spin = 1, n_spin

           change_of_first_order_DM =                &
           max( change_of_first_order_DM,             &
           dabs(first_order_density_matrix(i_basis, i_spin)  &
              - old_first_order_density_matrix(i_basis, i_spin)) )

        enddo ! i_spin
        enddo ! i_basis



         converged = (change_of_first_order_DM.lt.DFPT_sc_accuracy_dm).and.(number_of_loops.gt.1) 

         if (.not.converged) then
             ! wyj: TODO debug
             use_dfpt_pulay = .false.

             if(use_dfpt_pulay) then ! pulay mixing   

                 !---begin shanghui note------
                 !Here the pulay_mix coefficients depend on the size of the input matrix, 
                 !so the DFPT_polar_reduce_meomry result is a little bit diffirent from 
                 !the DFPT_polarizabily when using pulay_mix 
                 !---end shanghui note------
                 call pulay_mix(first_order_density_matrix, size(first_order_density_matrix), number_of_loops, dfpt_pulay_steps, DFPT_mixing)

                 old_first_order_density_matrix =  first_order_density_matrix

             else ! linear mixing

                 first_order_density_matrix =       &
                     (1.0d0-DFPT_mixing)*old_first_order_density_matrix+  &
                     DFPT_mixing*first_order_density_matrix

                 old_first_order_density_matrix = first_order_density_matrix
             endif ! mixing
 

          !-------shanghui begin debug_mode------
          if (module_is_debugged("DFPT")) then
          write(info_str,'(A)') "((((((((((((((((((((((((((((((((("
          call localorb_info(info_str, use_unit,'(A)', OL_norm )
          write(info_str,*) change_of_first_order_DM
          call localorb_info(info_str, use_unit,'(A)', OL_norm )
          write(info_str,'(A)') ")))))))))))))))))))))))))))))))))"
          call localorb_info(info_str, use_unit,'(A)', OL_norm )
          endif
          !-------shanghui end debug_mode--------
          !wyj :debug
          write(info_str,'(A)') "((((((((((((((((((((((((((((((((("
          call localorb_info(info_str, use_unit,'(A)', OL_norm )
          write(info_str,*) change_of_first_order_DM
          call localorb_info(info_str, use_unit,'(A)', OL_norm )
          write(info_str,'(A)') ")))))))))))))))))))))))))))))))))"
          call localorb_info(info_str, use_unit,'(A)', OL_norm )

          write(info_str,'(A)') ''
          call localorb_info(info_str, use_unit,'(A)', OL_norm  )

          write(info_str,'(2X,A)') &
          "CPSCF convergence accuracy:"
          call localorb_info(info_str, use_unit,'(A)', OL_norm  )
          write(info_str,'(2X,A,1X,E10.4,1X,E10.4)') &
                  "| Change of first_order_density_matrix     :", change_of_first_order_DM
          call localorb_info(info_str, use_unit,'(A)', OL_norm  )


       call get_timestamps(time_first_order_density, clock_time_first_order_density)

       ! wyj local_index: rho
       if (use_local_index .and. use_load_balancing) then
           if(first_iteration) then
               ! We use the batch permutation from the initial integration,
               ! which should be better than the default batch permutation,
               ! to calculate the weights for the density calculation
               get_batch_weights = .true.
               use_batch_permutation = n_bp_integ
           else
               ! Use the already-calculated batch permutation for density
               get_batch_weights = .false.
               !use_batch_permutation = n_bp_density
               ! wyj:TODO debug
               use_batch_permutation = n_bp_integ
           endif
           ! wyj:TODO debug
           !get_batch_weights = .false.
           !use_batch_permutation = n_bp_integ
       endif


      ! print *, 'wyj_1', myid, 'mat_size=', n_local_matrix_size, 'full_points=', n_full_points, &
      !     batch_perm(1)%n_local_matrix_size, batch_perm(1)%n_full_points, &
      !     batch_perm(2)%n_local_matrix_size, batch_perm(2)%n_full_points, &
      !     batch_perm(3)%n_local_matrix_size, batch_perm(3)%n_full_points

       !if (n_local_matrix_size .ne. batch
       !---(1.2) calculate first_order_rho
       do i_spin = 1, n_spin
       call integrate_first_order_rho_polar_reduce_memory(partition_tab, l_shell_max,  &
           first_order_density_matrix(:,i_spin), &
           first_order_rho(i_spin,:), size(first_order_density_matrix, dim=1))
       enddo ! i_spin
      ! print *, 'wyj_2', myid, 'mat_size=', n_local_matrix_size, 'full_points=', n_full_points, &
      !     batch_perm(1)%n_local_matrix_size, batch_perm(1)%n_full_points, &
      !     batch_perm(2)%n_local_matrix_size, batch_perm(2)%n_full_points, &
      !     batch_perm(3)%n_local_matrix_size, batch_perm(3)%n_full_points

      !if (use_local_index .and. use_load_balancing) then
      !    !use_batch_permutation = n_bp_hpot
      !    ! wyj: TODO debug
      !    get_batch_weights = .false.
      !    use_batch_permutation = n_bp_integ
      !endif
      !call integrate_polar_reduce_memory &
      !    (partition_tab,first_order_rho,polarizability(1:3,j_coord))
      !if(myid == 0) print *, myid, 'wyj_polar=', polarizability(1:3,j_coord)

       !print *, 'CPSCF_rho =', first_order_rho(i_spin,:)
       !print *, 'exit first_order_rho-------------'
       !if(use_local_index .and. use_load_balancing .and. first_iteration) then
       !    call reset_batch_permutation(n_bp_density)
       !    call compute_balanced_batch_distribution(n_bp_density)
       !endif
       ! wyj: reset, avoid affect after function
       get_batch_weights = .false.
       use_batch_permutation = 0

       !print *, 'wyj_3', myid, 'mat_size=', n_local_matrix_size, 'full_points=', n_full_points, &
       !    batch_perm(1)%n_local_matrix_size, batch_perm(1)%n_full_points, &
       !    batch_perm(2)%n_local_matrix_size, batch_perm(2)%n_full_points, &
       !    batch_perm(3)%n_local_matrix_size, batch_perm(3)%n_full_points

       call get_times(time_first_order_density, clock_time_first_order_density, &
           &              tot_time_first_order_density, tot_clock_time_first_order_density)

!------------(1) end first-order-density update and mixing--------


!------------(2)begain to calculate first_order_H-----------------
       call get_timestamps(time_first_order_potential, clock_time_first_order_potential)

           first_order_total_rho = 0.0d0
           do i_spin = 1, n_spin
              first_order_total_rho(1:n_full_points) = & 
              first_order_total_rho(1:n_full_points) + &
              first_order_rho(i_spin, 1:n_full_points)  
           enddo
           !print *, myid, 'SUM(rho)=', sum(first_order_rho(1,:)), first_order_rho(1,1:5)
           !call init_prof(0)
           !call start_prof()

           time1 = mpi_wtime()

        if (use_rho_multipole_shmem) then
           call update_hartree_potential_p2_shanghui &
               ( hartree_partition_tab,first_order_total_rho(1:n_full_points),& 
               delta_v_hartree_part_at_zero, &
               delta_v_hartree_deriv_l0_at_zero, &
               multipole_moments, multipole_radius_sq, &
               l_hartree_max_far_distance, &
               outer_potential_radius )
        else
           call update_hartree_potential_p2_shanghui_no_shmem &
               ( hartree_partition_tab,first_order_total_rho(1:n_full_points),& 
               delta_v_hartree_part_at_zero, &
               delta_v_hartree_deriv_l0_at_zero, &
               multipole_moments, multipole_radius_sq, &
               l_hartree_max_far_distance, &
               outer_potential_radius )
        endif

           time_c = mpi_wtime() - time1
           if(myid .eq. 0) print*, "myid=", myid, " time of update_hartree_potential_p2_shanghui = ", time_c

           ! wyj: local index begin sum_up_whole_potential...
           if (use_local_index .and. use_load_balancing) then
               if (first_iteration) then
                   get_batch_weights = .true.
                   use_batch_permutation = n_bp_integ
               else
                   get_batch_weights = .false.
                   use_batch_permutation = n_bp_hpot
               endif
           ! wyj: TODO debug
           !get_batch_weights = .false.
           !use_batch_permutation = n_bp_integ
           endif


           !call sum_up_whole_potential_p2_shanghui &
           !     ( delta_v_hartree_part_at_zero, &
           !     delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
           !     partition_tab, first_order_total_rho(1:n_full_points), &
           !     first_order_total_potential(1:n_full_points),  & !<--------get first_order_DM_potential
           !     .false., multipole_radius_sq, &
           !     l_hartree_max_far_distance, &
           !     outer_potential_radius)
           call sum_up_whole_potential_shanghui_dielectric &
                ( delta_v_hartree_part_at_zero, &
                delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
                partition_tab, first_order_total_rho(1:n_full_points), &
                first_order_total_potential(1:n_full_points),  & !<--------get first_order_DM_potential
                multipole_radius_sq, &
                l_hartree_max_far_distance, &
                outer_potential_radius)

            if(use_local_index .and. use_load_balancing .and. first_iteration) then
                call reset_batch_permutation(n_bp_hpot)
                call compute_balanced_batch_distribution(n_bp_hpot)
            endif

            get_batch_weights = .false.
            use_batch_permutation = 0

            !call pause_prof()
            !call stop_prof()
       call get_times(time_first_order_potential, clock_time_first_order_potential, &
        &              tot_time_first_order_potential, tot_clock_time_first_order_potential)
    first_iteration = .false.

       call get_timestamps(time_first_order_H, clock_time_first_order_H)

       if (use_local_index .and. use_load_balancing) then
           use_batch_permutation = n_bp_integ
       endif

       time1 = mpi_wtime()

       call integrate_first_order_H_polar_reduce_memory &
           (hartree_potential,first_order_total_potential, rho, rho_gradient,&
           first_order_rho, &
           partition_tab, l_shell_max, &
           j_coord, & 
           first_order_density_matrix, &
           first_order_H, n_matrix_size )

       time_c = mpi_wtime() - time1
       if(myid .eq. 0) print*, "myid=", myid, " time of integrate_first_order___H_polar_reduce_memory = ", time_c


       !if (myid == 0) print *, myid, '2_debug_H=' 
       !! wyj: 
       !do i_spin = 1, n_spin
       !if (use_local_index .and. use_load_balancing) then
       !    call print_sparse_to_dense_local_index_cpscf(use_batch_permutation,&
       !        first_order_H(:,i_spin), 1)
       !else
       !    if (myid == 0) call print_sparse_to_dense_global_index_cpscf(first_order_H(:,i_spin))
       !endif
       !enddo

       time1 = mpi_wtime()

    if(use_local_index .and. use_load_balancing ) then
        call init_comm_full_local_matrix(&
            batch_perm(n_bp_integ)%n_basis_local, &
            batch_perm(n_bp_integ)%i_basis_local)

        ! wyj: TODO BUG?
        do i_spin = 1, n_spin
        call get_set_full_local_matrix_scalapack_cpscf(first_order_H(:,i_spin), 1, i_spin)
        enddo
        !wyj: scalapack is equal to first_order_H
        !print *, 'print_ham_cpscf'
        !call print_ham_cpscf(first_order_ham_polar_reduce_memory_scalapack)

        get_batch_weights = .false.
        use_batch_permutation = 0
    endif

       time_c = mpi_wtime() - time1
       if(myid .eq. 0) print*, "myid=", myid, " time of get_set_full_local_matrix_scalapack_cpscf = ", time_c

       call get_times(time_first_order_H, clock_time_first_order_H, &
        &              tot_time_first_order_H, tot_clock_time_first_order_H)

!------------(2) end to calculate first_order_H-----------------


!------------(3) begin calculation of first_order_DM-----------------
       call get_timestamps(time_first_order_DM, clock_time_first_order_DM)

       if(use_local_index .and. use_load_balancing ) then
           call init_comm_full_local_matrix(&
               batch_perm(n_bp_integ)%n_basis_local, &
               batch_perm(n_bp_integ)%i_basis_local)
       endif

       if(myid .eq. 0) print*, "use_elsi_dm_cpscf=", use_elsi_dm_cpscf
       if(use_elsi_dm_cpscf) then
         call evaluate_first_order_DM_elsi_dm_cpscf_polar_reduce_memory(  &
              first_order_H,  &
              first_order_density_matrix, n_matrix_size)
       else ! not use_elsi_dm_cpscf
         call evaluate_first_order_DM_polar_reduce_memory(  &
              first_order_H,  & 
              KS_eigenvalue, occ_numbers,   &
              first_order_density_matrix, n_matrix_size)
       endif 
      !! wyj: [NOPASS]
      !if (myid == 0) print *, myid, '2_density_debug'
      !do i_spin = 1, n_spin
      !if (use_local_index .and. use_load_balancing) then
      !    !call print_sparse_to_dense_local_index_cpscf(use_batch_permutation,&
      !    call print_sparse_to_dense_local_index_cpscf(n_bp_integ,&
      !        first_order_density_matrix(:,i_spin), 2)
      !else
      !    if (myid == 0) call print_sparse_to_dense_global_index_cpscf(first_order_density_matrix(:,i_spin))
      !endif
      !enddo

       call get_times(time_first_order_DM, clock_time_first_order_DM, &
        &              tot_time_first_order_DM, tot_clock_time_first_order_DM)
!--------------(3) end first_order_DM ----------------

! --------- Check convergence ----->>

!         Check convergence.
!         Continue with density update and new Hartree potential.
!         Get total energy pieces for next iteration.


!         check convergence of self-consistency loop
         
       end if ! if (.not.converged)

        !converged = (change_of_first_order_DM.lt.DFPT_sc_accuracy_dm).and.(number_of_loops.gt.1) 

!  ---------- Update electron density and perform mixing ---->>

        if (converged) then
!           We are done - no further evaluation of density / potential needed

          write(info_str,'(A)') ''
          call localorb_info(info_str, use_unit,'(A)',OL_norm)
          write(info_str,'(2X,A,I5,A)') "CP-self-consistency cycle converged in", number_of_loops-1, " iterations"
          call localorb_info(info_str, use_unit,'(A)',OL_norm)
          write(info_str,'(A)') ''
          call localorb_info(info_str, use_unit,'(A)',OL_norm)
 


        else if (number_of_loops.ge.sc_iter_limit) then
!           This was the last self-consistency cycle - we do not need any more potential / density evaluations

          below_it_limit = .false.
        end if

        call get_times(time_cpscf_loop, clock_time_cpscf_loop)

        ! current SCF loop ends here

! ----- Printing out time data -- >>
        ! shanghui note: the time of the last iteration is zero, so do not rely on it.    
        write(info_str,'(A,I5)') &
        & "End CPSCF iteration # ", number_of_loops 
        call output_timeheader(deffmt, info_str, OL_norm)
        call output_times(deffmt, "Time for this iteration", &
        &                 time_cpscf_loop, clock_time_cpscf_loop, OL_norm)
        call output_times(deffmt, "first_order_density", &
          &                 time_first_order_density, clock_time_first_order_density, OL_norm)
        call output_times(deffmt, "first_order_potential", &
        &                 time_first_order_potential, clock_time_first_order_potential, OL_norm)
        call output_times(deffmt, "first_order_H", &
        &                 time_first_order_H, clock_time_first_order_H, OL_norm)
        call output_times(deffmt, "first_order_DM", &
        &                 time_first_order_DM, clock_time_first_order_DM, OL_norm)
        call output_times(deffmt, "Solution of Sternheimer eqns.", &
        &                 time_Sternheimer, clock_time_Sternheimer, OL_norm)
        write(info_str,'(A)') &
        "------------------------------------------------------------"
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

! << ---- end printing out data---------

        load_balance_finished = .true.
        if(load_balance_finished) opencl_util_init = .true.
        
!       this is the end of the self-consistency cycle.

  end do SCF_LOOP
! << ------ end self consistent cycle--------

      total_number_of_loops = total_number_of_loops + number_of_loops

      if (use_local_index .and. use_load_balancing) then
          get_batch_weights = .false.
          use_batch_permutation = n_bp_hpot
      endif

     call  integrate_polar_reduce_memory &
           (partition_tab,first_order_rho,polarizability(1:3,j_coord))

      if (use_local_index .and. use_load_balancing) then
          get_batch_weights = .false.
          use_batch_permutation = 0
      endif

     if(use_dfpt_pulay) then
     call cleanup_pulay_mixing()
     endif
   
 
   enddo ! j_coord

!    print *, myid, 'wyj_output_cube'
!    call setup_cube_defaults
!    call output_cube_files_p2_cpscf(first_order_density_matrix(:,1))

   !call pause_prof()
   !call stop_prof()

   call get_timestamps(time_Hessian, clock_time_Hessian)

        polar_mean = (polarizability(1,1) + polarizability(2,2) + polarizability(3,3))/3.0

        write(info_str,'(A)') 'DFPT for polarizability:--->'
        call localorb_info(info_str, use_unit,'(A)', OL_norm)
        do j_coord=1,3
         write(info_str,*) polarizability(j_coord,1:3)
         call localorb_info(info_str, use_unit,'(A)', OL_norm)
        enddo
        ! write(info_str,'(A,F11.6)') "The mean polarizability is ", polar_mean
        write(info_str,*) "The mean polarizability is ", polar_mean
        call localorb_info(info_str, use_unit,'(A)', OL_norm )
        write (info_str,'(A)') 'DFPT polarizability (Bohr^3)        xx        yy        zz        xy        xz        yz'
        call localorb_info(info_str, use_unit,'(A)', OL_norm)
        ! write (info_str,'(2X,A,1X,6F10.5)') &
        write (info_str,*) &
        '| Polarizability:--->          ', polarizability(1,1), polarizability(2,2), &
        polarizability(3,3), polarizability(1,2), polarizability(1,3), polarizability(2,3)
        
        ! if(myid .eq. 0) then

        !     print*, polar_mean
        !     print*, polarizability(1,1), polarizability(2,2), polarizability(3,3), &
        !             polarizability(1,2), polarizability(1,3), polarizability(2,3)
        ! endif

        call localorb_info(info_str, use_unit,'(A)', OL_norm)
    call get_times(time_Hessian, clock_time_Hessian, &
        &              tot_time_Hessian, tot_clock_time_Hessian)

    write(info_str,'(A)') ''
    call localorb_info(info_str, use_unit,'(A)', OL_norm  )
    call output_timeheader(deffmt, info_str, OL_norm)
    call output_times(deffmt, "Time for polarizability calculation", &
        &                 time_Hessian, clock_time_Hessian, OL_norm)

    write(info_str,'(A)') "==========================================================================="
    call localorb_info(info_str, use_unit,'(A)', OL_norm  )

      ! Perform any post-processing that is required after every scf cycle:
      ! * scaled ZORA, if required
      ! * output of a band structure
      ! * Direct calculation of a binding energy
      !
      !!! Please add any other postprocessing right here in the future !!!

!      require_post_prc = ( (flag_rel.eq.REL_zora) .and. (force_potential.eq.0) ) &
!                         .or. out_band .or. use_qmmm &
!                         .or. out_hirshfeld .or. out_hirshfeld_iterative &
!                         .or. use_vdw_correction_hirshfeld &
!                         .or. use_ll_vdwdf.or. flag_compute_kinetic .or. use_meta_gga &
!                         .or. use_nlcorr_post .or. use_vdw_post
!
!      
!      if (require_post_prc) then
!      end if

  !-------------------begin deallocate------------------------------------
     !----------------(1) grid---------------------------------------
     deallocate(first_order_rho)
     deallocate(first_order_total_rho)
     deallocate(first_order_total_potential)
     !----------------(2) matrix---------------------------------------
     deallocate(first_order_density_matrix)
     deallocate(old_first_order_density_matrix)
     deallocate(first_order_H)
  !-------------------end deallocate------------------------------------



    end subroutine cpscf_solver_polar_reduce_memory
!******
