!****s* FHI-aims/cpscf_solver_polarizability
!  NAME
!    cpscf_solver_polarizability
!  SYNOPSIS

    subroutine cpscf_solver_polarizability &
    (converged)

!  PURPOSE
!  an cpscf process.
!  USES

      use dimensions
      use timing
      use physics
      use species_data
      use localorb_io
      use geometry
      use scalapack_wrapper, only : eigenvec, my_row, my_col,  &
                                    construct_first_order_ham_polar_scalapack, &
                                    evaluate_first_order_U_polar_scalapack, & 
                                    construct_hamiltonian_real_for_elsi_scalapack, & 
                                    construct_overlap_real_for_elsi_scalapack,     &
                                    set_full_matrix_real, &
                                    mxld, mxcol, &
                                    construct_first_order_hamiltonian_polar_for_elsi_scalapack, & 
                                    get_first_order_dm_polar_for_elsi_scalapack

      use synchronize_mpi_basic, only: sync_vector
      use debugmanager, only: module_is_debugged, debugprint
      use runtime_choices
      use DFPT_pulay_mixing,    only: pulay_mix, cleanup_pulay_mixing
      use mpi_tasks, only: aims_stop, myid
   
      use aims_memory_tracking, only: aims_allocate, aims_deallocate 
      use elsi, only: elsi_dm_real, elsi_dm_real_cpscf
      use elsi_wrapper, only: eh_scf




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
      integer :: i_spin
      integer :: i_atom
      integer :: i_basis, j_basis, i_coord,i_coord_2, i_point,i_state

      character(*), parameter :: func = 'cpscf_solver'

      logical :: use_elsi_dm_cpscf
      !logical, parameter ::  use_elsi_dm_cpscf = .true.
      !logical, parameter ::  use_elsi_dm_cpscf = .false.

      real*8  polarizability(3,3), polar_mean

      real*8  change_of_first_order_DM
      real*8  time_start,time_end
      integer  max_occ_number(n_spin)

      real*8, external :: ddot

  !-------------------begin define------------------------------------
     !----------------(1) grid---------------------------------------
     real*8, allocatable :: first_order_rho(:,:,:)
     real*8, allocatable :: first_order_total_rho(:)
     real*8, allocatable :: first_order_total_potential(:,:)
     !----------------(2) matrix---------------------------------------
     real*8, allocatable :: first_order_density_matrix(:,:,:,:)
     real*8, allocatable :: old_first_order_density_matrix(:,:,:,:)
     real*8, allocatable :: first_order_H(:,:,:,:)
     real*8, allocatable :: first_order_U(:,:,:,:)
     real*8, allocatable :: first_order_E(:,:,:)
     real*8, allocatable :: matrix_tmp(:,:,:,:)
  !-------------------end define------------------------------------- 

     real*8, dimension(:,:,:), allocatable   :: ham
     real*8, dimension(:,:),   allocatable   :: ovlp
     real*8, dimension(:,:,:), allocatable   :: dm
     real(8)  :: energy
     real*8, dimension(:,:,:), allocatable   :: first_order_ham
     real*8, dimension(:,:),   allocatable   :: first_order_ovlp
     real*8, dimension(:,:,:), allocatable   :: first_order_dm



  !-------------------begin allocate------------------------------------
     !----------------(1) grid---------------------------------------
     allocate(first_order_rho(3,n_spin,n_full_points))
     allocate(first_order_total_rho(n_full_points))
     allocate(first_order_total_potential(3,n_full_points))
     !----------------(2) matrix---------------------------------------
     allocate(first_order_density_matrix(3, n_basis,n_basis, n_spin))
     allocate(old_first_order_density_matrix(3, n_basis,n_basis, n_spin))
     allocate(first_order_H(3, n_basis,n_basis, n_spin))
     allocate(first_order_U(3, n_states, n_states, n_spin))
     allocate(first_order_E(3, n_states, n_spin))
     allocate(matrix_tmp(3,size(eigenvec,1),size(eigenvec,2),n_spin))
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



!---------------------------------------------------------------
!------------shanghui begin test call for NTPoly in cpscf-------
!   if(use_elsi_dm_cpscf) then 
!    if(use_scalapack) then ! we only have scalapack version
!
!     call aims_allocate(ham,  mxld, mxcol, n_spin, "+ham")
!     call aims_allocate(ovlp, mxld, mxcol, "+ovlp")
!     call aims_allocate(dm,  mxld, mxcol, n_spin, "+dm")
!
!     if(real_eigenvectors) then
!        write(6,*) '=====================begin test_in_polar_before_cpscf=====================================' 
!        write(6,*) 'shanghui test in polar before elsi(H):', hamiltonian 
!        write(6,*) 'shanghui test in polar before elsi(S):', overlap_matrix
!        write(6,*) '-----------------------------------'
!        write(6,*) 'occ_number:', occ_numbers
!        write(6,*) 'KS_eigen:',  KS_eigenvalue       
!
!
!        call construct_hamiltonian_real_for_elsi_scalapack(hamiltonian, ham)
!        call construct_overlap_real_for_elsi_scalapack(overlap_matrix, ovlp)
!        do i_spin = 1, n_spin   
!           call set_full_matrix_real(ham(:,:,i_spin))
!           call set_full_matrix_real(ovlp)
!           call elsi_dm_real(eh_scf,ham(:,:,i_spin),ovlp,dm(:,:,i_spin),energy)
!        enddo
!
!     end if
!     write(6,*) 'shanghui test in DFPT_polar after elsi:',energy 
!     write(6,*) 'test:',converged, real_eigenvectors
!     write(6,*) '=====================end test_in_polar_before_cpscf=====================================' 
!     write(6,*) '      ' 
!     if(allocated(ham))  call aims_deallocate(ham, "+ham")
!     if(allocated(ovlp)) call aims_deallocate(ovlp, "+ovlp")
!     if(allocated(dm))   call aims_deallocate(dm, "+dm")
!    endif  
!   endif
!------------shanghui end call for NTPoly in cpscf--------------
!---------------------------------------------------------------



      !-------check if we should use ntpoly, we keep the same as DFT part
      if(elsi_solver.eq.6) then
        use_elsi_dm_cpscf = .true.  ! use NTpoly DM1, O(N)
      else
        use_elsi_dm_cpscf = .false. ! use normal DM1, O(N3)
      endif



      ! begin work
      number_of_loops = 0

      below_it_limit = (number_of_loops.lt.sc_iter_limit)
     
      !check if this is a degenerate system.
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


      !find the max_occ_number
        do i_spin = 1, n_spin, 1
            max_occ_number(i_spin) = 0
            do i_state = n_states, 1, -1
             if (dabs(occ_numbers(i_state,i_spin,1)).gt.0.d0) then
              max_occ_number(i_spin) = i_state
              exit
             endif
            enddo
        enddo ! i_spin



    first_order_density_matrix=0.0d0
    old_first_order_density_matrix=0.0d0
    first_order_H=0.0d0
    first_order_U=0.0d0 
    first_order_E=0.0d0
    first_order_rho=0.0d0
    first_order_total_potential=0.0d0

!Starts calculation of U1, which at this point only contains -r from H1

        call  integrate_first_order_H_polarizability &
        (hartree_potential,first_order_total_potential, rho, rho_gradient,&
         first_order_rho, &
         partition_tab, l_shell_max, & 
         first_order_density_matrix, &
         first_order_H &
       )

       !-------shanghui begin parallel------
     if (module_is_debugged("DFPT")) then
       if(myid.eq.0) then
       do i_spin = 1, n_spin
       write(use_unit,*) '************shanghui begain first_order_H(X)****************'
        do i_basis=1,n_basis
        write(use_unit,'(7f15.9)') (first_order_H(1,i_basis,j_basis,i_spin),j_basis=1,n_basis )
        enddo
       write(use_unit,*) '************shanghui end first_order_H****************'
       write(use_unit,*) ''
       enddo
       endif 
     endif
       !-------shanghui end parallel--------

!---------#ELSI----------------------------
    if(use_elsi_dm_cpscf) then 

     if(use_scalapack) then ! we only have scalapack version
 
      call aims_allocate(ham,  mxld, mxcol, n_spin, "+ham")
      call aims_allocate(ovlp, mxld, mxcol, "+ovlp")
      call aims_allocate(dm,  mxld, mxcol, n_spin, "+dm")
      call aims_allocate(first_order_ham,  mxld, mxcol, n_spin, "+first_order_ham")
      call aims_allocate(first_order_ovlp, mxld, mxcol, "+first_order_ovlp")
      call aims_allocate(first_order_dm,  mxld, mxcol, n_spin, "+first_order_dm")

      if(real_eigenvectors) then

         call construct_overlap_real_for_elsi_scalapack(overlap_matrix, ovlp)
         call construct_hamiltonian_real_for_elsi_scalapack(hamiltonian, ham)

        do i_coord = 1, 3 

          !-------1. global matrix to scalapack matrix--------
            first_order_ovlp = 0.0d0
            call construct_first_order_hamiltonian_polar_for_elsi_scalapack(first_order_H(i_coord,:,:,:), first_order_ham)
         do i_spin = 1, n_spin   
            call set_full_matrix_real(ovlp)
            call set_full_matrix_real(ham(:,:,i_spin))
            call set_full_matrix_real(first_order_ham(:,:,i_spin))
          !-------2. get first_order_DM----------
            call elsi_dm_real_cpscf(eh_scf,ham(:,:,i_spin),ovlp,dm(:,:,i_spin), energy, & 
               first_order_ham(:,:,i_spin), first_order_ovlp, first_order_dm(:,:,i_spin)  )
          !-------3. scalapack matrix to global matrix---------
            call get_first_order_dm_polar_for_elsi_scalapack( first_order_dm(:,:,i_spin), first_order_density_matrix(i_coord,:,:,i_spin))
         enddo ! i_spin

        enddo ! i_coord
      endif   ! real_eigen

      if(allocated(ham))  call aims_deallocate(ham, "+ham")
      if(allocated(ovlp)) call aims_deallocate(ovlp, "+ovlp")
      if(allocated(dm))   call aims_deallocate(dm, "+dm")
      if(allocated(first_order_ham))  call aims_deallocate(first_order_ham, "+first_order_ham")
      if(allocated(first_order_ovlp)) call aims_deallocate(first_order_ovlp, "+first_order_ovlp")
      if(allocated(first_order_dm))   call aims_deallocate(first_order_dm, "+first_order_dm")

     endif    ! use_scalapack

    else ! not use_elsi_dm_cpscf

        if (use_scalapack) then
           ! i) first transform dense first-order H lapack to dense first-order H scalapack

           call construct_first_order_ham_polar_scalapack(first_order_H)
             
           !----for debug--------
           !  call  get_first_order_ham_polar_scalapack(first_order_H,1)

           ! ii) evaluate first-order U scalapack

           call evaluate_first_order_U_polar_scalapack(occ_numbers,KS_eigenvalue)

           ! iii) for debug 
           ! call get_first_order_U_polar_scalapack(first_order_U, 1)


        else ! lapack version
           call evaluate_first_order_U_polarizability(first_order_H,  &
             KS_eigenvector, KS_eigenvalue, occ_numbers,  max_occ_number, &
             first_order_U,first_order_E)
        end if

     !call evaluate_first_order_DM_polarizability( &
     !     KS_eigenvector, KS_eigenvalue, occ_numbers,max_occ_number,   &
     !     first_order_U,old_first_order_density_matrix)

    endif ! use_elsi_dm_cpscf 




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
       call get_timestamps(time_first_order_DM, clock_time_first_order_DM)


!---------#ELSI----------------------------
    if(use_elsi_dm_cpscf) then 

     if(use_scalapack) then ! we only have scalapack version
 
      call aims_allocate(ham,  mxld, mxcol, n_spin, "+ham")
      call aims_allocate(ovlp, mxld, mxcol, "+ovlp")
      call aims_allocate(dm,  mxld, mxcol, n_spin, "+dm")
      call aims_allocate(first_order_ham,  mxld, mxcol, n_spin, "+first_order_ham")
      call aims_allocate(first_order_ovlp, mxld, mxcol, "+first_order_ovlp")
      call aims_allocate(first_order_dm,  mxld, mxcol, n_spin, "+first_order_dm")

      if(real_eigenvectors) then

         call construct_overlap_real_for_elsi_scalapack(overlap_matrix, ovlp)
         call construct_hamiltonian_real_for_elsi_scalapack(hamiltonian, ham)

        do i_coord = 1, 3 
         
          !-------1. global matrix to scalapack matrix--------
            first_order_ovlp = 0.0d0
            call construct_first_order_hamiltonian_polar_for_elsi_scalapack(first_order_H(i_coord,:,:,:), first_order_ham)
         do i_spin = 1, n_spin   
            call set_full_matrix_real(ovlp)
            call set_full_matrix_real(ham(:,:,i_spin))
            call set_full_matrix_real(first_order_ham(:,:,i_spin))
          !-------2. get first_order_DM----------
            call elsi_dm_real_cpscf(eh_scf,ham(:,:,i_spin),ovlp,dm(:,:,i_spin), energy, & 
               first_order_ham(:,:,i_spin), first_order_ovlp, first_order_dm(:,:,i_spin)  )
          !-------3. scalapack matrix to global matrix---------
            call get_first_order_dm_polar_for_elsi_scalapack( first_order_dm(:,:,i_spin), first_order_density_matrix(i_coord,:,:,i_spin))
         enddo ! i_spin
        enddo ! i_coord
            call sync_vector(first_order_density_matrix,3*n_basis*n_basis*n_spin) 
      endif   ! real_eigen

      if(allocated(ham))  call aims_deallocate(ham, "+ham")
      if(allocated(ovlp)) call aims_deallocate(ovlp, "+ovlp")
      if(allocated(dm))   call aims_deallocate(dm, "+dm")
      if(allocated(first_order_ham))  call aims_deallocate(first_order_ham, "+first_order_ham")
      if(allocated(first_order_ovlp)) call aims_deallocate(first_order_ovlp, "+first_order_ovlp")
      if(allocated(first_order_dm))   call aims_deallocate(first_order_dm, "+first_order_dm")


     endif    ! use_scalapack




    else ! not use_elsi_dm_cpscf

       call evaluate_first_order_DM_polarizability(  &
             KS_eigenvector, KS_eigenvalue, occ_numbers,max_occ_number,   &
             first_order_U,first_order_density_matrix)

    endif 

        call get_times(time_first_order_DM, clock_time_first_order_DM, &
        &              tot_time_first_order_DM, tot_clock_time_first_order_DM)

        !-------shanghui begin parallel------
     if (module_is_debugged("DFPT")) then
       if(myid.eq.0) then
       do i_spin = 1,n_spin
       write(use_unit,*) '************shanghui begain first_order_DM(X)****************'
       do i_basis=1,n_basis
       write(use_unit,'(7f15.9)') (first_order_density_matrix(1,i_basis,j_basis,i_spin), & 
                            j_basis=1,n_basis )
       enddo
       write(use_unit,*) '************shanghui end first_order_DM****************'
       write(use_unit,*) ''
       enddo 
       endif
     endif
        !-------shanghui end parallel--------


         !shanghui: we need to perform density matrix mixing here
    
          change_of_first_order_DM =0.0d0         

         do i_basis = 1,n_basis
         do j_basis = 1,n_basis

            do i_coord = 1, 3
            do i_spin = 1, n_spin

         change_of_first_order_DM =                &
         max( change_of_first_order_DM,             &
         dabs(first_order_density_matrix(i_coord,i_basis,j_basis, i_spin)  &
         - old_first_order_density_matrix(i_coord,i_basis,j_basis, i_spin)) )
   

          !if(number_of_loops.eq.1) then                  
          !  !-------begin linear mixing for benchmark-----------
          !  first_order_density_matrix(i_coord,i_basis,j_basis, i_spin) =       &
          !  (1.0d0-DFPT_mixing)*old_first_order_density_matrix(i_coord,i_basis,j_basis, i_spin)+  &
          !  DFPT_mixing*first_order_density_matrix(i_coord,i_basis,j_basis,i_spin)
          !  !-------end linear mixing for benchmark-----------
      
          !  old_first_order_density_matrix(i_coord,i_basis,j_basis,i_spin) =   &
          !  first_order_density_matrix(i_coord,i_basis,j_basis,i_spin)
          !endif

            enddo ! i_spin
            enddo ! i_coord
       
         enddo
         enddo

         converged = (change_of_first_order_DM.lt.DFPT_sc_accuracy_dm).and.(number_of_loops.gt.1) 

       if (.not.converged) then
 
          !if(number_of_loops > 0) then

           if(use_dfpt_pulay) then 

              if(.not.use_scalapack) then  
             ! Former way of calling pulay_mix:
              call pulay_mix(first_order_density_matrix, size(first_order_density_matrix), & 
                 number_of_loops, dfpt_pulay_steps, &
                 DFPT_mixing )
              else 
              ! Here matrix_tmp has local dimensions
                matrix_tmp(:,:,:,:) = first_order_density_matrix(:,my_row,my_col,:)
                call pulay_mix(matrix_tmp, size(matrix_tmp), number_of_loops, dfpt_pulay_steps, DFPT_mixing)
                first_order_density_matrix = 0d0
                first_order_density_matrix(:,my_row,my_col,:) = matrix_tmp(:,:,:,:)
                call sync_vector(first_order_density_matrix, size(first_order_density_matrix))
              endif
             old_first_order_density_matrix =  first_order_density_matrix
           else ! linear mixing

             first_order_density_matrix =       &
             (1.0d0-DFPT_mixing)*old_first_order_density_matrix+  &
             DFPT_mixing*first_order_density_matrix

             old_first_order_density_matrix = first_order_density_matrix
           endif ! mixing
 
          !endif  ! number_of_loops > 1 

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

        write(info_str,'(A)') ''
        call localorb_info(info_str, use_unit,'(A)', OL_norm  )

        write(info_str,'(2X,A)') &
        "CPSCF convergence accuracy:"
        call localorb_info(info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(2X,A,1X,E10.4,1X,E10.4)') &
                "| Change of first_order_density_matrix     :", change_of_first_order_DM
        call localorb_info(info_str, use_unit,'(A)', OL_norm  )


        call get_timestamps(time_first_order_density, clock_time_first_order_density)

        do i_spin = 1, n_spin
        call integrate_first_order_rho_polarizability(partition_tab, l_shell_max,  &
             KS_eigenvector(:,:,i_spin,1),first_order_density_matrix(:,:,:,i_spin), &
             first_order_rho(:,i_spin,:))
        enddo ! i_spin

        call get_times(time_first_order_density, clock_time_first_order_density, &
        &              tot_time_first_order_density, tot_clock_time_first_order_density)

!      call  integrate_polarizability &
!        (partition_tab,first_order_rho,polarizability)
!
!       !-------shanghui begin debug_mode------
!      if (module_is_debugged("DFPT")) then
!       write(info_str,'(A)') "polarizability:---> "
!       call localorb_info(info_str, use_unit,'(A)', OL_norm )
!       write(info_str,*)  polarizability(1:3)
!       call localorb_info(info_str, use_unit,'(A)', OL_norm )
!       write(info_str,'(A)') " "
!       call localorb_info(info_str, use_unit,'(A)', OL_norm )
!      endif 
!       !-------shanghui end debug_mode--------

!------------(1) end first-order-density update and mixing--------


!------------(2)begain to calculate first_order_H-----------------

        call get_timestamps(time_first_order_potential, clock_time_first_order_potential)

        do i_coord= 1, 3
           
           first_order_total_rho = 0.0d0
           do i_spin = 1, n_spin
              first_order_total_rho(1:n_full_points) = & 
              first_order_total_rho(1:n_full_points) + &
              first_order_rho(i_coord,i_spin,1:n_full_points)  
           enddo

        call update_hartree_potential_p2_shanghui &
            ( hartree_partition_tab,first_order_total_rho(1:n_full_points),& 
            delta_v_hartree_part_at_zero, &
            delta_v_hartree_deriv_l0_at_zero, &
            multipole_moments, multipole_radius_sq, &
            l_hartree_max_far_distance, &
            outer_potential_radius )

        !-------shanghui begin debug_mode------
       if (module_is_debugged("DFPT")) then
         write(info_str,'(A)') '{----shanghui in cpscf_solver.f:-----for first_order_hartree_potential:'
         call localorb_info(info_str, use_unit,'(A)', OL_norm )
         write(info_str,*) 'i_coord:',i_coord,forces_on
         call localorb_info(info_str, use_unit,'(A)', OL_norm )
         write(info_str,'(A)') '----shanghui in cpscf_solver.f:-----for first_order_hartree_potential}'
         call localorb_info(info_str, use_unit,'(A)', OL_norm )
       endif
        !-------shanghui end debug_mode--------

        call sum_up_whole_potential_p2_shanghui &
                ( delta_v_hartree_part_at_zero, &
                delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
                partition_tab, first_order_total_rho(1:n_full_points), &
                first_order_total_potential(i_coord,1:n_full_points),  & !<--------get first_order_DM_potential
                .false., multipole_radius_sq, &
                l_hartree_max_far_distance, &
                outer_potential_radius)

        enddo !i_coord

        call get_times(time_first_order_potential, clock_time_first_order_potential, &
        &              tot_time_first_order_potential, tot_clock_time_first_order_potential)



        call get_timestamps(time_first_order_H, clock_time_first_order_H)

        call  integrate_first_order_H_polarizability &
        (hartree_potential,first_order_total_potential, rho, rho_gradient,&
         first_order_rho, &
         partition_tab, l_shell_max, & 
         first_order_density_matrix, &
         first_order_H &
       )

        if( use_hartree_fock ) then
          call add_HF_to_first_order_H_polarizability &
                 (first_order_density_matrix, first_order_H)
        endif   

        call get_times(time_first_order_H, clock_time_first_order_H, &
        &              tot_time_first_order_H, tot_clock_time_first_order_H)



       !-------shanghui begin parallel------
     if (module_is_debugged("DFPT")) then
       if(myid.eq.0) then
       do i_spin = 1, n_spin
       write(use_unit,*) '************shanghui begain first_order_H(X)****************'
        do i_basis=1,n_basis
        write(use_unit,'(7f15.9)') (first_order_H(1,i_basis,j_basis,i_spin),j_basis=1,n_basis )
        enddo
       write(use_unit,*) '************shanghui end first_order_H****************'
       write(use_unit,*) ''
       enddo
       endif 
     endif
       !-------shanghui end parallel--------

! ------------(2) end to calculate first_order_H-----------------




!--------------(3) begin calculation of first_order_U-----------------
     if(.not.use_elsi_dm_cpscf) then 
        call get_timestamps(time_Sternheimer, clock_time_Sternheimer)

        if (use_scalapack) then
           ! i) first transform dense first-order H lapack to dense first-order H scalapack

           call construct_first_order_ham_polar_scalapack(first_order_H)
             
           !----for debug--------
           !  call  get_first_order_ham_polar_scalapack(first_order_H,1)

           ! ii) evaluate first-order U scalapack

           call evaluate_first_order_U_polar_scalapack(occ_numbers,KS_eigenvalue)

           ! iii) for debug 
           ! call get_first_order_U_polar_scalapack(first_order_U, 1)


        else ! lapack version
           call evaluate_first_order_U_polarizability(first_order_H,  &
             KS_eigenvector, KS_eigenvalue, occ_numbers,  max_occ_number, &
             first_order_U,first_order_E)
        end if

        call get_times(time_Sternheimer, clock_time_Sternheimer, &
        &              tot_time_Sternheimer, tot_clock_time_Sternheimer)
        !-------shanghui begin parallel------
!debug       if(myid.eq.0) then
!debug       do i_spin = 1, n_spin
!debug       write(use_unit,*) '************shanghui begain first_order_U(X)****************'
!debug        do i_basis=1,n_basis
!debug        write(use_unit,'(6f15.9)') (first_order_U(1,i_basis,j_basis,i_spin),j_basis=1,n_basis )
!debug        enddo
!debug       write(use_unit,*) '************shanghui end first_order_U****************'
!debug       enddo
!debug       write(use_unit,*) ''
!debug       endif 
        !-------shanghui end parallel--------

     endif 
!--------------(3) end solve first_order_U problem----------------

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


!       this is the end of the self-consistency cycle.

  end do SCF_LOOP
! << ------ end self consistent cycle--------

      total_number_of_loops = total_number_of_loops + number_of_loops


   call get_timestamps(time_Hessian, clock_time_Hessian)
   call  integrate_polarizability &
        (partition_tab,first_order_rho,polarizability)

        polar_mean = (polarizability(1,1) + polarizability(2,2) + polarizability(3,3))/3.0

        write(info_str,'(A)') 'DFPT for polarizability:--->'
        call localorb_info(info_str, use_unit,'(A)', OL_norm)
        do i_coord=1,3
         write(info_str,*) polarizability(i_coord,1:3)
         call localorb_info(info_str, use_unit,'(A)', OL_norm)
        enddo
        write(info_str,'(A,F11.6)') "The mean polarizability is ", polar_mean
        call localorb_info(info_str, use_unit,'(A)', OL_norm )
        write (info_str,'(A)') 'DFPT polarizability (Bohr^3)        xx        yy        zz        xy        xz        yz'
        call localorb_info(info_str, use_unit,'(A)', OL_norm)
        write (info_str,'(2X,A,1X,6F10.5)') &
        '| Polarizability:--->          ', polarizability(1,1), polarizability(2,2), &
        polarizability(3,3), polarizability(1,2), polarizability(1,3), polarizability(2,3)
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
     deallocate(first_order_U)
     deallocate(first_order_E)
     deallocate(matrix_tmp)
     if(use_dfpt_pulay) then
     call cleanup_pulay_mixing()
     endif
  !-------------------end deallocate------------------------------------



    end subroutine cpscf_solver_polarizability
!******
