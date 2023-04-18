!****s* FHI-aims/calculate_first_order_S_and_H_phonon_reduce_memory
!  NAME
!    calculate_first_order_S_and_H_phonon_reduce_memory
!  SYNOPSIS

    subroutine calculate_first_order_S_and_H_phonon_reduce_memory &
              (i_q_point, j_atom, j_coord, & 
               first_order_S_sparse, first_order_S_complex, &
               first_order_H_sparse, first_order_H_complex, & 
               first_order_density_matrix_sparse, first_order_U_complex, &
               converged)

!  PURPOSE
!  a module to get first_order_S and first_order_H for phonon_reduce_memory 
!  shanghui 2017 @Berlin.
!  
!  USES

      use runtime_choices
      use dimensions
      use timing
      use physics
      use species_data
      use localorb_io
      use geometry
      use scalapack_wrapper  , only : my_k_point, & 
                                      construct_first_order_overlap_scalapack, & 
                                      construct_first_order_hamiltonian_scalapack, & 
                                      evaluate_first_order_U_scalapack
                                      !get_first_order_overlap_sparse_matrix_scalapack, & just for debug
                                      !get_first_order_hamiltonian_sparse_matrix_scalapack & just for debug


      use pbc_lists , only : kweight_occs
      use mpi_tasks, only : aims_stop, myid, n_tasks
      use synchronize_mpi, only : sync_vector_complex

      use debugmanager, only: module_is_debugged

      !---------begin add for pulay_mixing--------
      use runtime_choices, only : use_dfpt_pulay, dfpt_pulay_steps 
      !this dfpt_pulay_steps by default is set to 8, if set in control.in to 1, then equals to linear mix. 
      use DFPT_pulay_mixing,    only: pulay_mix, cleanup_pulay_mixing
      !--------end add for pulay_mixing-----------


      implicit none

!  ARGUMENTS

     integer , intent(in) :: i_q_point
     integer , intent(in) :: j_atom
     integer , intent(in) :: j_coord
     complex*16, dimension(n_hamiltonian_matrix_size), intent(inout) :: first_order_S_sparse
     complex*16, dimension(n_basis, n_basis,n_k_points_task), intent(inout) :: first_order_S_complex
     complex*16, dimension(n_hamiltonian_matrix_size,n_spin), intent(inout) :: first_order_H_sparse
     complex*16, dimension(n_basis, n_basis,n_k_points_task), intent(inout) :: first_order_H_complex
     complex*16, dimension(n_basis, n_basis,n_k_points_task), intent(inout) :: first_order_U_complex
     complex*16, dimension(n_hamiltonian_matrix_size), intent(inout) ::  first_order_density_matrix_sparse

     
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


!-------------------shanghui begin define variables------------------------------------
!------------------------(1) grid ---------------------------------------
     complex*16, allocatable :: first_order_rho(:)
     complex*16, allocatable :: first_order_potential(:)

     real*8, allocatable :: first_order_rho_Re(:)
     real*8, allocatable :: first_order_rho_for_Born(:,:)
     real*8, allocatable :: first_order_rho_Im(:)
     real*8, allocatable :: first_order_potential_Re(:)
     real*8, allocatable :: first_order_potential_Im(:)

     complex*16, allocatable :: rho_free_gradient(:)
     complex*16, allocatable :: v_free_gradient(:)

!------------------------(2) matrix -----------------------------------
     real*8, allocatable ::  density_matrix_sparse(:)
     real*8, allocatable ::  first_order_density_matrix_sparse_real(:)
     complex*16, allocatable ::  old_first_order_density_matrix_sparse(:)


!------------------------(4) counters------------------------------
     integer :: i_basis, j_basis
     integer :: i_k_point,i_k_task
     integer :: i_point

     character*100 :: info_str
     character*8  :: cdate
     character*10 :: ctime
     real*8  time_start,time_end
     character(*), parameter :: deffmt = '2X'

     real*8  change_of_first_order_DM
     logical :: below_it_limit
     real*8, dimension(3, n_atoms)           :: hellman_feynman_term_at_gamma
!-------------------shanghui end define variables------------------------------------




!-------------------shanghui begin allocate variables------------------------------------
!----------------------(1)grid-------------------------------------------
       allocate(first_order_rho(n_full_points))
       allocate(first_order_potential(n_full_points))

       allocate(first_order_rho_Re(n_full_points))
       allocate(first_order_rho_for_Born(3,n_full_points))
       allocate(first_order_rho_Im(n_full_points))
       allocate(first_order_potential_Re(n_full_points))
       allocate(first_order_potential_Im(n_full_points))

       allocate(rho_free_gradient(n_full_points))
       allocate(v_free_gradient(n_full_points))

!----------------------(2)matrix-----------------------------------------


       allocate(density_matrix_sparse(n_hamiltonian_matrix_size))
       allocate(first_order_density_matrix_sparse_real(n_hamiltonian_matrix_size))
       allocate(old_first_order_density_matrix_sparse(n_hamiltonian_matrix_size))

!-------------------shanghui end allocate variables-----------------------------------


       if(packed_matrix_format.ne.PM_index) then
       call aims_stop('shanghui only use sparse matrix for DFPT_phonon_reduce_memory', & 
                      'cpscf_solver_phonon_reduce_memory')
       endif


!----------- becase DFPT is working for the second derivative of E, so we use j_atom,j_coord

      write(info_str,'(A)') ''
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
      write(info_str,'(A)') ''
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
      write(info_str,'(A)') "==========================================================================="
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
      write(info_str,'(2X,A,1X,I4,5X,A,1X,I4)') 'CPSCF working for i_q_point =',i_q_point
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      write(info_str,'(2X,A,1X,I4,5X,A,1X,I4)') 'CPSCF working for j_atom =',j_atom,'j_coord =',j_coord
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
  
    !-------shanghui begin parallel------
    !if(myid.eq.0) then
    ! write(use_unit,*) '-----------------------------------------------'
    ! write(use_unit,*) 'shanghui test n_basis:',n_basis
    ! write(use_unit,*) 'shanghui test n_hamiltonian_matrix_size',n_hamiltonian_matrix_size
    ! write(use_unit,*) '-----------------------------------------------'
    !endif
    !-------shanghui end parallel------



!--------------CPSCF (0.1) begin first_order_S-------------------
     call get_timestamps(time_first_order_S, clock_time_first_order_S)

     call integrate_first_order_S_phonon_reduce_memory(partition_tab, l_shell_max, & 
                                                       i_q_point, j_atom, j_coord, & 
                                                       first_order_S_sparse)

       !if(myid.eq.0) then
       !write(use_unit,*) '************first_order_S(atom1,X),before scalapck****************'
       !write(use_unit,*) (first_order_S_sparse(i_basis),i_basis=1,n_hamiltonian_matrix_size-1)
       !write(use_unit,*) '************first_order_S(atom1,X),before scalapck****************'
       !endif

     if(use_scalapack) then
       call construct_first_order_overlap_scalapack(first_order_S_sparse)
       !-------begin debug for scalapack------------ 
        !call get_first_order_overlap_sparse_matrix_scalapack(first_order_S_sparse)
        !if(myid.eq.0) then
        !write(use_unit,*) '************first_order_S(atom1,X),after scalapck****************'
        !write(use_unit,*) (first_order_S_sparse(i_basis),i_basis=1,n_hamiltonian_matrix_size-1)
        !write(use_unit,*) '************first_order_S(atom1,X),after scalapck****************'
        !endif
       !-------end debug for scalapck---------------- 
     else 
       i_k_task=0
       do i_k_point=1,n_k_points
         if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
            i_k_task = i_k_task + 1
            call construct_first_order_matrix_phonon_reduce_memory(first_order_S_sparse, &
                 first_order_S_complex(1:n_basis,1:n_basis,i_k_task), &
                 i_k_point)

         endif


       enddo  !i_k_point

     endif  !scalapck if.




     !-------to be delected soon--------
     !call construct_first_order_S_phonon_reduce_memory(first_order_S_sparse,first_order_S_complex)

     call get_times(time_first_order_S, clock_time_first_order_S, &
     &              tot_time_first_order_S, tot_clock_time_first_order_S)
!--------------CPSCF (0.1) end first_order_S---------------   
    


!--------------CPSCF (0.2) begin calculate gradient_free_V-------------------------------
    call  integrate_free_atom_sum_gradient_phonon_reduce_memory &
          (partition_tab, i_q_point, j_atom, j_coord, rho_free_gradient,v_free_gradient) 
!--------------CPSCF (0.2) end calculate gradient_free_V-------------------------------

!---------------CPSCF (0.3) shanghui test how to calcualte DM_sparse---------------------------------

    call  evaluate_zero_order_DM_phonon_reduce_memory &
          (KS_eigenvector, KS_eigenvector_complex, occ_numbers,density_matrix_sparse)
    !-------shanghui begin parallel------
    if (module_is_debugged("DFPT")) then
     if(myid.eq.0) then
      write(use_unit,*) '************zero_order_dm(atom1,X)****************'
      write(use_unit,'(40f20.15)') (density_matrix_sparse(i_basis),i_basis=1,n_hamiltonian_matrix_size-1)
      write(use_unit,*) '************zero_order_dm(atom1,X)****************'
     endif
    endif
    !-------shanghui end parallel------
!---------------CPSCF (0.3) shanghui end test how to calcualte DM_sparse-----------------------------



!--------------CPSCF (0.4) begin calculate init_first_order_DM1-------------------------------
    first_order_U_complex = (0.0d0, 0.0d0)
    call evaluate_first_order_DM_phonon_reduce_memory(first_order_S_complex,  &
         KS_eigenvector, KS_eigenvector_complex, occ_numbers,   &
         first_order_U_complex,old_first_order_density_matrix_sparse)
    !-------shanghui begin parallel------
    !if(myid.eq.0) then
    ! write(use_unit,*) '************old_first_order_dm(atom1,X)****************'
    ! write(use_unit,*) (old_first_order_density_matrix_sparse(i_basis),i_basis=1,n_hamiltonian_matrix_size-1)
    ! write(use_unit,*) '************old_first_order_dm(atom1,X)****************'
    !endif
    !-------shanghui end parallel------
!--------------CPSCF (0.4) end calculate init_first_order_DM1-------------------------------

! ------------------------ self-consistency loop -------------->>
  SCF_LOOP: do while ( (.not.converged) .and.  &
  &                    below_it_limit )
        number_of_loops = number_of_loops + 1


        write(info_str,'(A)') ''
        call localorb_info(info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "---------------------shanghui@@@CPSCF---------------------------------------"
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


!----------CPSCF: (1) Begin  update first_order_rho-----------------------------------
       call cpu_time(time_start)
       call get_timestamps(time_first_order_DM, clock_time_first_order_DM)

       call evaluate_first_order_DM_phonon_reduce_memory(first_order_S_complex,  &
            KS_eigenvector, KS_eigenvector_complex, occ_numbers,   &
            first_order_U_complex,first_order_density_matrix_sparse)

       call get_times(time_first_order_DM, clock_time_first_order_DM, &
        &              tot_time_first_order_DM, tot_clock_time_first_order_DM)


    !-------shanghui begin parallel------
    if (module_is_debugged("DFPT")) then
     if(myid.eq.0) then
      write(use_unit,*) '************first_order_dm(atom1,X)****************'
      write(use_unit,'(40f20.15)') (first_order_density_matrix_sparse(i_basis),i_basis=1,n_hamiltonian_matrix_size-1)
      write(use_unit,*) '************first_order_dm(atom1,X)****************'
     endif
    endif 
   !-------shanghui end parallel------

         change_of_first_order_DM =0.0d0         

         do i_basis = 1, n_hamiltonian_matrix_size - 1 
            change_of_first_order_DM =                &
            max( change_of_first_order_DM,             &
            dabs( dble(first_order_density_matrix_sparse(i_basis)  &
               - old_first_order_density_matrix_sparse(i_basis))) )
         enddo
     

         if(use_dfpt_pulay) then  
         !---------begin add for pulay_mixing--------
            ! As pulay_mix only works for real type, I first translate the complex one to real type.
            do i_basis = 1, n_hamiltonian_matrix_size  
               first_order_density_matrix_sparse_real(i_basis) = & 
               dble(first_order_density_matrix_sparse(i_basis))
            enddo 

            if(number_of_loops.eq.1) then
            ! shanghui add note: for loop=1, we need to inital X_in to  
            ! first_order_density_matrix_sparse, otherwise, the X_in 
            ! will set to 0.0d0, this is fine for electric_field response, 
            ! however, in phonon, X_in is not 0.0d0, we need to set here.   
            call pulay_mix(first_order_density_matrix_sparse_real, n_hamiltonian_matrix_size, &
                 number_of_loops, dfpt_pulay_steps, &
                 1.0d0 ) 
            else ! > 1
            call pulay_mix(first_order_density_matrix_sparse_real, n_hamiltonian_matrix_size, &
                 number_of_loops-1, dfpt_pulay_steps, &
                 DFPT_mixing )
            endif
 
            ! Then change back to complex type  
            do i_basis = 1, n_hamiltonian_matrix_size
              old_first_order_density_matrix_sparse(i_basis) = & 
              dcmplx(first_order_density_matrix_sparse_real(i_basis), 0.0d0)
              first_order_density_matrix_sparse(i_basis) = & 
              dcmplx(first_order_density_matrix_sparse_real(i_basis), 0.0d0)
            enddo
         !---------end add for pulay_mixing--------
         else ! linear mixing 
            do i_basis = 1, n_hamiltonian_matrix_size
              first_order_density_matrix_sparse(i_basis) =       &
              (1.0d0-DFPT_mixing)*old_first_order_density_matrix_sparse(i_basis)+  &
              DFPT_mixing*first_order_density_matrix_sparse(i_basis)
          
              old_first_order_density_matrix_sparse(i_basis) =   &
              first_order_density_matrix_sparse(i_basis)
            enddo
         endif ! mixing 


        !-------shanghui begin parallel------
        if(myid.eq.0) then
        write(use_unit,*) "((((((((((((((((((((((((((((((((("
        write(use_unit,*) change_of_first_order_DM
        write(use_unit,*) ")))))))))))))))))))))))))))))))))"
        endif 
        !-------shanghui end parallel--------

!--------CPSCF : (1) end first-order-density update and mixing--------

!--------CPSCF : (2) begain to calculate first_order_H-----------------
        !----in pbc case, we have used density matrix instead of orbitals, it looks neat. 
        call get_timestamps(time_first_order_density, clock_time_first_order_density)

        call integrate_first_order_rho_phonon_reduce_memory(partition_tab, l_shell_max,  &
             density_matrix_sparse,first_order_density_matrix_sparse, &
             i_q_point, j_atom, j_coord, &
             first_order_rho)

             first_order_rho(1:n_full_points)=  &           !@ 
             first_order_rho(1:n_full_points)+  &           !@
             rho_free_gradient(1:n_full_points)             !@

          do i_point =1 ,n_full_points
             first_order_rho_Re(i_point)=  &           
             dble(first_order_rho(i_point)) 

             first_order_rho_Im(i_point)=  &           
             dimag(first_order_rho(i_point))  

          enddo  


        call get_times(time_first_order_density, clock_time_first_order_density, &
        &              tot_time_first_order_density, tot_clock_time_first_order_density)


        call get_timestamps(time_first_order_potential, clock_time_first_order_potential)

        call update_hartree_potential_shanghui_phonon_reduce_memory &
            (hartree_partition_tab,first_order_rho_Re(1:n_full_points),& 
             delta_v_hartree_part_at_zero, &
             delta_v_hartree_deriv_l0_at_zero, &
             multipole_moments, multipole_radius_sq, &
             l_hartree_max_far_distance, &
             outer_potential_radius )

        call sum_up_whole_potential_shanghui_phonon_reduce_memory &
            (delta_v_hartree_part_at_zero, &
             delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
             partition_tab, first_order_rho_Re(1:n_full_points), &
             first_order_potential_Re(1:n_full_points),  & 
             .false., multipole_radius_sq, &
             l_hartree_max_far_distance, &
             outer_potential_radius, & 
             hellman_feynman_term_at_gamma) 

        if(.false.) then ! we only turned on when we need q.ne.0  
        call update_hartree_potential_shanghui_phonon_reduce_memory &
            (hartree_partition_tab,first_order_rho_Im(1:n_full_points),& 
             delta_v_hartree_part_at_zero, &
             delta_v_hartree_deriv_l0_at_zero, &
             multipole_moments, multipole_radius_sq, &
             l_hartree_max_far_distance, &
             outer_potential_radius )

        call sum_up_whole_potential_shanghui_phonon_reduce_memory &
            (delta_v_hartree_part_at_zero, &
             delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
             partition_tab, first_order_rho_Im(1:n_full_points), &
             first_order_potential_Im(1:n_full_points),  & 
             .false., multipole_radius_sq, &
             l_hartree_max_far_distance, &
             outer_potential_radius, & 
             hellman_feynman_term_at_gamma)
         else
            first_order_potential_Im = 0.0d0
         endif  


          do i_point =1 ,n_full_points
             first_order_potential(i_point)=  &           
             dcmplx(first_order_potential_Re(i_point),first_order_potential_Im(i_point)) 
          enddo  
 
               first_order_potential(1:n_full_points)= &      !@
               first_order_potential(1:n_full_points)   &     !@
               +(-v_free_gradient(1:n_full_points))           !@

        call get_times(time_first_order_potential, clock_time_first_order_potential, &
        &              tot_time_first_order_potential, tot_clock_time_first_order_potential)


        call get_timestamps(time_first_order_H, clock_time_first_order_H)

        !call integrate_first_order_rho_phonon_reduce_memory(partition_tab, l_shell_max,  &
        !     density_matrix_sparse,first_order_density_matrix_sparse, &
        !     i_q_point, j_atom, j_coord, &
        !     first_order_rho)
             first_order_rho(1:n_full_points)=  &           !@ 
             first_order_rho(1:n_full_points)-  &           !@
             rho_free_gradient(1:n_full_points)             !@

        call  integrate_first_order_H_phonon_reduce_memory &
            (hartree_potential,first_order_potential, & 
             rho, rho_gradient, first_order_rho, &
             partition_tab, l_shell_max,    &
             first_order_density_matrix_sparse, density_matrix_sparse, &
             i_q_point, j_atom, j_coord,    &
             first_order_H_sparse(:,1))
!   if(myid.eq.0) then
!    write(use_unit,*) '************first_order_H(atom1,X)****************'
!    write(use_unit,'(40f20.15)') (first_order_H_sparse(i_basis,1),i_basis=1,n_hamiltonian_matrix_size-1)
!    write(use_unit,*) '************first_order_H(atom1,X)****************'
!   endif

        call get_times(time_first_order_H, clock_time_first_order_H, &
         &              tot_time_first_order_H, tot_clock_time_first_order_H)
        
!--------CPSCF : (2) end to calculate first_order_H-----------------


!--------CPSCF : (3) begain to calculate first_order_U-----------------
       call get_timestamps(time_Sternheimer, clock_time_Sternheimer)

      !if(myid.eq.0) then
      !write(use_unit,*) '************first_order_H(atom1,X),before scalapck****************'
      !write(use_unit,*) (first_order_H_sparse(i_basis,1),i_basis=1,n_hamiltonian_matrix_size-1)
      !write(use_unit,*) '************first_order_H(atom1,X),before scalapck****************'
      !endif


    if(use_scalapack)then

      call  construct_first_order_hamiltonian_scalapack(first_order_H_sparse)  
       !-------begin debug for scalapack------------ 
       !call get_first_order_hamiltonian_sparse_matrix_scalapack(first_order_H_sparse(:,1))
       ! if(myid.eq.0) then
       ! write(use_unit,*) '************first_order_H(atom1,X),after scalapck****************'
       ! write(use_unit,*) (first_order_H_sparse(i_basis,1),i_basis=1,n_hamiltonian_matrix_size-1)
       ! write(use_unit,*) '************first_order_H(atom1,X),after scalapck****************'
       ! endif
       !-------end debug for scalapck---------------- 



      call  evaluate_first_order_U_scalapack(occ_numbers, KS_eigenvalue)
 
    else ! lapack version 
      i_k_task = 0
      do i_k_point = 1, n_k_points, 1
      if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
         i_k_task = i_k_task + 1   
      call construct_first_order_matrix_phonon_reduce_memory(first_order_H_sparse(:,1), &
           first_order_H_complex(1:n_basis,1:n_basis,i_k_task), i_k_point)
      call evaluate_first_order_U_phonon_reduce_memory(& 
           first_order_H_complex(:,:,i_k_task),& 
           first_order_S_complex(:,:,i_k_task),  &
           KS_eigenvector_complex(:,:,:,i_k_task), KS_eigenvalue(:,:,i_k_point), &  
           occ_numbers(:,:,i_k_point),  &
           first_order_U_complex(:,:,i_k_task))
      endif
      enddo ! n_k_point

    endif
 


        call cpu_time(time_end)
        if(myid.eq.0) then 
        write(use_unit,*) '@@@@@@@@time for CPSCF@@@@@@@@@@@@@@@@@@',time_end-time_start
        endif

        call get_times(time_Sternheimer, clock_time_Sternheimer, &
        &              tot_time_Sternheimer, tot_clock_time_Sternheimer)
!-------CPSCF: (3) end solve first_order_U problem----------------


! --------- Check convergence ----->>


!         check convergence of self-consistency loop
         

        converged = (change_of_first_order_DM.lt.DFPT_sc_accuracy_dm).and.(number_of_loops.ne.1)

!  ---------- Update electron density and perform mixing ---->>

        if (converged) then
!           We are done - no further evaluation of density / potential needed

          write(info_str,'(A)') ''
          call localorb_info(info_str, use_unit,'(A)',OL_norm)
          write(info_str,'(2X,A)') "CP-self-consistency cycle converged."
          call localorb_info(info_str, use_unit,'(A)',OL_norm)
          write(info_str,'(A)') ''
          call localorb_info(info_str, use_unit,'(A)',OL_norm)
 

        else if (number_of_loops.ge.sc_iter_limit) then
!           This was the last self-consistency cycle - we do not need any more potential / density evaluations

          below_it_limit = .false.
        end if

        call get_times(time_cpscf_loop, clock_time_cpscf_loop)

        ! current CPSCF loop ends here

! ----- Printing out time data -- >>
        write(info_str,'(A,I5)') &
        & "End CPSCF iteration # ", number_of_loops
        call output_timeheader(deffmt, info_str, OL_norm)
        call output_times(deffmt, "Time for this iteration", &
        &                 time_cpscf_loop, clock_time_cpscf_loop, OL_norm)
        call output_times(deffmt, "first_order_DM", &
        &                 time_first_order_DM, clock_time_first_order_DM, OL_norm)
        call output_times(deffmt, "first_order_density", &
          &                 time_first_order_density, clock_time_first_order_density, OL_norm)
        call output_times(deffmt, "first_order_potential", &
        &                 time_first_order_potential, clock_time_first_order_potential, OL_norm)
        call output_times(deffmt, "first_order_H", &
        &                 time_first_order_H, clock_time_first_order_H, OL_norm)
        call output_times(deffmt, "Solution of Sternheimer eqns.", &
        &                 time_Sternheimer, clock_time_Sternheimer, OL_norm)
        write(info_str,'(A)') &
        "------------------------------------------------------------"
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

! << ---- end printing out data---------




  end do SCF_LOOP
! << ------ end self consistent cycle--------

     if(use_dfpt_pulay) then  
    !------begin add for pulay_mixing-----------
     call cleanup_pulay_mixing()
    !------end add for pulay_mixing-----------
     endif



!-------------------shanghui begin deallocate------------------------------------
!----------------------(1) grid------------------------------------
       deallocate(first_order_rho)
       deallocate(first_order_potential)

       deallocate(first_order_rho_Re)
       deallocate(first_order_rho_Im)
       deallocate(first_order_potential_Re)
       deallocate(first_order_potential_Im)

       deallocate(rho_free_gradient)
       deallocate(v_free_gradient)

!----------------------(2) matrix-----------------------------------------
       deallocate(density_matrix_sparse)
       deallocate(old_first_order_density_matrix_sparse)
!-------------------shanghui end deallocate------------------------------------



    end subroutine calculate_first_order_S_and_H_phonon_reduce_memory
