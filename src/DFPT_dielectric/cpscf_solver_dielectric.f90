!****s* FHI-aims/cpscf_solver_dielectric
!  NAME
!    cpscf_solver_dielectric
!  SYNOPSIS

    subroutine cpscf_solver_dielectric &
    (converged)

!  PURPOSE
!  an cpscf process for dielectric 
!  shanghui 2015.07
!  
!  USES

      use constants, only: pi, pi4_inv
      use runtime_choices
      use dimensions
      use timing
      use physics  ! overlap_matrix, hamiltonian
      use species_data
      use localorb_io
      use geometry
      use runtime_choices, only : DFPT_sc_accuracy_dm, DFPT_mixing
      !---------begin add for pulay_mixing--------
      use DFPT_pulay_mixing, only: pulay_mix, cleanup_pulay_mixing
      !--------end add for pulay_mixing-----------
      use scalapack_wrapper  , only : construct_momentum_matrix_dielectric_scalapack, &
                                      construct_first_order_hamiltonian_dielectric_scalapack, &
                                      evaluate_first_order_U_dielectric_scalapack, & 
                                      first_order_U_scalapack, first_order_U_complex_scalapack,&
                                      first_order_ham_scalapack, first_order_ham_complex_scalapack, &
                                      construct_hamiltonian_real_for_elsi_scalapack, &   
                                      construct_overlap_real_for_elsi_scalapack,     & 
                                      set_full_matrix_real, &
                                      mxld, mxcol, &
                                      construct_momentum_matrix_dielectric_for_elsi_scalapack, &
                                      construct_first_order_hamiltonian_dielectric_for_elsi_scalapack, &
                                      get_first_order_dm_sparse_matrix_dielectric_for_elsi_scalapack
      use pbc_lists , only : kweight_occs
      use synchronize_mpi, only : sync_vector_complex 
      use mpi_tasks, only: aims_stop, myid, n_tasks

      use aims_memory_tracking, only: aims_allocate, aims_deallocate 
      use elsi, only: elsi_dm_real, elsi_dm_real_cpscf
      use elsi_wrapper, only: eh_scf
      use synchronize_mpi_basic, only: sync_vector

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


!-------------------shanghui begin define variables------------------------------------
!------------------------(1) grid ---------------------------------------
     real*8, allocatable :: first_order_rho(:)
     real*8, allocatable :: first_order_potential(:)


!------------------------(2) matrix -----------------------------------
!   now is mpi version.
!   for scalapck version in the future, we just need smaller (la,lb) for first_order_S_complex
     real*8, allocatable     :: first_order_H_sparse(:)
     complex*16, allocatable :: first_order_H_complex(:,:,:) 
     complex*16, allocatable :: first_order_U_complex(:,:,:,:)
     
     complex*16, allocatable :: Omega_MO(:,:,:,:)    ! = < i(k)|-r| j(k) > 
    !<1> = (C^+ momentum_matrix_complex C)_ij/( Ei(k) - Ej(k) ) 
     real*8, allocatable     :: momentum_matrix_sparse(:)
     complex*16, allocatable :: momentum_matrix_complex(:,:,:,:)



     real*8, allocatable ::  first_order_density_matrix_sparse(:)
     real*8, allocatable ::  old_first_order_density_matrix_sparse(:)

     real*8  dielectric_constant(3,3)
     real*8  Born_effective_charges(3,n_atoms,3)
     real*8  Born_effective_charges_HF(3,n_atoms)
     real*8  Born_effective_charges_MP(3,n_atoms)
     real*8  Born_effective_charges_Pulay(3,n_atoms)

     complex*16 dP_dE(3,3)

!------------------------(3) counters------------------------------
     integer :: i_atom, i_coord, j_coord 
     integer :: i_basis, j_basis
     integer :: i_k_point,i_k_task
     integer :: i_point
     integer :: i_spin

     !----------------(4) for elsi_ntpoly, the scalapack matrix------------------------
     real*8, dimension(:,:,:), allocatable   :: ham
     real*8, dimension(:,:),   allocatable   :: ovlp
     real*8, dimension(:,:,:), allocatable   :: dm
     real(8)  :: energy
     real*8, dimension(:,:,:), allocatable   :: mm_scalapack
     real*8, dimension(:,:,:), allocatable   :: first_order_ham
     real*8, dimension(:,:),   allocatable   :: first_order_ovlp
     real*8, dimension(:,:,:), allocatable   :: first_order_dm


     character*1000 :: info_str
     character*8  :: cdate
     character*10 :: ctime
      character(*), parameter :: deffmt = '2X'
     real*8  time_start,time_end

     real*8  change_of_first_order_DM
     logical :: below_it_limit
     logical, parameter :: calculate_Born_effective_charges = .false.
     logical, parameter :: test_Born_effective_charges_by_print_force = .false.
     logical, parameter :: test_Born_effective_charges_at_every_coord = .false.
     logical, parameter ::  use_elsi_dm_cpscf = .false.

!-------------------shanghui end define variables------------------------------------


  if(myid.eq.0) then 
  write(use_unit,*) "-------------------------------------------------------------------------------"
  write(use_unit,*) "|              ENTERING DFPT_DIELECTRIC (PERIODIC CALCULATION)                |"
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

!-------------------shanghui begin allocate variables------------------------------------
!----------------------(1)grid-------------------------------------------
       call aims_allocate(first_order_rho,n_full_points, "+first_order_rho")
       call aims_allocate(first_order_potential,n_full_points, "+first_order_potential")
       

!----------------------(2)matrix-----------------------------------------
       call aims_allocate(first_order_H_sparse, n_hamiltonian_matrix_size,"+first_order_H_sparse")
       call aims_allocate(momentum_matrix_sparse, n_hamiltonian_matrix_size_no_symmetry,"+momentum_matrix_sparse")
       call aims_allocate(first_order_density_matrix_sparse, n_hamiltonian_matrix_size,"+first_order_density_matrix_sparse")
       call aims_allocate(old_first_order_density_matrix_sparse, n_hamiltonian_matrix_size,"+old_first_order_density_matrix_sparse")

      if(use_scalapack) then
       allocate(first_order_U_complex(1,1,1,3))
      else 
       allocate(first_order_U_complex(n_states,n_states,n_k_points_task,3))
      endif


      call aims_allocate(mm_scalapack,  mxld, mxcol, 3, "+mm_scalpack")

!-------------------shanghui end allocate variables-----------------------------------



       if(packed_matrix_format.ne.PM_index) then
       call aims_stop('shanghui only use sparse matrix for DFPT_dielectric', & 
                      'cpscf_solver_dielectric')
       endif

       if(use_scalapack) then
       !  in lapack version we will add k_weights by ourself, so only scalapck version need this.
       call kweight_occs('cpscf_solver_dielectric', occ_numbers)
       endif


!---------------CPSCF (0.1) shanghui calculate momentum_matrix---------------------------------
      write(info_str,'(A)') ''
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
      write(info_str,'(A)') ''
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
      write(info_str,'(A)') "==========================================================================="
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
      write(info_str,'(2X,A,1X,I4)') 'before CPSCF, get Omega_MO = <i(k)|-r|j(k)> in MO basis'
      call localorb_info(info_str, use_unit,'(A)',OL_norm)
      write(info_str,'(A)') ''
      call localorb_info(info_str, use_unit,'(A)', OL_norm )
      write(info_str,'(A)') "=========================================================================="
      call localorb_info(info_str, use_unit,'(A)', OL_norm )

    if(.not.use_scalapack)then !lapack version

       allocate(first_order_H_complex(n_basis, n_basis,n_k_points_task))
       allocate(Omega_MO(n_states,n_states,n_k_points_task,3))
       allocate(momentum_matrix_complex(n_basis,n_basis,n_k_points_task,3))

      !---------------<1> begin from momentum_matrix_complex---------------
      momentum_matrix_sparse = 0.0d0
      momentum_matrix_complex = (0.0d0,0.0d0)
      Omega_MO = (0.0d0, 0.0d0)
 
      do j_coord = 1,3

        call  integrate_momentum_matrix_sparse &
             ( partition_tab, l_shell_max,    &
               j_coord, momentum_matrix_sparse)
 
         i_k_task = 0
         do i_k_point = 1,n_k_points, 1
         if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
      
            i_k_task = i_k_task + 1   
            call construct_matrix_complex_no_symmetry(momentum_matrix_sparse, &
            momentum_matrix_complex(1:n_basis,1:n_basis,i_k_task,j_coord),i_k_point)
 
         endif ! i_k_task
         enddo ! n_k_point
 
        call evaluate_Omega_MO_v1(&
             momentum_matrix_complex(:,:,:,j_coord),&
             Omega_MO(:,:,:,j_coord))

      enddo
 
!      !---------------<1> end from momentum_matrix_complex---------------


!---------------CPSCF (0.1) shanghui end calcualte Omega_MO------------------------------
   else ! scalapack
     momentum_matrix_sparse = 0.0d0
     do j_coord = 1, 3
      call  integrate_momentum_matrix_sparse &
           ( partition_tab, l_shell_max,    &
             j_coord, momentum_matrix_sparse)
       if(.not.use_elsi_dm_cpscf) then
        call construct_momentum_matrix_dielectric_scalapack(momentum_matrix_sparse,j_coord)
       else 
        call construct_momentum_matrix_dielectric_for_elsi_scalapack(momentum_matrix_sparse, mm_scalapack(:,:,j_coord))
       endif
     enddo
   endif 


      !initial  

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

      if(use_scalapack) then
        if(real_eigenvectors) then  
        first_order_U_scalapack = 0.0d0
        else 
        first_order_U_complex_scalapack = (0.0d0,0.0d0)
        endif
      else
        first_order_U_complex(:,:,:,j_coord) = (0.0d0, 0.0d0)
      endif

    !-------shanghui begin parallel------
    if(myid.eq.0) then
     write(use_unit,*) '--------------------------------------------------------------------'
     write(use_unit,*) 'shanghui test n_basis:',n_basis
     write(use_unit,*) 'shanghui test n_states:',n_states
     write(use_unit,*) 'shanghui test n_hamiltonian_matrix_size',n_hamiltonian_matrix_size
     write(use_unit,*) '--------------------------------------------------------------------'
    endif
    !-------shanghui end parallel------
 

!--------------CPSCF (0.3) begin calculate init_first_order_DM1-------------------------------
!    call evaluate_first_order_DM_dielectric(  &
!         KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers,   &
!         first_order_U_complex(1:n_states,1:n_states,1:n_k_points_task,j_coord), & 
!         old_first_order_density_matrix_sparse)
!--------------CPSCF (0.3) end calculate init_first_order_DM1-------------------------------

    old_first_order_density_matrix_sparse=0.d0
    first_order_H_sparse = 0.d0

! Start with a nonzero guess for U1 (at this point it will only contain the contribution from Omega_MO)
    if(use_scalapack)then

      !call  construct_first_order_hamiltonian_dielectric_scalapack(first_order_H_sparse)

      if(real_eigenvectors) then
        first_order_ham_scalapack=0.d0
      else
        first_order_ham_complex_scalapack=(0.d0,0.d0)
      endif

      call  evaluate_first_order_U_dielectric_scalapack(occ_numbers, KS_eigenvalue,dP_dE,j_coord)

    else ! lapack

     first_order_H_complex = (0.d0,0.d0)

     call evaluate_first_order_U_dielectric(& 
          Omega_MO(:,:,:,j_coord),first_order_H_complex(:,:,:),& 
          first_order_U_complex(:,:,:,j_coord)) 

    endif ! lapack/scalapack

! ------------------------ self-consistency loop -------------->>
  SCF_LOOP: do while ( (.not.converged) .and.  &
  &                    below_it_limit )
        number_of_loops = number_of_loops + 1

        write(info_str,'(A)') ''
        call localorb_info(info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "--------------------- CPSCF ---------------------------------------"
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
       call get_timestamps(time_first_order_DM, clock_time_first_order_DM)


!---------#ELSI----------------------------
    if(use_elsi_dm_cpscf) then 

     if(use_scalapack.and.real_eigenvectors) then ! we only have scalapack version
 
      call aims_allocate(ham,  mxld, mxcol, n_spin, "+ham")
      call aims_allocate(ovlp, mxld, mxcol, "+ovlp")
      call aims_allocate(dm,  mxld, mxcol, n_spin, "+dm")
      call aims_allocate(first_order_ham,  mxld, mxcol, n_spin, "+first_order_ham")
      call aims_allocate(first_order_ovlp, mxld, mxcol, "+first_order_ovlp")
      call aims_allocate(first_order_dm,  mxld, mxcol, n_spin, "+first_order_dm")
 
         call construct_overlap_real_for_elsi_scalapack(overlap_matrix, ovlp)
         call construct_hamiltonian_real_for_elsi_scalapack(hamiltonian, ham)
 
         
          !-------1. global matrix to scalapack matrix--------
            first_order_ovlp = 0.0d0
         do i_spin = 1, n_spin   
            call construct_first_order_hamiltonian_dielectric_for_elsi_scalapack(first_order_H_sparse,first_order_ham(:,:,i_spin))
            call set_full_matrix_real(ovlp)
            call set_full_matrix_real(ham(:,:,i_spin))
            call set_full_matrix_real(first_order_ham(:,:,i_spin))
           !------------this is wrong, so we remove----------
           !do i_basis = 1, mxld
           !   do j_basis =1, mxcol
           !      first_order_ham(i_basis,j_basis,i_spin) = first_order_ham(i_basis,j_basis,i_spin) + mm_scalapack(i_basis,j_basis, j_coord) 
           !   enddo
           !enddo
          !-------2. get first_order_DM----------
            call elsi_dm_real_cpscf(eh_scf,ham(:,:,i_spin),ovlp,dm(:,:,i_spin), energy, & 
               first_order_ham(:,:,i_spin), first_order_ovlp, first_order_dm(:,:,i_spin)  )
          !-------3. scalapack matrix to global matrix---------
            call get_first_order_dm_sparse_matrix_dielectric_for_elsi_scalapack( first_order_dm(:,:,i_spin), first_order_density_matrix_sparse)
         enddo ! i_spin


         call sync_vector(first_order_density_matrix_sparse,n_hamiltonian_matrix_size)
 
      if(allocated(ham))  call aims_deallocate(ham, "+ham")
      if(allocated(ovlp)) call aims_deallocate(ovlp, "+ovlp")
      if(allocated(dm))   call aims_deallocate(dm, "+dm")
      if(allocated(first_order_ham))  call aims_deallocate(first_order_ham, "+first_order_ham")
      if(allocated(first_order_ovlp)) call aims_deallocate(first_order_ovlp, "+first_order_ovlp")
      if(allocated(first_order_dm))   call aims_deallocate(first_order_dm, "+first_order_dm")
 

     endif    ! use_scalapack and real



    else ! not use_elsi_dm_cpscf

       call evaluate_first_order_DM_dielectric( &
            KS_eigenvector,KS_eigenvector_complex, KS_eigenvalue, occ_numbers,   &
            first_order_U_complex(:,:,1:n_k_points_task,j_coord), & 
            first_order_density_matrix_sparse)
    endif 

        call get_times(time_first_order_DM, clock_time_first_order_DM, &
        &              tot_time_first_order_DM, tot_clock_time_first_order_DM)


        change_of_first_order_DM =0.0d0         

        do i_basis = 1, n_hamiltonian_matrix_size - 1 

         change_of_first_order_DM =                &
         max( change_of_first_order_DM,             &
         dabs( dble(first_order_density_matrix_sparse(i_basis)  &
            - old_first_order_density_matrix_sparse(i_basis))) )
         

         !if(number_of_loops.eq.1) then                  
         ! first_order_density_matrix_sparse(i_basis) =       &
         ! (1.0d0-DFPT_mixing)*old_first_order_density_matrix_sparse(i_basis)+  &
         ! DFPT_mixing*first_order_density_matrix_sparse(i_basis)
     
         ! old_first_order_density_matrix_sparse(i_basis) =   &
         ! first_order_density_matrix_sparse(i_basis)
         !endif

        enddo

        if(number_of_loops > 0) then
         if(use_dfpt_pulay) then  
         !---------begin add for pulay_mixing--------
            call pulay_mix(first_order_density_matrix_sparse, n_hamiltonian_matrix_size, & 
                 number_of_loops-0, dfpt_pulay_steps, &
                 DFPT_mixing )
         old_first_order_density_matrix_sparse = first_order_density_matrix_sparse
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
        endif  ! number_of_loops > 1 


        !-------shanghui begin parallel------
        if(myid.eq.0) then
        write(use_unit,*) "((((((((((((((((((((((((((((((((("
        write(use_unit,*) change_of_first_order_DM
        write(use_unit,*) ")))))))))))))))))))))))))))))))))"
        endif 
        !-------shanghui end parallel--------

        write(info_str,'(A)') ''
        call localorb_info(info_str, use_unit,'(A)', OL_norm  )

        write(info_str,'(2X,A)') &
        "CPSCF convergence accuracy:"
        call localorb_info(info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(2X,A,1X,E10.4,1X,E10.4)') &
                "| Change of first_order_density_matrix     :", change_of_first_order_DM
        call localorb_info(info_str, use_unit,'(A)', OL_norm  )


    !---------shanghui test compile@thnec (1)---------------
    ! if(myid.eq.0) then
    !  write(use_unit,*) '************first_order_dm(atom1,X)****************'
    !  write(use_unit,'(40f20.15)') (first_order_density_matrix_sparse(i_basis),i_basis=1,n_hamiltonian_matrix_size-1)
    !  write(use_unit,*) '************first_order_dm(atom1,X)****************'
    ! endif

!--------CPSCF : (1) end first-order-density update and mixing--------

!--------CPSCF : (2) begain to calculate first_order_H-----------------
        !----in pbc case, we have used density matrix instead of orbitals, it looks neat. 
        call get_timestamps(time_first_order_density, clock_time_first_order_density)

        call integrate_first_order_rho_dielectric(partition_tab, l_shell_max,  &
             first_order_density_matrix_sparse, &
             first_order_rho)

        call get_times(time_first_order_density, clock_time_first_order_density, &
        &              tot_time_first_order_density, tot_clock_time_first_order_density)


        call get_timestamps(time_first_order_potential, clock_time_first_order_potential)

        call update_hartree_potential_shanghui_dielectric &
            (hartree_partition_tab,first_order_rho(1:n_full_points),& 
             delta_v_hartree_part_at_zero, &
             delta_v_hartree_deriv_l0_at_zero, &
             multipole_moments, multipole_radius_sq, &
             l_hartree_max_far_distance, &
             outer_potential_radius )

        call sum_up_whole_potential_shanghui_dielectric &
            (delta_v_hartree_part_at_zero, &
             delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
             partition_tab, first_order_rho(1:n_full_points), &
             first_order_potential(1:n_full_points),  & 
              multipole_radius_sq, &
             l_hartree_max_far_distance, &
             outer_potential_radius)

        call get_times(time_first_order_potential, clock_time_first_order_potential, &
        &              tot_time_first_order_potential, tot_clock_time_first_order_potential)


        call get_timestamps(time_first_order_H, clock_time_first_order_H)

        call  integrate_first_order_H_dielectric &
             (hartree_potential,first_order_potential, & 
             rho, rho_gradient, first_order_rho, &
             partition_tab, l_shell_max,    &
             j_coord,    &
             first_order_density_matrix_sparse, &
             first_order_H_sparse)

        call get_times(time_first_order_H, clock_time_first_order_H, &
        &              tot_time_first_order_H, tot_clock_time_first_order_H)
!--------CPSCF : (2) end to calculate first_order_H-----------------


!--------CPSCF : (3) begin to calculate first_order_U-----------------
   
    if(.not.use_elsi_dm_cpscf) then 
        call get_timestamps(time_Sternheimer, clock_time_Sternheimer)
     if(use_scalapack)then

      call  construct_first_order_hamiltonian_dielectric_scalapack(first_order_H_sparse)
      call  evaluate_first_order_U_dielectric_scalapack(occ_numbers, KS_eigenvalue,dP_dE,j_coord)
      !call  sync_vector_complex(dP_dE, 3*3)
  
      !write(info_str,'(A)') 'Warning: only dP_dE is done at every cycle:--->'
      !call localorb_info(info_str, use_unit,'(A)', OL_norm)
      !write(info_str,*) dP_dE
      !call localorb_info(info_str, use_unit,'(A)', OL_norm)

      !do i_coord=1,3
      ! write(use_unit,*) 'dP_TEST', dP_dE(1,i_coord),dP_dE(2,i_coord),dP_dE(3,i_coord)
      !enddo

      !write(info_str,'(A)') 'dP_dE (Bohr^3) at every cycle:--->'
      !call localorb_info(info_str, use_unit,'(A)', OL_norm)
      !write(info_str,*) dP_dE(1,1:3)
      !call localorb_info(info_str, use_unit,'(A)', OL_norm)
      !write(info_str,*) dP_dE(2,1:3)
      !call localorb_info(info_str, use_unit,'(A)', OL_norm)
      !write(info_str,*) dP_dE(3,1:3)
      !call localorb_info(info_str, use_unit,'(A)', OL_norm)

      !do i_coord = 1, 3
      !   if(i_coord .eq. j_coord ) then
      !     dielectric_constant(j_coord,j_coord) = 1.0d0 + dble(4*pi*dP_dE(i_coord,j_coord)/(cell_volume))
      !   else
      !     !shanghui note: the scalapack could only get diagonal terms, need to be extended later 
      !     !dielectric_constant(i_coord,j_coord) =  0.0d0
      !     dielectric_constant(i_coord,j_coord) = dble(4*pi*dP_dE(i_coord,j_coord)/(cell_volume))
      !   endif
      !enddo     ! i_coord 



     else ! lapack version 

      i_k_task = 0
      do i_k_point = 1,n_k_points, 1
      if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
      
         i_k_task = i_k_task + 1   
 
         call construct_matrix_complex(first_order_H_sparse, &
          first_order_H_complex(1:n_basis,1:n_basis,i_k_task), i_k_point)
 

      endif
     enddo ! n_k_point

     call evaluate_first_order_U_dielectric(& 
          Omega_MO(:,:,:,j_coord),first_order_H_complex(:,:,:),& 
          first_order_U_complex(:,:,:,j_coord)) 
          ! here i_k_point is just to cheat compile@thnec to use -O3, because it need to write out something,
          ! otherwise give NaN. 


     endif ! scalapack 

        call get_times(time_Sternheimer, clock_time_Sternheimer, &
        &              tot_time_Sternheimer, tot_clock_time_Sternheimer)


    endif 

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

      total_number_of_loops = total_number_of_loops + number_of_loops



     if(use_dfpt_pulay) then  
     !------begin add for pulay_mixing-----------
     call cleanup_pulay_mixing()
     !------end add for pulay_mixing-----------
     endif 
  enddo    ! j_coord 



   call get_timestamps(time_Hessian, clock_time_Hessian)

   if(.not.use_scalapack) then
!------------begin postprocessing(mpi version): dielectric_constant(all coords)----------------
     !call get_timestamps(time_Hessian, clock_time_Hessian)
     call evaluate_dielectric_constant &
        (occ_numbers, Omega_MO, first_order_U_complex, &  
         dielectric_constant)
!------------end postprocessing(mpi version): dielectric_constant(all coords)----------------
   else
      call  sync_vector_complex(dP_dE, 3*3)

      write(info_str,'(A)') 'dP_dE (Bohr^3) at every cycle:--->'
      call localorb_info(info_str, use_unit,'(A)', OL_norm)
      write(info_str,*) dP_dE(1,1:3)
      call localorb_info(info_str, use_unit,'(A)', OL_norm)
      write(info_str,*) dP_dE(2,1:3)
      call localorb_info(info_str, use_unit,'(A)', OL_norm)
      write(info_str,*) dP_dE(3,1:3)
      call localorb_info(info_str, use_unit,'(A)', OL_norm)
     !if (myid.eq.0) then
     ! do i_coord=1,3
     !   !write(use_unit,*)  (dP_dE(i_coord,j_coord),j_coord=1,3)
     !   write(use_unit,*)  dP_dE(i_coord,:)
     ! enddo
      write (info_str,'(A)') 'DFPT polarizability (Bohr^3)        xx        yy        zz        xy        xz        yz'
      call localorb_info(info_str, use_unit,'(A)', OL_norm)
      write (info_str,'(2X,A,1X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3,2X,F8.3)') &
      '| Polarizability:--->          ', real(dP_dE(1,1)), real(dP_dE(2,2)), &
      & real(dP_dE(3,3)), real(dP_dE(1,2)), real(dP_dE(1,3)), real(dP_dE(2,3))
      call localorb_info(info_str, use_unit,'(A)', OL_norm)

     !endif

     do i_coord = 1, 3 
     do j_coord = 1, 3
        if(i_coord .eq. j_coord ) then 
          dielectric_constant(i_coord,j_coord) = 1 + dble(4*pi*dP_dE(i_coord,j_coord)/(cell_volume))
        else 
          dielectric_constant(i_coord,j_coord) =  dble(4*pi*dP_dE(i_coord,j_coord)/(cell_volume))
        endif
     enddo     ! i_coord 
     enddo     ! j_coord 

   endif ! lapack/scalapack
 
   write(info_str,'(A)') ''
   call localorb_info(info_str, use_unit,'(A)', OL_norm  )
   write(info_str,'(A)') 'DFPT for dielectric_constant:--->'
   call localorb_info(info_str, use_unit,'(A)', OL_norm)
   write(info_str,*) dielectric_constant(1,1:3)
   call localorb_info(info_str, use_unit,'(A)', OL_norm)
   write(info_str,*) dielectric_constant(2,1:3)
   call localorb_info(info_str, use_unit,'(A)', OL_norm)
   write(info_str,*) dielectric_constant(3,1:3)
   call localorb_info(info_str, use_unit,'(A)', OL_norm)

    call get_times(time_Hessian, clock_time_Hessian, &
        &              tot_time_Hessian, tot_clock_time_Hessian)

    write(info_str,'(A)') ''
    call localorb_info(info_str, use_unit,'(A)', OL_norm  )
    call output_timeheader(deffmt, info_str, OL_norm)
    call output_times(deffmt, "Time for  dielectric calculation", &
        &                 time_Hessian, clock_time_Hessian, OL_norm)

    write(info_str,'(A)') "==========================================================================="
    call localorb_info(info_str, use_unit,'(A)', OL_norm  )


!-------------------shanghui begin deallocate------------------------------------
!----------------------(1) grid------------------------------------
       if(allocated( first_order_rho ))    call aims_deallocate(first_order_rho,  "first_order_rho" )
       if(allocated( first_order_potential ))    call aims_deallocate(first_order_potential,  "first_order_potential" )

!----------------------(2) matrix-----------------------------------------
       if(allocated( first_order_H_sparse ))    call aims_deallocate(first_order_H_sparse,  "first_order_H_sparse" )
       if(allocated( momentum_matrix_sparse ))  call aims_deallocate(momentum_matrix_sparse,  "momentum_matrix_sparse" )
       if(allocated( first_order_density_matrix_sparse )) & 
          &          call aims_deallocate(first_order_density_matrix_sparse,  "first_order_density_matrix_sparse" )
       if(allocated( old_first_order_density_matrix_sparse ))  & 
          &          call aims_deallocate(old_first_order_density_matrix_sparse,  "old_first_order_density_matrix_sparse" )
               
       deallocate(first_order_U_complex)

       if(.not.use_scalapack)then
       deallocate(first_order_H_complex)
       deallocate(Omega_MO)
       deallocate(momentum_matrix_complex)
       endif 
      if(allocated(mm_scalapack))   call aims_deallocate(mm_scalapack, "+mm_scalapack")

!-------------------shanghui end deallocate------------------------------------


    end subroutine cpscf_solver_dielectric
!******
