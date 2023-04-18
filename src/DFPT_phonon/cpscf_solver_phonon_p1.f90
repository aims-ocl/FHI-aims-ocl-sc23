!****s* FHI-aims/cpscf_solver_p1
!  NAME
!    cpscf_solver
!  SYNOPSIS

    subroutine cpscf_solver_phonon_p1 &

    (converged)

!  PURPOSE
!  cpscf for phonon, using sparse storage of matrix for S,H and DM
!  shanghui 2014.11.11
!  
!  USES

      use constants, only: pi
      use dimensions
      use timing
      use physics
      use species_data
      use localorb_io
      use geometry
      use pbc_lists
      use runtime_choices
      use scalapack_wrapper, only : initialize_scalapack_DFPT_phonon, & 
                                    construct_overlap_supercell_scalapack, & 
                                    construct_hamiltonian_supercell_scalapack, & 
                                    solve_evp_supercell_scalapack, & 
                                    mxld_DFPT_phonon, mxcol_DFPT_phonon
      use synchronize_mpi

      use DFPT_phonon_supercell
      use debugmanager, only: module_is_debugged

      !---------begin add for pulay_mixing--------
      use runtime_choices, only : use_dfpt_pulay, dfpt_pulay_steps 
      !this dfpt_pulay_steps by default is set to 8, if set in control.in to 1, then equals to linear mix. 
      use DFPT_pulay_mixing,    only: pulay_mix, cleanup_pulay_mixing
      !--------end add for pulay_mixing-----------
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
      character*100 :: info_str
      logical :: below_it_limit

      character*8  :: cdate
      character*10 :: ctime
      character(*), parameter :: deffmt = '2X'

!---------------shanghui begin for vib-cal-------------
     real*8, parameter:: const_u          = 1.66053886d-27         
     real*8, parameter:: const_eV         = 1.60217653d-19
     real*8, parameter:: const_c          = 299792458d0
     real*8, parameter:: const_Angstr     = 1d-10
     
     integer  lwork, info  
     complex*16, allocatable :: workspace(:)
     real*8, allocatable :: r_workspace(:)
     real*8, allocatable :: eigenvalues(:)
     real*8, allocatable :: mass_vector(:)
     real*8, allocatable :: reduced_mass(:)
     real*8, allocatable :: frequencies(:)
     complex*16, allocatable :: dynamical_matrix(:,:)     
     real*8   hessian_factor
     real*8, external :: ddot
     complex*16     :: buf   

!---------------shanghui end for vib-cal----------------



!-------------------shanghui add for DFPT------------------------------------
!------------------------(1) grid ---------------------------------------
     real*8, allocatable :: first_order_rho_circle(:,:,:)
     real*8, allocatable :: first_order_potential_circle(:,:,:)
     real*8, allocatable :: first_order_rho_cubic(:,:,:)
     real*8, allocatable :: first_order_potential_cubic(:,:,:)

     real*8, allocatable :: rho_free_gradient(:,:,:)
     real*8, allocatable :: v_free_gradient(:,:,:)

!------------------------(2) matrix -----------------------------------
     real*8, allocatable :: first_order_S_sparse(:,:,:)
     real*8, allocatable :: first_order_H_sparse(:,:,:)

     real*8, allocatable ::  density_matrix_sparse(:)
     real*8, allocatable ::  first_order_density_matrix_sparse(:,:,:)
     real*8, allocatable ::  old_first_order_density_matrix_sparse(:,:,:)

     real*8, allocatable ::  energy_density_matrix_sparse(:) 
     real*8, allocatable ::  first_order_energy_density_matrix_sparse(:,:,:)

!------------------shanghui end add for DFPT------------------------------------- 





      ! counters

      integer :: i_spin
      integer :: i_atom,j_atom, i_coord, j_coord,i_coord_2, i_center,i_center_read,j_center, i_center_in_hamiltonian
      integer :: i_basis, j_basis,  i_point,i_state, i_k_point,i_k_task,i_q_point, i_compute


      character(*), parameter :: func = 'cpscf_solver'
      
      logical first_iteration


     !--------shanghui begin for Hessian------------------------------------------------------------------------
      logical, parameter :: only_free_hartree_test = .false.
      logical, parameter :: read_fd_DM = .false.
      logical, parameter :: output_phonon_band = .true.  !.false.


      logical DM1_restart_file_exists

      integer :: i_basis_read, i_k_point_read
      real*8  change_of_first_order_DM
      real*8  time_start,time_end
      real*8  temp1, temp2
      real*8 :: left(n_hamiltonian_matrix_size),right(n_hamiltonian_matrix_size)

      real*8  hessian(3,n_centers_in_sc_DFPT,3,n_atoms)
      real*8  hellman_feynman_hessian_free_part(3,n_centers_in_sc_DFPT,3,n_atoms)
      real*8  hellman_feynman_hessian_delta_part(3,n_centers_in_sc_DFPT,3,n_atoms)
      real*8  pulay_hessian(3,n_centers_in_sc_DFPT,3,n_atoms)
     !--------shanghui end for Hessian------------------------------------------------------------------------

     !--------shanghui begin for q-point------------------
      real*8 q_x, q_y, q_z, delta_band(1:3), band_distance, delta_band_length
      integer i_cell, i_band, i_cell_x, i_cell_y, i_cell_z
      complex*16, dimension(:), allocatable:: q_phase_band 
     !--------shanghui end for q-point-------------------


!-------------------shanghui begin add for DFPT------------------------------------
!----------------------(1)grid-------------------------------------------
       allocate(first_order_rho_circle(3,n_centers_in_sc_DFPT,n_full_points))
       allocate(first_order_potential_circle(3,n_centers_in_sc_DFPT,n_full_points))
       allocate(first_order_rho_cubic(3,n_centers_in_sc_DFPT,n_full_points))
       allocate(first_order_potential_cubic(3,n_centers_in_sc_DFPT,n_full_points))

       allocate(rho_free_gradient(3,n_centers_in_sc_DFPT,n_full_points))
       allocate(v_free_gradient(3,n_centers_in_sc_DFPT,n_full_points))

!----------------------(2)matrix-----------------------------------------
       allocate(first_order_S_sparse(3, n_centers_in_sc_DFPT,n_hamiltonian_matrix_size))
       allocate(first_order_H_sparse(3, n_centers_in_sc_DFPT,n_hamiltonian_matrix_size))

       allocate(density_matrix_sparse(n_hamiltonian_matrix_size))
       allocate(first_order_density_matrix_sparse(3, n_centers_in_sc_DFPT,n_hamiltonian_matrix_size))
       allocate(old_first_order_density_matrix_sparse(3, n_centers_in_sc_DFPT, n_hamiltonian_matrix_size))

       allocate(energy_density_matrix_sparse(n_hamiltonian_matrix_size))
       allocate(first_order_energy_density_matrix_sparse(3,n_centers_in_sc_DFPT,n_hamiltonian_matrix_size))


       if(use_scalapack_DFPT_phonon) then
         call initialize_scalapack_DFPT_phonon()    ! all the local array are allocatable in this routine

         allocate(overlap_supercell(1,1))
         allocate(hamiltonian_supercell(1,1))
         allocate(KS_eigenvalue_supercell(1))
         allocate(KS_eigenvector_supercell(1,1))
       else 
         allocate(overlap_supercell(n_basis_sc_DFPT,n_basis_sc_DFPT))
         allocate(hamiltonian_supercell(n_basis_sc_DFPT,n_basis_sc_DFPT))
         allocate(KS_eigenvalue_supercell(n_basis_sc_DFPT))
         allocate(KS_eigenvector_supercell(n_basis_sc_DFPT,n_basis_sc_DFPT))
       endif
!-------------------shanghui end add for DFPT-----------------------------------


!-------------------shanghui begin for vib-cal-------------
       lwork=9*n_atoms
       allocate( workspace(lwork) ) 
       allocate( r_workspace(lwork) ) 
       allocate( eigenvalues(3*n_atoms) )
       allocate( mass_vector(3*n_atoms) )
       allocate( reduced_mass(3*n_atoms) )
       allocate( frequencies(3*n_atoms) )
       allocate( dynamical_matrix(3*n_atoms,3*n_atoms) )   
!-------------------shanghui end for vib-cal-------------


       if(packed_matrix_format.ne.PM_index) then
       call aims_stop('shanghui only use sparse matrix for DFPT_phonon', 'cpscf_solver_phonon_p1')
       endif 
   
    if(myid.eq.0) then 
    write(info_str, *)'===============begin memory(GB)/core needs for DFPT_phonon=================='
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    write(info_str,'(a,i6,i6)') 'n_centers_in_sc_DFPT,n_hamiltonian_matrix_size ', &
                      n_centers_in_sc_DFPT, n_hamiltonian_matrix_size
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )


    write(info_str,'(a,f15.6,a)') 'first_order_S_sparse like matrixes                      (total_number = 5):', & 
          5*8*3*dble(n_centers_in_sc_DFPT)*dble(n_hamiltonian_matrix_size)/2**30,  ' GB'
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    write(info_str,'(a,f15.6,a)') 'first_order_rho like matrixes                           (total_number = 6):', & 
          6*8*3*dble(n_centers_in_sc_DFPT)*dble(n_full_points)/2**30,  ' GB'
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    write(info_str,'(a,f15.6,a)') 'lapack version overlap_supercell like matrixes          (total_nubmer = 9):', & 
          9*8*3*dble(n_basis_sc_DFPT)*dble(n_basis_sc_DFPT)/2**30,  ' GB'
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    write(info_str,'(a,f15.6,a)') 'scalapack version ovlp_supercell_scalapack like matrixes(total_number =16):', & 
          16*8*3*dble(mxld_DFPT_phonon)*dble(mxcol_DFPT_phonon)/2**30,  ' GB'
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )

!    write(info_str,'(a,f12.3,a)') 'old second_order_S_sparse: ', &
!     8*3*dble(n_centers_in_sc_DFPT)*3*dble(n_centers_in_sc_DFPT)*dble(n_hamiltonian_matrix_size)/2**30,  'GB'
!    call localorb_info ( info_str,use_unit,'(a)', ol_norm )
!   write(info_str,'(a,2(f12.3,a))') ' temp_second_order_S: ', &
!               8*3*dble(n_centers_in_sc_DFPT)*3*dble(n_centers_in_sc_DFPT)/2**30,  'GB'
!    call localorb_info ( info_str,use_unit,'(a)', ol_norm )
    write(info_str, *)'===============end memory(GB)/core needs for DFPT_phonon=================='
    call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    endif

       
      ! begin work
      converged= .false.
      number_of_loops = 0
      below_it_limit = (number_of_loops.lt.sc_iter_limit)
      first_iteration = .true.
 
      if(output_phonon_band.and.myid.eq.0) then
      open(unit=87,file='shanghui-FHI-aims-Hessian.dat')
      open(unit=88,file='shanghui-FHI-aims-band_structure.dat')
      endif


  if (module_is_debugged("DFPT")) then
  if(myid.eq.0) then
   !-----begin write down data for my get_right_hessian_to_band.f90, need to be remove later---
   write(use_unit,*) 'mass:'
   do i_atom = 1 ,n_atoms 
   write(use_unit,*)  dsqrt(species_m(species(i_atom)))
   enddo 
   write(use_unit,*) 'cell_index_sc_DFPT:'
   do i_cell = 1, n_cells_in_sc_DFPT
     write(use_unit,*)  cell_index_sc_DFPT(i_cell,1:3)
   enddo
   write(use_unit,*) '' 

   write(use_unit,*) 'center_in_sc_DFPT_to_atom:'
   do i_center = 1, n_centers_in_sc_DFPT
     write(use_unit,*) center_in_sc_DFPT_to_atom(i_center)
   enddo
   write(use_unit,*) '' 

   write(use_unit,*) 'center_in_sc_DFPT_to_cell_in_sc_DFPT.dat'
   do i_center = 1, n_centers_in_sc_DFPT
     write(use_unit,*) center_in_sc_DFPT_to_cell_in_sc_DFPT(i_center)
   enddo
   write(use_unit,*) '' 



   write(use_unit,*) 'shanghui begin test before S1----------:'
   write(use_unit,*) '(1)----------n_Cbasis-------------'
   write(use_unit,*) '============C_basis      center_of_Cbasis     index[1:3]=========='
   do i_compute = 1 ,n_Cbasis
     write(use_unit,'(2i5,a,3i5,a)') i_compute, Cbasis_to_center(i_compute), '[',cell_index(center_to_cell(Cbasis_to_center(i_compute)),1:3) , ']'
   enddo

   write(use_unit,*) '(2)----------n_centers_in_sc_DFPT-------------'
   write(use_unit,*) '============i_center_sc_DFPT      i_center     index[1:3]=========='
   do i_compute = 1 ,n_centers_in_sc_DFPT
     write(use_unit,'(2i5,a,3i5,a)') i_compute,center_in_sc_DFPT_to_center(i_compute), '[', cell_index_sc_DFPT( &
                   center_in_sc_DFPT_to_cell_in_sc_DFPT(i_compute) , 1:3), ']'
   enddo

   write(use_unit,*) '------------------------------------------------------------'
   write(use_unit,*) '(3)          DFPT-sc atomic structure '

    do i_cell = 1, 3, 1
          write(use_unit,'(2X,A,3(2X,F16.8),2X)') &
               "lattice_vector ", &
               (nsc_DFPT(i_cell)*lattice_vector(i_coord,i_cell)*bohr, i_coord=1,3,1)
    enddo
          write(use_unit,*) '  '

    do i_center = 1, n_centers_in_sc_DFPT, 1
      write(use_unit,'(12X,A,3(2X,F16.8),2X,A)') &
            "atom ", &
            (coords_center(i_coord,center_in_sc_DFPT_to_center(i_center))*bohr, i_coord=1,3,1), &
             species_name(species(center_in_sc_DFPT_to_atom(i_center)))
    enddo
   write(use_unit,*) '------------------------------------------------------------'

   write(use_unit,*) '-----------------------------------------------------------'
   write(use_unit,*) ' (4)          k_weights(i_k_point)                           '
       do i_k_point = 1,n_k_points 
          write(use_unit,'(i5,3F9.5,F9.5)') i_k_point, k_point_list (i_k_point,1:3), k_weights(i_k_point) 
       enddo 
   write(use_unit,*) '------------------------------------------------------------'

   write(use_unit,*) 'shanghui end test before S1----------:'

   endif

   endif  ! debugmanager
   !-----end write down data for my get_right_hessian_to_band.f90, need to be remove later---
 

  
!--------------CPSCF (0.1) begin first_order_S-------------------
     call get_timestamps(time_first_order_S, clock_time_first_order_S)
     call integrate_first_order_S_p1(partition_tab, l_shell_max,& 
          first_order_S_sparse)

     if (module_is_debugged("DFPT")) then
     if(myid.eq.0) then
       write(use_unit,*) 'shanghui for python: 1x'
       write(use_unit,'(39f20.15)') (first_order_S_sparse(1,1,i_basis),i_basis=1,39)
       write(use_unit,*) ' '
       write(use_unit,*) 'shanghui for python: 2x'
       write(use_unit,'(39f20.15)') (first_order_S_sparse(1,2,i_basis),i_basis=1,39)
       write(use_unit,*) ' '
       write(use_unit,*) 'shanghui for python: 15x'
       write(use_unit,'(39f20.15)') (first_order_S_sparse(1,15,i_basis),i_basis=1,39)
     endif
     endif
 
     call get_times(time_first_order_S, clock_time_first_order_S, &
     &              tot_time_first_order_S, tot_clock_time_first_order_S)
!--------------CPSCF (0.1) end first_order_S---------------   



!--------------CPSCF (0.2) begin calculate gradient_free_V-------------------------------
    call  integrate_free_atom_sum_gradient_p1 &
          (partition_tab, rho_free_gradient,v_free_gradient) 
!--------------CPSCF (0.2) end calculate gradient_free_V-------------------------------

!---------------CPSCF (0.3) shanghui test how to calcualte DM_sparse---------------------------------
    if(use_scalapack) then  
      !  (1) in lapack version we will add k_weights by ourself, so only scalapck version need this.
      !  (2) Here we use 'use_scalapack' is because this kweight_occs only matters for DM0 and EDM0,
      !  which is related to the KS_solver of DFT. 
      call kweight_occs('cpscf_solver_phonon_p1', occ_numbers) 
    endif 

    call evaluate_zero_order_DM_p1 & 
          (KS_eigenvector, KS_eigenvector_complex, occ_numbers,density_matrix_sparse)

    !---------shanghui begin test DM0 supercell----------------------------
    !call trans_sparse_to_matrix(density_matrix_sparse,overlap_supercell)
    ! write(info_str, *)'===============(1) DM0 from sparse ================='
    !call localorb_info ( info_str,use_unit,'(A)', OL_norm )
    ! if(myid.eq.0) then 
    ! do i_basis=1,n_basis_sc_DFPT
    ! write(use_unit,'(90f9.5)') (overlap_supercell(i_basis,j_basis),j_basis=1,n_basis_sc_DFPT)
    ! enddo 
    ! write(use_unit,'(40f20.15)') (density_matrix_sparse(i_basis),i_basis=1,n_hamiltonian_matrix_size-1)
    ! endif
    !---------shanghui end test DM0 supercell----------------------------

 
  if (module_is_debugged("DFPT")) then
     if(myid.eq.0) then
      write(use_unit,*) '************zero_order_dm(atom1,X)****************'
      write(use_unit,'(39f20.15)') (density_matrix_sparse(i_basis),i_basis=1,39)
      write(use_unit,*) '************zero_order_dm(atom1,X)****************'
     endif
  endif 
!---------------CPSCF (0.3) shanghui end test how to calcualte DM_sparse-----------------------------

!--------------CPSCF (0.4) begin calculate init_first_order_DM1-------------------------------
     call cpu_time(time_start)

     if(use_scalapack_DFPT_phonon) then
    ! scalapack version 
       call construct_overlap_supercell_scalapack(overlap_matrix)
       call construct_hamiltonian_supercell_scalapack(hamiltonian)
       call solve_evp_supercell_scalapack() 
     else 
    ! lapack version 
       call trans_sparse_to_matrix(overlap_matrix,overlap_supercell)
       call trans_sparse_to_matrix(hamiltonian(:,1),hamiltonian_supercell)
       call solve_eigen_supercell(overlap_supercell,hamiltonian_supercell, &   
            KS_eigenvalue_supercell,KS_eigenvector_supercell)

       if (module_is_debugged("DFPT")) then
        if(myid.eq.0) then 
           write(use_unit,*) 'lapack KS_eigenvalue_supercell:', KS_eigenvalue_supercell(:)   
        endif
       endif

     endif 

     call cpu_time(time_end)
     if(myid.eq.0) write(use_unit,*) '@@@@@@time for solve_eigen_supercell@@@@@@@',time_end-time_start

!-------------------shanghui begin old_first_order_density_matrix---------------------------
        first_order_H_sparse = 0.0d0
        call evaluate_first_order_DM_supercell_p1(first_order_S_sparse,first_order_H_sparse,  &  
             old_first_order_density_matrix_sparse)
!-------------------shanghui end calculate old_first_order_density_matrix-------------------
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

        write(info_str,*) 'DFPT_mixing = ', DFPT_mixing
        call localorb_info(info_str, use_unit,'(A)', OL_norm )
        write(info_str,*) 'DFPT_sc_accuracy_dm = ', DFPT_sc_accuracy_dm
        call localorb_info(info_str, use_unit,'(A)', OL_norm )
        call date_and_time(cdate, ctime)
        write(info_str,'(2X,A,A,A,A)') "Date     :  ", cdate, ", Time     :  ", ctime
        call localorb_info(info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(A)') "------------------------------------------------------------"
        call localorb_info(info_str, use_unit,'(A)', OL_norm )

        call get_timestamps ( time_cpscf_loop, clock_time_cpscf_loop )

!----------CPSCF: (1) Begin  update first_order_rho-----------------------------------
       call get_timestamps(time_first_order_DM, clock_time_first_order_DM)
       call cpu_time(time_start)

       call evaluate_first_order_DM_supercell_p1(first_order_S_sparse,first_order_H_sparse,  &  
             first_order_density_matrix_sparse)

      if (module_is_debugged("DFPT")) then
       if(myid.eq.0) then 
       write(use_unit,*) '************CPSCF: first_order_dm(atom1,X)****************'
       write(use_unit,*) 'shanghui for python: 1x'
       write(use_unit,'(39f20.15)') (first_order_density_matrix_sparse(1,1,i_basis),i_basis=1,39)
       write(use_unit,*) '************CPSCF: first_order_dm(atom1,X)****************'
       endif
      endif 

    !shanghui: we need to perform density matrix mixing here
         change_of_first_order_DM =0.0d0         

         do i_basis = 1,n_hamiltonian_matrix_size
            do i_center = 1,n_centers_in_sc_DFPT
            do i_coord  = 1,3
            change_of_first_order_DM =                &
         max( change_of_first_order_DM,             &
         abs( first_order_density_matrix_sparse(i_coord,i_center,i_basis)  &
            - old_first_order_density_matrix_sparse(i_coord,i_center,i_basis))  )
            enddo
            enddo
         enddo

         if(use_dfpt_pulay) then  

            !------begin add for pulay_mixing-----------
            if(number_of_loops.eq.1) then
            ! shanghui add note: for loop=1, we need to inital X_in to  
            ! first_order_density_matrix_sparse, otherwise, the X_in 
            ! will set to 0.0d0, this is fine for electric_field response, 
            ! however, in phonon, X_in is not 0.0d0, we need to set here.   
            call pulay_mix(first_order_density_matrix_sparse, 3*n_centers_in_sc_DFPT*n_hamiltonian_matrix_size, &
                 number_of_loops, dfpt_pulay_steps, &
                 1.0d0 ) 
            else ! > 1
            call pulay_mix(first_order_density_matrix_sparse, 3*n_centers_in_sc_DFPT*n_hamiltonian_matrix_size, &
                 number_of_loops-1, dfpt_pulay_steps, &
                 DFPT_mixing )
            endif
            old_first_order_density_matrix_sparse = first_order_density_matrix_sparse
            !------end add for pulay_mixing-----------

         else
 
         do i_basis = 1,n_hamiltonian_matrix_size
            do i_center = 1,n_centers_in_sc_DFPT
            do i_coord  = 1,3

         first_order_density_matrix_sparse(i_coord,i_center,i_basis) =       &
         (1.0d0-DFPT_mixing)*old_first_order_density_matrix_sparse(i_coord,i_center,i_basis)+  &
         DFPT_mixing*first_order_density_matrix_sparse(i_coord,i_center,i_basis)
     
         old_first_order_density_matrix_sparse(i_coord,i_center,i_basis) =   &
         first_order_density_matrix_sparse(i_coord,i_center,i_basis)
            enddo
            enddo
         enddo

         endif ! use_dfpt_pulay

      !-------shanghui begin parallel------          
        if(myid.eq.0) then 
        write(use_unit,*) "((((((((((((((((((((((((((((((((("
        write(use_unit,*) change_of_first_order_DM
        write(use_unit,*) ")))))))))))))))))))))))))))))))))"
        endif
      !-------shanghui end parallel------          
      


!-------------------restart_read_DM1 code-----------------------
     !ref: restart_scalapack_read
     inquire(FILE='restart_DM1.dat',EXIST=DM1_restart_file_exists)
     call bcast_logical(DM1_restart_file_exists, 0)
 
     if( restart_read_DM1 .and.first_iteration.and.DM1_restart_file_exists) then

        if(myid.eq.0) then 
        write(info_str,'(2X,2A)') 'reading restart information from file'
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

        open(file = 'restart_DM1.dat', unit = 58,  & 
             status = 'old', form = 'unformatted',action='read')

        do i_center=1,n_centers_in_sc_DFPT
        do i_basis =1 ,n_hamiltonian_matrix_size
           read(58) first_order_density_matrix_sparse(1:3,i_center,i_basis)

           old_first_order_density_matrix_sparse(1:3,i_center,i_basis) =   &
               first_order_density_matrix_sparse(1:3,i_center,i_basis)
        enddo
        enddo
 
        close(58) 
        else
        first_order_density_matrix_sparse = 0.0d0 
        old_first_order_density_matrix_sparse = 0.0d0
        endif

        do i_coord = 1,3
        do i_center = 1, n_centers_in_sc_DFPT
           call sync_sparse_matrix( first_order_density_matrix_sparse(i_coord,i_center,:) )
           call sync_sparse_matrix( old_first_order_density_matrix_sparse(i_coord,i_center,:) )
        enddo
        enddo
 
     endif
 
!-------------------restart_read_DM1 code-----------------------

!-------------------restart_write_DM1 code-----------------------
     if( restart_write_DM1 .and. myid.eq.0) then

        write(info_str,'(2X,2A)') 'writing restart information to file'
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

        open(file = 'restart_DM1.dat', unit = 58,  & 
             status = 'replace', form = 'unformatted',action='write')

        do i_center=1,n_centers_in_sc_DFPT
        do i_basis =1 ,n_hamiltonian_matrix_size
           write(58) first_order_density_matrix_sparse(1:3,i_center,i_basis)
        enddo
        enddo
 
        close(58) 

     endif 
!-------------------restart_write_DM1 code-----------------------


   
!---------------check if first_order_density_matrix_sparse is right for every_center---------- 
      if(read_fd_DM.and.(.not.first_iteration)) then
         if(myid.eq.0) write(use_unit,*) 'check if my first_order_density_matrix_sparse is right compared with read'
      do i_center =1,n_centers_in_sc_DFPT 
         change_of_first_order_DM = 0 
         do i_basis =1 ,n_hamiltonian_matrix_size
         change_of_first_order_DM =                &
         max( change_of_first_order_DM,             &
         abs( first_order_density_matrix_sparse(1,i_center,i_basis)  &
            - old_first_order_density_matrix_sparse(1,i_center,i_basis))  )
         enddo 
         if(myid.eq.0) write(use_unit,*) i_center,change_of_first_order_DM 
      enddo

         call aims_stop('The  H1 --> DM1 test is finished at cpscf cycle 2') 
      endif 
!---------------end check if first_order_density_matrix_sparse is right for every_center---------- 
 

     if( read_fd_DM .and. first_iteration) then
     !-------for read_fd_DM-----------------
     open(unit=57,file='fd_DM')
     first_order_density_matrix_sparse = 0.0d0
     do i_center=1,n_centers_in_sc_DFPT
        read(57,*) i_center_read
        if(myid.eq.0)  write(use_unit,*) i_center_read,i_center
        read(57,*) left(1:n_hamiltonian_matrix_size-1)
        read(57,*) right(1:n_hamiltonian_matrix_size-1)

        old_first_order_density_matrix_sparse(1,i_center,1:n_hamiltonian_matrix_size-1) =  &
        (right(1:n_hamiltonian_matrix_size-1)-left(1:n_hamiltonian_matrix_size-1))/(0.02d0/0.52917721)

        first_order_density_matrix_sparse(1,i_center,1:n_hamiltonian_matrix_size-1) =  &
        (right(1:n_hamiltonian_matrix_size-1)-left(1:n_hamiltonian_matrix_size-1))/(0.02d0/0.52917721)
     enddo

      if(myid.eq.0) then 
      write(use_unit,*) '************shanghui readed first_order_dm(atom1,X)****************'
      write(use_unit,*) 'shanghui for python: 1x'
      write(use_unit,'(39f20.15)') (first_order_density_matrix_sparse(1,1,i_basis),i_basis=1,39)
      write(use_unit,*) '************shanghui readed first_order_dm(atom1,X)****************'
      endif
     !-------end for read_fd_DM-----------------
     endif 

     call get_times(time_first_order_DM, clock_time_first_order_DM, &
        &           tot_time_first_order_DM, tot_clock_time_first_order_DM)

!--------CPSCF : (1) end first-order-density update and mixing--------

!--------CPSCF : (2) begain to calculate first_order_H-----------------
        !----in pbc case, we have used density matrix instead of orbitals, it looks neat. 

        first_iteration = .false.

        
        call get_timestamps(time_first_order_density, clock_time_first_order_density)

        call integrate_first_order_rho_p1(partition_tab, l_shell_max,  &
             density_matrix_sparse,first_order_density_matrix_sparse, &
             first_order_rho_cubic)
        call trans_cubic_to_circle(partition_tab, & 
             first_order_rho_cubic, first_order_rho_circle)

          do i_center =1, n_centers_in_sc_DFPT
          do i_coord=1,3
 
             first_order_rho_circle(i_coord,i_center,1:n_full_points)=  &       !@   delta_rho
             first_order_rho_circle(i_coord,i_center,1:n_full_points)+  &       !@ = rho
             rho_free_gradient(i_coord,i_center,1:n_full_points)         !@ - free_rho
 
          enddo 
          enddo 


        call get_times(time_first_order_density, clock_time_first_order_density, &
        &              tot_time_first_order_density, tot_clock_time_first_order_density)


        call get_timestamps(time_first_order_potential, clock_time_first_order_potential)

 
          first_order_potential_circle=0.0d0 
 

         do i_center = 1, n_centers_in_sc_DFPT
              
          do i_coord=1, 3

          if(myid.eq.0) write(use_unit,*) 'i_center and i_coord:',i_center ,i_coord
           
            !less center test: if(i_center.eq.1.or.(i_center.ge.9.and.i_center.le.20)) then  
          !---------give rho^(1)(atom)-------------------------------------- 
          call update_hartree_potential_shanghui_p1 &
             ( i_center,hartree_partition_tab, & 
               first_order_rho_circle(i_coord,1:n_centers_in_sc_DFPT,1:n_full_points)) 

          !---------evevry rho^(1)(atom) give potential^(1)(atom+R)-------  
          call sum_up_whole_potential_shanghui_p1 &
             ( partition_tab, first_order_rho_circle(i_coord,i_center,1:n_full_points), &
               first_order_potential_circle(i_coord,i_center,1:n_full_points))  
 
           !endif ! end less i_center test
          enddo 
          enddo 
 
 
          do i_center =1, n_centers_in_sc_DFPT
          do i_coord=1,3
 
             first_order_potential_circle(i_coord,i_center, 1:n_full_points)= &    !@    pot  
             first_order_potential_circle(i_coord,i_center,1:n_full_points)   &    !@ = delta_pot
             +(-v_free_gradient(i_coord,i_center,1:n_full_points))          !@ + free_pot
 
          enddo 
          enddo 
 
        call get_times(time_first_order_potential, clock_time_first_order_potential, &
        &              tot_time_first_order_potential, tot_clock_time_first_order_potential)


        call get_timestamps(time_first_order_H, clock_time_first_order_H)

     !----------------begin free_part test--------------------------------------
     if(only_free_hartree_test) then 
        !---used for (1) hartree_free DFPT
        do i_center=1,n_centers_in_sc_DFPT
           do i_coord=1,3
              first_order_potential_circle(i_coord,i_center,1:n_full_points) =  &     
              -v_free_gradient(i_coord,i_center,1:n_full_points)   
           enddo
        enddo

        !---used for  (1)hartree_free fd benchmark (2) hartree_free DFPT 
         hartree_potential(:)= free_hartree_superpos(:)
     endif
     !----------------end free_part test---------------------------------------------


     call trans_circle_to_cubic(partition_tab, & 
          first_order_potential_circle, first_order_potential_cubic)

     call  integrate_first_order_H_p1 &
          (hartree_potential,first_order_potential_cubic, rho, rho_gradient,&
           partition_tab, l_shell_max, en_xc, en_pot_xc,    &
           density_matrix_sparse,first_order_density_matrix_sparse, & 
           first_order_H_sparse &
           )
      
     if (module_is_debugged("DFPT")) then
      if(myid.eq.0) then
      write(use_unit,*) '************first_order_H_sparse****************' 
      write(use_unit,*) 'shanghui for python: 1x'
      write(use_unit,'(39f20.15)') (first_order_H_sparse(1,1,i_basis),i_basis=1,39)
      write(use_unit,*) ' '
      write(use_unit,*) 'shanghui for python: 2x'
      write(use_unit,'(39f20.15)') (first_order_H_sparse(1,2,i_basis),i_basis=1,39)
      write(use_unit,*) ' '
      write(use_unit,*) 'shanghui for python: 15x'
      write(use_unit,'(39f20.15)') (first_order_H_sparse(1,15,i_basis),i_basis=1,39)
      endif
     endif

        call get_times(time_first_order_H, clock_time_first_order_H, &
         &              tot_time_first_order_H, tot_clock_time_first_order_H)
        
!--------CPSCF : (2) end to calculate first_order_H-----------------


!--------CPSCF : (3) begain to calculate first_order_U-----------------
      !call evaluate_first_order_U_supercell_p1(KS_eigenvalue_supercell,  &  
      !     KS_eigenvalue_supercell, first_order_S_sparse, first_order_H_sparse, & 
      !     first_order_U_supercell)


     !-----------shanghui begin for DFPT--------------------------
!      call evaluate_first_order_C_p1(first_order_H_sparse,  & 
!           first_order_S_sparse,  &
!           KS_eigenvector_complex, KS_eigenvalue, & 
!           occ_numbers,first_order_density_matrix_sparse) 
!      if(myid.eq.0) then 
!      write(use_unit,*) '************DFPT: first_order_dm(atom1,X)****************'
!      write(use_unit,*) 'shanghui for python: 1x'
!      write(use_unit,'(40f20.15)') (first_order_density_matrix_sparse(1,1,i_basis),i_basis=1,n_hamiltonian_matrix_size-1)
!      write(use_unit,*) '************DFPT: first_order_dm(atom1,X)****************'
!      endif
!     !-----------shanghui end for DFPT--------------------------

     !-----------shanghui begin for CPSCF---------------------- 
!     i_k_task = 0 
!      do i_k_point = 1,n_k_points, 1
!      if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
      
!         i_k_task = i_k_task + 1   

!      call construct_first_order_H_p1(   &
!            first_order_H_sparse,                       & 
!            first_order_H_complex(:,:,:,:,i_k_task),i_k_point)
!
!      if(myid.eq.0) then
!      write(use_unit,*) 'k=',i_k_point
      !write(use_unit,*) '************shanghui begain first_order_S_complex(atom1,X);****************'
      !do i_basis=1,n_basis
      !   write(use_unit,'(6f12.6)') (first_order_S_complex(1,1,i_basis,j_basis,i_k_point),j_basis=1,n_basis)
      !enddo
      !write(use_unit,*) '************shanghui enddo first_order_S_complex(atom1,X)****************'


      !write(use_unit,*) '************shanghui begin first_order_H_complex(atom1,X)****************'
      !do i_basis=1,n_basis
      ! write(use_unit,'(6f12.6)') (first_order_H_complex(1,1,i_basis,j_basis,i_k_point),j_basis=1,n_basis )
      !enddo
      !write(use_unit,*) '************shanghui end first_order_H_complex(atom1,X)****************'


      !write(use_unit,*) '************shanghui begain first_order_U(atom1,X)****************'
!     call evaluate_first_order_U_p1(first_order_H_complex(:,:,:,:,i_k_task),&
!          first_order_S_complex(:,:,:,:,i_k_task),    & 
!          KS_eigenvector_complex(:,:,:,i_k_task), KS_eigenvalue(:,:,i_k_point), &  
!         occ_numbers(:,:,i_k_point),  &
!          first_order_U_complex(:,:,:,:,i_k_task))
!
!      do i_basis=1,n_basis
!       write(use_unit,'(6f12.6)') (first_order_U_complex(1,1,i_basis,j_basis,i_k_point),j_basis=1,n_basis )
!      enddo
!      write(use_unit,*) '************shanghui end first_order_U(atom1,X)****************'
!      write(use_unit,*) ''
!
!      endif
!
!     endif
!     enddo

        call cpu_time(time_end)
        write(use_unit,*) '@@@@@@@@time for CPSCF@@@@@@@@@@@@@@@@@@',time_end-time_start

!--------------(3) end solve first_order_U problem----------------



! --------- Check convergence ----->>



!         check convergence of self-consistency loop
          if(read_fd_DM) then 
          converged = (number_of_loops.eq.2)
          else 
          converged = (change_of_first_order_DM.lt.DFPT_sc_accuracy_dm).and.(number_of_loops.ne.1)
          endif               

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

!       this is the end of the self-consistency cycle.

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
        write(info_str,'(A)') &
        "------------------------------------------------------------"
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

! << ---- end printing out data---------

  end do SCF_LOOP
! << ------ end self consistent cycle--------


        call get_timestamps(time_Hessian, clock_time_Hessian)

!-------------(1) for hellman_feynman_hessian--------------
         hellman_feynman_hessian_free_part=0.0d0
         hellman_feynman_hessian_delta_part=0.0d0

        !------(1.1) free_part------------
        call integrate_free_atom_hessian_p1(hellman_feynman_hessian_free_part)

        if (module_is_debugged("DFPT")) then
          if(myid.eq.0) then
           write(use_unit,*) ''
           write(use_unit,*) 'shanghui in cpscf.f90, for hellman_feynman_hessian(a.u.)------->free_part:'
           do i_center =1,n_centers_in_sc_DFPT
           write(use_unit,*) i_center,' x 1 x: ',hellman_feynman_hessian_free_part(1,i_center,1,1) 
           enddo  
           write(use_unit,*) 'shanghui in cpscf.f90, for hellman_feynman_hessian(a.u.)------->free_part:'
           do i_center =1,n_centers_in_sc_DFPT
           write(use_unit,*) i_center,' y 1 y: ',hellman_feynman_hessian_free_part(2,i_center,2,1)
           enddo
          endif
        endif


        !------(1.2) delta_part------------
        call integrate_first_order_rho_p1(partition_tab, l_shell_max,  &
             density_matrix_sparse,first_order_density_matrix_sparse, &
             first_order_rho_cubic)
        call trans_cubic_to_circle(partition_tab, & 
             first_order_rho_cubic, first_order_rho_circle)

          do i_center =1, n_centers_in_sc_DFPT
          do i_coord=1,3
 
             first_order_rho_circle(i_coord,i_center,1:n_full_points)=  &    !@   delta_rho 
             first_order_rho_circle(i_coord,i_center,1:n_full_points)+  &    !@ = rho  
             rho_free_gradient(i_coord,i_center,1:n_full_points)      !@ - rho_free
 
          enddo 
          enddo 


          do i_center=1, n_atoms
          do i_coord=1, 3

          !---------(1.2.1) get rho_multipole_supercell^(1)-------------------------------------- 
          call update_hartree_potential_shanghui_p1 &
             ( i_center,hartree_partition_tab, & 
               first_order_rho_circle(i_coord,1:n_centers_in_sc_DFPT,1:n_full_points)) 
          !---------(1.2.2) get V(0)^(1)(n_centers_in_sc_DFPT)-----------------  
          call update_hartree_potential_at_zero() 

          !---------(1.2.3) get hellman_feynman_hessian_delta_part-------------------- 
          call integrate_delta_hellman_hessian_p1  & 
            ( hellman_feynman_hessian_delta_part(1:3,1:n_centers_in_sc_DFPT,i_coord,i_center))
 
          enddo
          enddo

       
        if (module_is_debugged("DFPT")) then
         if(myid.eq.0) then
          write(use_unit,*) ''
          write(use_unit,*) 'shanghui in cpscf.f90, for hellman_feynman_hessian(a.u.)------->delta_part:'
          do i_center =1,n_centers_in_sc_DFPT
          write(use_unit,*) i_center,' x 1 x: ',hellman_feynman_hessian_delta_part(1,i_center,1,1) 
          enddo  
          write(use_unit,*) 'shanghui in cpscf.f90, for hellman_feynman_hessian(a.u.)------->delta_part:'
          do i_center =1,n_centers_in_sc_DFPT
          write(use_unit,*) i_center,' y 1 y: ',hellman_feynman_hessian_delta_part(2,i_center,2,1)
          enddo
         endif
        endif
 



!-------------(2) for pulay_hessian------------------------

      !------(2.1) prepare EDM(0),EDM(1) for Hessian----------

       call evaluate_zero_order_EDM_p1 & 
            (KS_eigenvector, KS_eigenvector_complex,KS_eigenvalue, occ_numbers, energy_density_matrix_sparse) 
      !-------shanghui begin parallel------
      ! if(myid.eq.0) then
      !  write(use_unit,*) '************zero_order_edm(atom1,X)****************'
      !  write(use_unit,'(40f20.15)') (energy_density_matrix_sparse(i_basis),i_basis=1,n_hamiltonian_matrix_size-1)
      !  write(use_unit,*) '************zero_order_edm(atom1,X)****************'
      ! endif
      !-------shanghui end parallel------

       call evaluate_first_order_EDM_supercell_p1(   & 
             first_order_S_sparse,first_order_H_sparse,  &  
             first_order_energy_density_matrix_sparse)

       !-------(2.2) get pulay_hessian------------------------

       call integrate_first_order_rho_p1(partition_tab, l_shell_max,  &
             density_matrix_sparse,first_order_density_matrix_sparse, &
             first_order_rho_cubic)

       call integrate_pulay_hessian_p1 &
           ( hartree_potential,  &  
             first_order_rho_cubic, first_order_potential_cubic, &
             rho, rho_gradient,  &
             partition_tab, l_shell_max,     &
             density_matrix_sparse, first_order_density_matrix_sparse,  &
             energy_density_matrix_sparse,first_order_energy_density_matrix_sparse,  &
             pulay_hessian &
            )
      

        if (module_is_debugged("DFPT")) then
        if(myid.eq.0) then
          write(use_unit,*) ''
          write(use_unit,*) 'shanghui in cpscf.f90, for pulay_hessian(a.u.)------->'
          do i_center =1,n_centers_in_sc_DFPT
          write(use_unit,*) i_center,' x 1 x: ', pulay_hessian(1,i_center,1,1) 
          enddo
          write(use_unit,*) 'shanghui in cpscf.f90, for pulay_hessian(a.u.)------->'
          do i_center =1,n_centers_in_sc_DFPT
          write(use_unit,*) i_center,' y 1 y: ', pulay_hessian(2,i_center,2,1)
          enddo
        endif
        endif 
         

     call get_times(time_Hessian, clock_time_Hessian, &
        &              tot_time_Hessian, tot_clock_time_Hessian)

     write(info_str,'(A)') ''
     call localorb_info(info_str, use_unit,'(A)', OL_norm  )
     call output_timeheader(deffmt, info_str, OL_norm)
     call output_times(deffmt, "Time for Hessian of this pertubation", &
        &                 time_Hessian, clock_time_Hessian, OL_norm)

     write(info_str,'(A)') "==========================================================================="
     call localorb_info(info_str, use_unit,'(A)', OL_norm  )

!---------------shanghui begin for vib-cal-------------
    allocate(q_phase_band(n_cells_in_sc_DFPT))


!  !-------shanghui begin output dynamical_matrix at q point------
!   do i_q_point = 1,n_k_points
! 
!    if(myid.eq.0) then
!       write(use_unit,*) '============================================='
!       write(use_unit,*) 'DFPT-Results: for q=',i_q_point
!       write(use_unit,*) '============================================='
!       write(use_unit,*)
!    endif
! 
!   q_phase_band(:) = k_phase(:,i_q_point) ! k_phase only works for n_cells_ham 
! 
! 
!      do i_center = 1, n_centers_in_sc_DFPT
!      do i_coord = 1, 3
!          do j_center = 1, n_atoms
!          do j_coord = 1 , 3 
! 
!             hessian(i_coord,i_center,j_coord,j_center) =  & 
!        (hellman_feynman_hessian_free_part(i_coord,i_center,j_coord,j_center))*hartree / (bohr*bohr)
! 
!          enddo  
!          enddo 
!      enddo 
!      enddo  
! 
!      call trans_hessian_to_dynamical_matrix(hessian, dynamical_matrix,q_phase_band)
! 
!    if(myid.eq.0) then
!      write(use_unit,*) ''
!      write(use_unit,*) 'helman_feynman_free in cpscf.f90:'
!      do i_atom = 1 , n_atoms
!      do i_coord = 1 , 3
! 
!         write(use_unit,*)  &
!             (((dynamical_matrix(3*i_atom+i_coord-3,3*j_atom+j_coord-3)) &
!                 , j_coord=1,1) ,j_atom=1,1)
!      enddo
!      enddo
!    endif
! 
! 
!      do i_center = 1, n_centers_in_sc_DFPT
!      do i_coord = 1, 3
!          do j_center = 1, n_atoms
!          do j_coord = 1 , 3 
! 
!             hessian(i_coord,i_center,j_coord,j_center) =  & 
!             !( hellman_feynman_hessian_free_part(i_coord,i_center,j_coord,j_center) + & 
!             !  hellman_feynman_hessian_delta_part(i_coord,i_center,j_coord,j_center)+ & 
!              ( pulay_hessian(i_coord,i_center,j_coord,j_center) )*hartree / (bohr*bohr)
! 
!          enddo  
!          enddo 
!      enddo 
!      enddo  
! 
!      call trans_hessian_to_dynamical_matrix(hessian, dynamical_matrix,q_phase_band)
! 
!    if(myid.eq.0) then
!      write(use_unit,*) ''
!      write(use_unit,*) 'pulay_hessian in cpscf.f90:'
!      do i_atom = 1 , n_atoms
!      do i_coord = 1 , 3
! 
!         write(use_unit,*)  &
!             (((dynamical_matrix(3*i_atom+i_coord-3,3*j_atom+j_coord-3)) &
!                 , j_coord=1,1) ,j_atom=1,1)
!      enddo
!      enddo
!    endif
! 
!   enddo  
!  !-------shanghui end output dynamical_matrix at q point------





        !-------shanghui begin parallel------
        if(myid.eq.0) then
      write(use_unit,*) ''
      write(use_unit,*) 'Total_hessian(ev/Ang*Ang) in cpscf.f90:'
        endif 
        !-------shanghui end parallel------

       do i_center = 1, n_centers_in_sc_DFPT
       do i_coord = 1, 3
          do j_center = 1, n_atoms
          do j_coord = 1 , 3 

             hessian(i_coord,i_center,j_coord,j_center) =  & 
             ( hellman_feynman_hessian_free_part(i_coord,i_center,j_coord,j_center) + & 
               hellman_feynman_hessian_delta_part(i_coord,i_center,j_coord,j_center)+ & 
               pulay_hessian(i_coord,i_center,j_coord,j_center) )*hartree / (bohr*bohr)

          enddo  
          enddo 
       enddo 
       enddo  

 
        !-------shanghui begin parallel------
        if(myid.eq.0) then
          write(use_unit,*) ''
          write(use_unit,*) 'Make Acoustic sum rule in cpscf.f90:'
        endif 
        !-------shanghui end parallel------


      !                                  ---
      !      d^2 E (supercell)           \         d^2 E (supercell)     
      !   -----------------------   =  - /      -----------------------  
      !   d u(J,alpha) d u(J,beta)       ----    d u(I,alpha) d u(J,beta) 
      !                                 (I.ne.J) 


      do i_coord = 1 , 3
      do j_coord = 1 , 3
 
        do j_center = 1 , n_atoms 
               hessian(i_coord,j_center, j_coord, j_center) = 0.0d0
        do i_center = 1, n_centers_in_sc_DFPT
 
             if(i_center.ne.j_center) then 
                ! here we directly compare i_center and j_center(j_atom) without using centers_in_hamiltonian 
                ! because:  
                ! in pbc_list, j_center always first run over n_atoms : 
                ! like this  centers_in_hamiltonian(i_atom)            = i_atom
               hessian(i_coord, j_center, j_coord, j_center)=  &
               hessian(i_coord, j_center, j_coord, j_center)-  &
               hessian(i_coord, i_center, j_coord, j_center)
            endif
 
        enddo   ! i_center
        enddo   ! j_center 
 
      enddo     ! j_coord
      enddo     ! i_coord
 
 
 
   !-------shanghui begin parallel------
   if(myid.eq.0) then

! code from phonopy, need to be remove later --------> 
!   f.write("# force constant matrix elements in eV/(Angstrom)^2 \n")
!    for (j_atom_count,j_atom) in enumerate(primitive.p2s_map):
!        for j_coord in range(3):
!            for i_atom in range(phonopy_obj.supercell.get_number_of_atoms()):
!                for i_coord in range(3):
!                    f.write("%18.10e " % Hessian[i_atom,j_atom,i_coord,j_coord])
!                f.write("  ")
!            f.write("\t # displaced atom %d, coord %d \n" % (j_atom+1,j_coord+1))
!    f.close()

     write(use_unit,*) ''
     write(use_unit,*) ''
     write(use_unit,*) '# output Hessian (ev/Ang*Ang), like phonopy-FHI-aims-Hessian.dat'

      do i_center = 1 , n_atoms
      do i_coord  = 1 , 3
         write(87,'(2500E18.10)')  &
            ((hessian(j_coord,j_center,i_coord,i_center), j_coord=1,3) ,j_center =1,n_centers_in_sc_DFPT)
         write(87,'(a,i5,a,i3)' )  '# displaced atom',i_center,'coord',i_coord
         write(87,*) ''
      enddo
      enddo



    write(use_unit,*) ''
    write(use_unit,*) ''
    write(use_unit,*) 'Get frequencies in cpscf.f90:'
   endif 
   !-------shanghui end parallel------

    ! factor to bring hessian into SI units .
    hessian_factor   = const_eV/(const_u*const_Angstr*const_Angstr)
    hessian(1:3, 1:n_centers_in_sc_DFPT, 1:3, 1:n_atoms) = & 
    hessian(1:3, 1:n_centers_in_sc_DFPT, 1:3, 1:n_atoms)  * hessian_factor 
 
    do i_atom=1, n_atoms 
       do i_coord = 1, 3
        mass_vector(3*i_atom+i_coord-3) = 1.0d0/dsqrt(species_m(species(i_atom)))
       end do
     end do
 
    
 
   !-------shanghui begin parallel------
   if(myid.eq.0) then
    write(use_unit,*) 'Solving eigenvalue system for Dynamical Matrix'
   endif 
   !-------shanghui end parallel------


 
     ! code from phonopy, need to be remove later --------> 
     !  dk = (kend-kstart)/(npoints-1)
     !   bands.append([(kstart + dk*n) for n in range(npoints)])
     !   dk_length = numpy.linalg.norm(dk)
     !for n in range(npoints):
     !       bands_distances.append(distance + dk_length*n)
     !
     ! distance += dk_length * (npoints-1)

    band_distance = 0.0d0 

    do i_band = 1, n_plot_band_DFPT_phonon

         delta_band(1:3) =  ( DFPT_phonon_band_end(i_band,1:3) - DFPT_phonon_band_begin(i_band,1:3)) & 
                            /real(n_points_in_band_DFPT_phonon(i_band) -1)   

         delta_band_length = dsqrt( delta_band(1)**2 + delta_band(2)**2 + delta_band(3)**2 ) 

         band_distance = band_distance - delta_band_length

    do i_q_point = 1 , n_points_in_band_DFPT_phonon(i_band)

             !----------begin q_phase----------------------
             q_x = 2.0d0*pi*(DFPT_phonon_band_begin(i_band,1) +  real(i_q_point-1) * delta_band(1))   
             q_y = 2.0d0*pi*(DFPT_phonon_band_begin(i_band,2) +  real(i_q_point-1) * delta_band(2))
             q_z = 2.0d0*pi*(DFPT_phonon_band_begin(i_band,3) +  real(i_q_point-1) * delta_band(3))

             band_distance = band_distance + delta_band_length

             do i_cell = 1, n_cells_in_sc_DFPT

                i_cell_x = cell_index_sc_DFPT(i_cell,1)
                i_cell_y = cell_index_sc_DFPT(i_cell,2)
                i_cell_z = cell_index_sc_DFPT(i_cell,3)

                ! q = 2*pi/R *i_q/N_q  , RN= N_cell*R ===> 
                ! phase_factor = i*q*RN = i*2*pi* [i_q/N_q] * [N_cell] 
                q_phase_band( i_cell) &
                = ( (exp((0.0d0,1.0d0)*q_x*i_cell_x))*   & 
                    (exp((0.0d0,1.0d0)*q_y*i_cell_y))*   &
                    (exp((0.0d0,1.0d0)*q_z*i_cell_z)) )

             end do ! i_cell
             !----------end q_phase----------------------


       call trans_hessian_to_dynamical_matrix(hessian, dynamical_matrix,q_phase_band)      
       
       ! symmetrize dynamical matrix 
       do i_coord = 1, 3*n_atoms
       do j_coord = 1, i_coord - 1 
          buf = (dynamical_matrix(i_coord,j_coord) + dconjg(dynamical_matrix(j_coord,i_coord)) )/2d0
          dynamical_matrix(i_coord,j_coord) = buf
          dynamical_matrix(j_coord,i_coord) = dconjg(buf)
       end do
       end do

       ! multiply with mass vector
       do i_coord = 1, 3*n_atoms, 1   
          dynamical_matrix(:,i_coord) = dynamical_matrix(:,i_coord)*mass_vector(i_coord)       
          dynamical_matrix(i_coord,:) = dynamical_matrix(i_coord,:)*mass_vector(i_coord) 
       end do


 
!--------------complex version works for all q point----------------------
      call ZHEEV('V','U',3*n_atoms,dynamical_matrix,3*n_atoms,eigenvalues,workspace,9*n_atoms,r_workspace,info)
!--------------end complex version works for all q point----------------------

    !---------shanghui begin output phonon band---------------
       do i_coord = 1, 3*n_atoms, 1
          if (eigenvalues(i_coord).gt.0d0) then
             frequencies(i_coord) = sqrt(eigenvalues(i_coord))
          else if (eigenvalues(i_coord).lt.0d0) then
             frequencies(i_coord) = -sqrt(-eigenvalues(i_coord))
          else 
             frequencies(i_coord) = 0d0
          end if
       enddo ! i_coord

    if(myid.eq.0) then

       write(88,'(100F15.6)') band_distance, (frequencies(i_coord)/(200*pi*const_c), i_coord = 1,3*n_atoms)

    endif 
    !---------shanghui end output phonon band---------------


   !-------shanghui begin test one point------
    if(myid.eq.0.and.i_q_point.eq.1.and.i_band.eq.1) then
       write(use_unit,*) '============================================='
       write(use_unit,*) 'DFPT-Results: for q1=',DFPT_phonon_band_begin(i_band,1:3) 
       write(use_unit,*) '============================================='
       write(use_unit,*)
       write(use_unit,*) 'List of all frequencies found:'
       write(use_unit,'(A13,A25)') 'Mode number','Frequency [cm^(-1)]'
     
       do i_coord = 1, 3*n_atoms, 1
          write(use_unit,'(I13,F25.8)') i_coord, frequencies(i_coord)/(200*pi*const_c)
       enddo ! i_coord

    endif 
   !-------shanghui end test one point------






  enddo ! i_q_point
     if(myid.eq.0)   write(88,*) '  '
  enddo ! i_band



    if( read_fd_DM ) then
     close(57)
    endif



    if(output_phonon_band.and.myid.eq.0) then
     close(87)
     close(88)
    endif


     if(use_dfpt_pulay) then  
    !------begin add for pulay_mixing-----------
     call cleanup_pulay_mixing()
    !------end add for pulay_mixing-----------
     endif

!---------------shanghui end for vib-cal-------------


!-------------------shanghui begin add for DFPT------------------------------------
!----------------------(1) grid------------------------------------
       deallocate(first_order_rho_circle)
       deallocate(first_order_potential_circle)
       deallocate(first_order_rho_cubic)
       deallocate(first_order_potential_cubic)

       deallocate(rho_free_gradient)
       deallocate(v_free_gradient)

!----------------------(2) matrix-----------------------------------------
!!!       deallocate(first_order_S_complex)
!!!       deallocate(first_order_H_complex)
!!!       deallocate(first_order_U_complex)


       deallocate(density_matrix_sparse)
       deallocate(first_order_density_matrix_sparse)
       deallocate(old_first_order_density_matrix_sparse)

       deallocate(energy_density_matrix_sparse)
       deallocate(first_order_energy_density_matrix_sparse)

       deallocate(overlap_supercell)
       deallocate(hamiltonian_supercell)
       deallocate(KS_eigenvalue_supercell)
       deallocate(KS_eigenvector_supercell)

!-------------------shanghui end add for DFPT------------------------------------

!-------------------shanghui begin for vib-cal-------------
       deallocate( workspace )
       deallocate( r_workspace )
       deallocate( eigenvalues )
       deallocate( mass_vector )
       deallocate( reduced_mass )
       deallocate( frequencies )
       deallocate(q_phase_band)
!-------------------shanghui end for vib-cal-------------
 

    end subroutine cpscf_solver_phonon_p1
!******
