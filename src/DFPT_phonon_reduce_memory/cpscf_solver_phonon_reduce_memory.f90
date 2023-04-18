!****s* FHI-aims/cpscf_solver_phonon_reduce_memory
!  NAME
!    cpscf_solver_phonon_reduce_memory
!  SYNOPSIS

    subroutine cpscf_solver_phonon_reduce_memory &
    (converged)

!  PURPOSE
!  an cpscf process for phonon_reduce_memory 
!  shanghui 2015.07
!  
!  USES

      use constants, only: pi
      use runtime_choices
      use dimensions
      use timing
      use physics
      use species_data
      use localorb_io
      use geometry
      use vdw_correction

      use scalapack_wrapper  , only : construct_first_order_hamiltonian_scalapack
      use pbc_lists , only : kweight_occs
      use synchronize_mpi, only : sync_vector_complex
      use debugmanager, only: module_is_debugged
      use mpi_tasks, only: aims_stop, myid, n_tasks

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
     complex*16, allocatable :: first_order_rho(:)
     complex*16, allocatable :: first_order_potential(:)

     real*8, allocatable :: first_order_rho_Re(:)
     real*8, allocatable :: first_order_rho_Im(:)
     real*8, allocatable :: first_order_potential_Re(:)
     real*8, allocatable :: first_order_potential_Im(:)

     complex*16, allocatable :: rho_free_gradient(:)
     complex*16, allocatable :: v_free_gradient(:)

!------------------------(2) matrix -----------------------------------
!   now is mpi version.
!   for scalapck version in the future, we just need smaller (la,lb) for first_order_S_complex
     complex*16, allocatable :: first_order_S_sparse(:)
     complex*16, allocatable :: first_order_S_complex(:,:,:)
     complex*16, allocatable :: first_order_H_sparse(:,:)
     complex*16, allocatable :: first_order_H_complex(:,:,:) 

     complex*16,allocatable :: first_order_U_complex(:,:,:)

     real*8, allocatable ::  density_matrix_sparse(:)
     complex*16, allocatable ::  first_order_density_matrix_sparse(:)

     real*8, allocatable ::  energy_density_matrix_sparse(:) 
     complex*16, allocatable ::  first_order_energy_density_matrix_sparse(:)


!------------------------(3) frequencies------------------------------
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
     complex*16, allocatable :: hellman_feynman_dynamical_matrix_free_part(:,:,:,:)
     complex*16, allocatable :: hellman_feynman_dynamical_matrix_delta_part(:,:,:,:)
     complex*16, allocatable :: pulay_dynamical_matrix(:,:,:,:)
     complex*16, allocatable :: dynamical_matrix(:,:)
     real*8  hessian_factor
     complex*16 buf
     real*8, dimension(3, n_atoms)           :: hellman_feynman_term_at_gamma
	 
	 real*8 :: vdw_energy 
	 real*8 :: vdw_forces(3,n_atoms)
	 real*8 :: vdW_Hessian(3*n_atoms,3*n_atoms)

!------------------------(4) Born-effective-charges------------------------------
     complex*16, allocatable :: Omega_MO(:,:,:,:)
    !    = < i(k)|-r|j(k) > 
    !<1> = (C^+ momentum_matrix_complex C)_ij/( Ei(k) - Ej(k) ) 
     real*8, allocatable     :: momentum_matrix_sparse(:)
     complex*16, allocatable :: momentum_matrix_complex(:,:,:,:)
     real*8, dimension(n_atoms, 3, 3) :: Born_effect_charges


!------------------------(5) counters------------------------------
     integer :: i_atom,j_atom, i_coord, j_coord 
     integer :: i_basis, j_basis
     integer :: i_q_point,  i_k_point,i_k_task
     integer :: i_point

     character*100 :: info_str
     character*8  :: cdate
     character*10 :: ctime
     real*8  time_start,time_end
     character(*), parameter :: deffmt = '2X'

     real*8  change_of_first_order_DM
     logical :: below_it_limit

     logical, parameter :: calculate_Born_effective_charges = .false.
!-------------------shanghui end define variables------------------------------------




!-------------------shanghui begin allocate variables------------------------------------
!----------------------(1)grid-------------------------------------------
       allocate(first_order_rho(n_full_points))
       allocate(first_order_potential(n_full_points))

       allocate(first_order_rho_Re(n_full_points))
       allocate(first_order_rho_Im(n_full_points))
       allocate(first_order_potential_Re(n_full_points))
       allocate(first_order_potential_Im(n_full_points))

       allocate(rho_free_gradient(n_full_points))
       allocate(v_free_gradient(n_full_points))

!----------------------(2)matrix-----------------------------------------
       allocate(first_order_S_sparse(n_hamiltonian_matrix_size))
       allocate(first_order_S_complex(n_basis, n_basis,n_k_points_task))
       allocate(first_order_H_sparse(n_hamiltonian_matrix_size,n_spin))
       allocate(first_order_H_complex(n_basis, n_basis,n_k_points_task))

       allocate(first_order_U_complex(n_basis,n_basis,n_k_points_task))

       allocate(density_matrix_sparse(n_hamiltonian_matrix_size))
       allocate(first_order_density_matrix_sparse(n_hamiltonian_matrix_size))
       allocate(energy_density_matrix_sparse(n_hamiltonian_matrix_size))
       allocate(first_order_energy_density_matrix_sparse(n_hamiltonian_matrix_size))

!----------------------(3) frequencies -------------------------------
       lwork=9*n_atoms
       allocate( workspace(lwork) )
       allocate( r_workspace(lwork) )
       allocate( eigenvalues(3*n_atoms) )
       allocate( mass_vector(3*n_atoms) )
       allocate( reduced_mass(3*n_atoms) )
       allocate( frequencies(3*n_atoms) )
       allocate( hellman_feynman_dynamical_matrix_free_part(3,n_atoms,3,n_atoms) )   
       allocate( hellman_feynman_dynamical_matrix_delta_part(3,n_atoms,3,n_atoms) )   
       allocate( pulay_dynamical_matrix(3,n_atoms,3,n_atoms) )   
       allocate( dynamical_matrix(3*n_atoms,3*n_atoms) )   

!------------------------(4) Born-effective-charges------------------------------
       if(calculate_Born_effective_charges) then 
       allocate(Omega_MO(n_basis,n_basis,n_k_points_task,3))
       allocate(momentum_matrix_sparse(n_hamiltonian_matrix_size_no_symmetry))
       allocate(momentum_matrix_complex(n_basis,n_basis,n_k_points_task,3))
       endif
!-------------------shanghui end allocate variables-----------------------------------




       if(packed_matrix_format.ne.PM_index) then
       call aims_stop('shanghui only use sparse matrix for DFPT_phonon_reduce_memory', & 
                      'cpscf_solver_phonon_reduce_memory')
       endif

       if(use_scalapack) then  
       !  in lapack version we will add k_weights by ourself, so only scalapck version need this.
       call kweight_occs('cpscf_solver_phonon_reduce_memory', occ_numbers) 
       endif 

!---------------shanghui begin calculate momentum_matrix---------------------------------
    if(calculate_Born_effective_charges) then 
    if(.not.use_scalapack)then !lapack version
      write(info_str,'(A)') ''
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
      write(info_str,'(A)') "==========================================================================="
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
      write(info_str,'(2X,A,1X,I4)') 'befor CPSCF, get Omega_MO = <i(k)|-r|j(k)> in MO basis'
      call localorb_info(info_str, use_unit,'(A)', OL_norm )
      write(info_str,'(A)') "=========================================================================="
      call localorb_info(info_str, use_unit,'(A)', OL_norm )

      !---------------<1> begin from momentum_matrix_complex---------------
      momentum_matrix_sparse = 0.0d0
      momentum_matrix_complex = (0.0d0,0.0d0)
      Omega_MO = (0.0d0, 0.0d0)

     do j_coord = 1,3
         call  integrate_momentum_matrix_sparse &
             ( partition_tab, l_shell_max,    &
               j_coord,    &
               momentum_matrix_sparse)
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
    endif
    endif
!---------------shanghui end calculate momentum_matrix---------------------------------





!----------- becase DFPT is working for the second derivative of E, so we use j_atom,j_coord
  do i_q_point = 1, 1 !n_k_points

     hellman_feynman_dynamical_matrix_free_part(1:3,1:n_atoms,1:3,1:n_atoms) = (0.0d0,0.0d0)
     hellman_feynman_dynamical_matrix_delta_part(1:3,1:n_atoms,1:3,1:n_atoms) = (0.0d0,0.0d0)
     pulay_dynamical_matrix(1:3,1:n_atoms,1:3,1:n_atoms) = (0.0d0,0.0d0)

     do j_atom = 1, n_atoms
     do j_coord =1 ,3

        !----begin DFPT cycle---------------
        call calculate_first_order_S_and_H_phonon_reduce_memory( & 
             i_q_point, j_atom, j_coord, & 
             first_order_S_sparse, first_order_S_complex, & 
             first_order_H_sparse, first_order_H_complex, & 
             first_order_density_matrix_sparse, first_order_U_complex, & 
             converged) 
        !----end DFPT cycle---------------


        !----begin Frequency--------------               
        call get_timestamps(time_Hessian, clock_time_Hessian)
        !------begin  delta_part for gamma point only ------------
        call  integrate_free_atom_sum_gradient_phonon_reduce_memory &
          (partition_tab, i_q_point, j_atom, j_coord, rho_free_gradient,v_free_gradient)

        call  evaluate_zero_order_DM_phonon_reduce_memory &
          (KS_eigenvector, KS_eigenvector_complex, occ_numbers,density_matrix_sparse)

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
 
             first_order_rho_Im(1:n_full_points)=  &           
             dimag(first_order_rho(1:n_full_points))  
          enddo  
 
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
             .true., multipole_radius_sq, &
             l_hartree_max_far_distance, &
             outer_potential_radius, & 
             hellman_feynman_term_at_gamma(1:3,1:n_atoms)) 
 
             do i_coord = 1,3 
				do i_atom = 1,n_atoms
					hellman_feynman_dynamical_matrix_delta_part(i_coord,i_atom,j_coord,j_atom) = & 
					dcmplx(hellman_feynman_term_at_gamma(i_coord,i_atom),0.0d0)
				enddo
             enddo
 

          do i_point =1 ,n_full_points
             first_order_potential(i_point)=  &           
             dcmplx(first_order_potential_Re(i_point),0.0d0) 
          enddo  
 
              first_order_potential(1:n_full_points)= &      !@
              first_order_potential(1:n_full_points)   &     !@
              +(-v_free_gradient(1:n_full_points))           !@

        !------end  delta_part for gamma point only ------------



        call integrate_first_order_rho_phonon_reduce_memory(partition_tab, l_shell_max,  &
             density_matrix_sparse,first_order_density_matrix_sparse, &
             i_q_point, j_atom, j_coord, &
             first_order_rho)

!-----------dynamical_matrix: (1)begin for hellman_feynman_dynamical_matrix--------------
      !------(1.1) free_part------------
        call integrate_free_atom_dynamical_matrix_phonon_reduce_memory( & 
             hellman_feynman_dynamical_matrix_free_part,i_q_point)

     !------(1.2) delta_part------------
     !   call integrate_hellam_dynamical_matrix_phonon_reduce_memory &
     !       (partition_tab, rho, first_order_rho,  &
     !        i_q_point, j_atom, j_coord, &
     !        hellman_feynman_dynamical_matrix_delta_part)

!-----------dynamical_matrix: (1)end for hellman_feynman_dynamical_matrix--------------


!-----------dynamical_matrix :(2)begin for pulay_dynamical_matrix------------------------

      !---here add DM1 in order to get the same result as DFPT_phonon_gamma, as shown in DFPT_reduced_mem,
      !---all DM1,EM1 need to be calcualted after cpscf convergence. 
      call evaluate_first_order_DM_phonon_reduce_memory(first_order_S_complex,  &
             KS_eigenvector, KS_eigenvector_complex, occ_numbers,   &
             first_order_U_complex,first_order_density_matrix_sparse)


     !----------prepare EDM(0),EDM(1)for Hessian----------
        call  evaluate_zero_order_EDM_phonon_reduce_memory &
             (KS_eigenvector, KS_eigenvector_complex,KS_eigenvalue, occ_numbers,energy_density_matrix_sparse)

        if(use_scalapack) then
        call  construct_first_order_hamiltonian_scalapack(first_order_H_sparse)  
        endif

        call  evaluate_first_order_EDM_phonon_reduce_memory( & 
              first_order_S_complex,& 
              first_order_H_complex,  &
              KS_eigenvector_complex, KS_eigenvalue, occ_numbers,  &
              first_order_U_complex, &
              first_order_energy_density_matrix_sparse)

    !-------shanghui begin parallel------
    ! if(myid.eq.0) then
    !  write(use_unit,*) '************first_order_edm(atom1,X)****************'
    !  write(use_unit,*) (first_order_energy_density_matrix_sparse(i_basis),i_basis=1,n_hamiltonian_matrix_size-1)
    !  write(use_unit,*) '************first_order_edm(atom1,X)****************'
    ! endif
    !-------shanghui end parallel------

        call integrate_dynamical_matrix_phonon_reduce_memory &
         ( hartree_potential, first_order_potential, &
           rho, rho_gradient, first_order_rho,  &
           partition_tab, l_shell_max,     &
           density_matrix_sparse, first_order_density_matrix_sparse,  &
           energy_density_matrix_sparse, first_order_energy_density_matrix_sparse,  &
           i_q_point, j_atom, j_coord, &
           hellman_feynman_dynamical_matrix_delta_part, &
           pulay_dynamical_matrix &
         )
!-----------dynamical_matrix :(2) end for pulay_dynamical_matrix------------------------

!------------------begin for debug comparing with DFPT_phonon_gamma-------------------------
!      if(myid.eq.0) then
!       write(use_unit,*) ''
!       write(use_unit,*) 'hellman_feynman_free:' 
!       do i_atom = 1 , n_atoms
!       do i_coord = 1 , 3
!         write(use_unit,*)  &
!              (hellman_feynman_dynamical_matrix_free_part(i_coord,i_atom,j_coord,j_atom))* &
!                 hartree / (bohr*bohr) 
!       enddo 
!       enddo
! 
!       write(use_unit,*) ''
!       write(use_unit,*) 'hellman_feynman_delta:' 
!       do i_atom = 1 , n_atoms
!       do i_coord = 1 , 3
!         write(use_unit,*)  &
!               (hellman_feynman_dynamical_matrix_delta_part(i_coord,i_atom,j_coord,j_atom))* & 
!                 hartree / (bohr*bohr) 
!       enddo 
!       enddo
! 
!      write(use_unit,*) ''
!      write(use_unit,*) 'pulay:' 
!      do i_atom = 1 , n_atoms
!      do i_coord = 1 , 3
!        write(use_unit,*)  &
!               (pulay_dynamical_matrix(i_coord,i_atom,j_coord,j_atom))* &
!                hartree / (bohr*bohr) 
!      enddo 
!      enddo
!       call aims_stop()
!      endif
!------------------end for debug comparing with DFPT_phonon_gamma-------------------------
 

     call get_times(time_Hessian, clock_time_Hessian, &
        &              tot_time_Hessian, tot_clock_time_Hessian)

     write(info_str,'(A)') ''
     call localorb_info(info_str, use_unit,'(A)', OL_norm  )
     call output_timeheader(deffmt, info_str, OL_norm)
     call output_times(deffmt, "Time for first_order_S of this pertubation", &
        &                 time_first_order_S, clock_time_first_order_S, OL_norm)
     call output_times(deffmt, "Time for Hessian of this pertubation", &
        &                 time_Hessian, clock_time_Hessian, OL_norm)

     write(info_str,'(A)') "==========================================================================="
     call localorb_info(info_str, use_unit,'(A)', OL_norm  )


!----------------begin Born effect charge----------------------------
     if(calculate_Born_effective_charges) then 
     call  evaluate_Born_effective_charges_phonon_reduce_memory &
          (occ_numbers, Omega_MO, first_order_U_complex, &  
           Born_effect_charges(j_atom,j_coord,1:3))
           !Born_effect_charges(j_atom,j_coord, E_field_i_coord) 
       
!            write(info_str,'(A,I5,2X,A,2X,I5)') 'DFPT for TOTAL Born effect charge: j_atom,j_coord', & 
!                                       j_atom,species_name(species(j_atom)),j_coord
!              call localorb_info(info_str, use_unit,'(A)', OL_norm)
! 
!            do i_coord=1,3
!               write(info_str,*) Born_effect_charges(i_coord,j_coord)
!               call localorb_info(info_str, use_unit,'(A)', OL_norm)
!            enddo
     endif    
!----------------end Born effect charge----------------------------
 

     enddo !j_coord =1 ,3
     enddo !j_atom = 1, n_atoms

!------------Begin print out Born effect charge----------------------------
     if(calculate_Born_effective_charges) then 
     write(info_str,'(A)') ''
     call localorb_info(info_str, use_unit,'(A)', OL_norm  )
     write(info_str,'(A)') 'DFPT for Born effective charges (electronic) in cartesian axis' 
     call localorb_info(info_str, use_unit,'(A)', OL_norm)

     do j_atom = 1, n_atoms 
        write(info_str,'(A,I5,2X,A)') 'atom', j_atom,species_name(species(j_atom))
        call localorb_info(info_str, use_unit,'(A)', OL_norm)
        
         write(info_str,'(6x,"Ex  (",3f15.5," )")') & 
            (Born_effect_charges(j_atom,j_coord,1), j_coord = 1,3)
         call localorb_info(info_str, use_unit,'(A)', OL_norm)
         write(info_str,'(6x,"Ey  (",3f15.5," )")') & 
            (Born_effect_charges(j_atom,j_coord,2), j_coord = 1,3)
         call localorb_info(info_str, use_unit,'(A)', OL_norm)
         write(info_str,'(6x,"Ez  (",3f15.5," )")') & 
            (Born_effect_charges(j_atom,j_coord,3), j_coord = 1,3)
         call localorb_info(info_str, use_unit,'(A)', OL_norm)

     enddo ! j_atom 
     endif
!------------End print out Born effect charge----------------------------


!----------dynamical_matrix: (3) begin for prepare total dynamical_matrix-------------
    if(myid.eq.0) then
      write(use_unit,*) ''
      write(use_unit,*) 'Total_hessian in cpscf.f90:'

      write(use_unit,*) ''
      write(use_unit,*) 'hellman_feynman_hessian_free_part in cpscf.f90:'
      do i_atom = 1 , n_atoms
      do i_coord = 1 , 3

         write(use_unit,'(100E20.12)')  &
             (((dble(hellman_feynman_dynamical_matrix_free_part(i_coord,i_atom,j_coord,j_atom)))* &
                 hartree / (bohr*bohr) , j_coord=1,3) ,j_atom=1,n_atoms)
      enddo
      enddo

      write(use_unit,*) ''
      write(use_unit,*) 'hellman_feynman_hessian_delta_part in cpscf.f90:'
      do i_atom = 1 , n_atoms
      do i_coord = 1 , 3

         write(use_unit,'(100E20.12)')  &
             (((dble(hellman_feynman_dynamical_matrix_delta_part(i_coord,i_atom,j_coord,j_atom)))* &
                 hartree / (bohr*bohr) , j_coord=1,3) ,j_atom=1,n_atoms)
      enddo
      enddo
      write(use_unit,*) ''
      write(use_unit,*) 'pulay_hessian in cpscf.f90:'
      do i_atom = 1 , n_atoms
      do i_coord = 1 , 3

         write(use_unit,'(100E20.12)')  &
             (((dble(pulay_dynamical_matrix(i_coord,i_atom,j_coord,j_atom)))* &
                 hartree / (bohr*bohr) , j_coord=1,3) ,j_atom=1,n_atoms)
      enddo
      enddo
    endif








      do i_atom = 1 , n_atoms
      do i_coord = 1 , 3

        do j_atom=1, n_atoms
        do j_coord =1, 3
        dynamical_matrix(3*i_atom+i_coord-3,3*j_atom+j_coord-3)= &
        ( hellman_feynman_dynamical_matrix_free_part(i_coord,i_atom,j_coord,j_atom)+ &
          hellman_feynman_dynamical_matrix_delta_part(i_coord,i_atom,j_coord,j_atom)+ &
          pulay_dynamical_matrix(i_coord,i_atom,j_coord,j_atom))* &
                 hartree / (bohr*bohr)
        enddo
        enddo

      enddo
      enddo

    ! symmetrize dynamical matrix 
    do i_coord = 1, 3*n_atoms
    do j_coord = 1, i_coord - 1
       buf = (dynamical_matrix(i_coord,j_coord) + dconjg(dynamical_matrix(j_coord,i_coord)) )/2d0
       dynamical_matrix(i_coord,j_coord) = buf
       dynamical_matrix(j_coord,i_coord) = dconjg(buf)
    end do
    end do

!------------------- shanghui make asr 2013.10.21----------------------

        !-------shanghui begin parallel------
        if(myid.eq.0) then
          write(use_unit,*) ''
          write(use_unit,*) 'Make Acoustic sum rule in cpscf.f90:'
        endif 
        !-------shanghui end parallel--------
    do i_coord = 1 , 3
    do j_coord = 1 , 3

        do i_atom = 1 , n_atoms
           dynamical_matrix(3*i_atom+i_coord-3,3*i_atom+j_coord-3)=(0.0d0,0.0d0)
        do j_atom=1, n_atoms

        if(j_atom.ne.i_atom) then
        dynamical_matrix(3*i_atom+i_coord-3,3*i_atom+j_coord-3)=  &
        dynamical_matrix(3*i_atom+i_coord-3,3*i_atom+j_coord-3)-  &
        dynamical_matrix(3*i_atom+i_coord-3,3*j_atom+j_coord-3)
        endif

        enddo
        enddo

    enddo
    enddo

    if(myid.eq.0) then
      do i_atom = 1 , n_atoms 
      do i_coord = 1 , 3

         write(use_unit,'(100E20.12)')  & 
             ((dble(dynamical_matrix(3*i_atom+i_coord-3,3*j_atom+j_coord-3)), j_coord=1,3) ,j_atom=1,n_atoms)
      enddo
      enddo
    endif
!----------dynamical_matrix: (3) end for prepare total dynamical_matrix-------------

	if (use_vdw_correction_hirshfeld) then
		call calc_vdw_hirshfeld(vdw_energy, vdw_forces, vdW_Hessian)
		dynamical_matrix = dynamical_matrix + vdW_Hessian
	endif

!----------dynamical_matrix: (4) begin for diagonal dynamical_matrix----------------
    if(myid.eq.0) then
    write(use_unit,*) ''
    write(use_unit,*) 'Get frequencies in cpscf.f90:'
    endif

    hessian_factor   = const_eV/(const_u*const_Angstr*const_Angstr)
    dynamical_matrix = dynamical_matrix * hessian_factor

    do i_atom=1, n_atoms 
       do i_coord = 1, 3
        mass_vector(3*i_atom+i_coord-3) = 1.0d0/dsqrt(species_m(species(i_atom)))
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

!----------dynamical_matrix: (4) end for diagonal dynamical_matrix----------------


!----------dynamical_matrix: (5) begin for print frequencies----------------
   if(myid.eq.0) then 
    write(use_unit,*) 'DFPT-Results: '
    write(use_unit,*)
    write(use_unit,*) 'List of all frequencies found:'
    write(use_unit,'(A13,A25)') 'Mode number','Frequency [cm^(-1)]'

    do i_coord = 1, 3*n_atoms, 1
       if (eigenvalues(i_coord).gt.0d0) then
          frequencies(i_coord) = sqrt(eigenvalues(i_coord))
       else if (eigenvalues(i_coord).lt.0d0) then
          frequencies(i_coord) = -sqrt(-eigenvalues(i_coord))
       else 
          frequencies(i_coord) = 0d0
       end if

       write(use_unit,'(I13,F25.8)') i_coord, frequencies(i_coord)/(200*pi*const_c)
    enddo
   endif 
!----------dynamical_matrix: (5) end for print frequencies----------------


        !----end Frequency--------------               



  enddo    !i_q_point = 1, 1 !n_k_points
    



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
       deallocate(first_order_S_sparse)
       deallocate(first_order_S_complex)
       deallocate(first_order_H_sparse)
       deallocate(first_order_H_complex)

       deallocate(first_order_U_complex)

       deallocate(density_matrix_sparse)
       deallocate(first_order_density_matrix_sparse)
       deallocate(energy_density_matrix_sparse)
       deallocate(first_order_energy_density_matrix_sparse)

!----------------------(3) frequencies--------------------------------
       deallocate( workspace )
       deallocate( r_workspace )
       deallocate( eigenvalues )
       deallocate( mass_vector )
       deallocate( reduced_mass )
       deallocate( frequencies )
       deallocate( hellman_feynman_dynamical_matrix_free_part )
       deallocate( hellman_feynman_dynamical_matrix_delta_part )
       deallocate( pulay_dynamical_matrix )
       deallocate( dynamical_matrix )

!------------------------(4) Born-effective-charges------------------------------
       if(calculate_Born_effective_charges) then 
       deallocate(Omega_MO)
       deallocate(momentum_matrix_sparse)
       deallocate(momentum_matrix_complex)
       endif
!-------------------shanghui end deallocate------------------------------------





    end subroutine cpscf_solver_phonon_reduce_memory
!******
