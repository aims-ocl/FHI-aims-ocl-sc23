!****s* FHI-aims/cpscf_solver_p0
!  NAME
!    cpscf_solver
!  SYNOPSIS

    subroutine cpscf_solver_phonon_p0 &

    (converged)

!  PURPOSE
!  an cpscf process for phonon_gamma 
!  shanghui 2013.12.30
!  
!  USES

      use constants, only: bohr, hartree, pi
      use dimensions
      use timing
      use physics
      use species_data
      use localorb_io
      use geometry
      use mpi_tasks, only: myid, n_tasks
      use runtime_choices, only: sc_iter_limit
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
     real*8, allocatable :: workspace(:)
     real*8, allocatable :: eigenvalues(:)
     real*8, allocatable :: mass_vector(:)
     real*8, allocatable :: reduced_mass(:)
     real*8, allocatable :: frequencies(:)
     real*8, allocatable :: hessian(:,:)     
     real*8   norm,buf,hessian_factor
     real*8, external :: ddot
!---------------shanghui end for vib-cal----------------

      real*8, parameter:: alpha=0.5d0


!-------------------shanghui add for DFPT------------------------------------
!------------------------(1) grid ---------------------------------------
     real*8, allocatable :: first_order_rho(:,:,:)
     real*8, allocatable :: first_order_potential(:,:,:)

     real*8, allocatable :: rho_free_gradient(:,:,:)
     real*8, allocatable :: v_free_gradient(:,:,:)

!------------------------(2) matrix -----------------------------------
     real*8, allocatable :: first_order_S(:,:,:,:)
     complex*16, allocatable :: first_order_S_complex(:,:,:,:,:)
     real*8, allocatable :: first_order_H(:,:,:,:)
     complex*16, allocatable :: first_order_H_complex(:,:,:,:,:) 

     real*8,allocatable :: first_order_U(:,:,:,:)
     complex*16,allocatable :: first_order_U_complex(:,:, :,:,:)
     real*8, allocatable ::  first_order_E(:, :,:)

     real*8, allocatable ::  density_matrix(:,:)
     real*8, allocatable ::  first_order_density_matrix(:, :, :,:)
     real*8, allocatable ::  old_first_order_density_matrix(:,:,:,:)

     real*8, allocatable ::  energy_density_matrix(:,:) 
     real*8, allocatable ::  first_order_energy_density_matrix(:,:,:,:)

!------------------shanghui end add for DFPT------------------------------------- 





      ! counters

      integer :: i_spin
      integer :: i_atom,j_atom, i_coord, j_coord,i_coord_2
      integer :: i_basis, j_basis,  i_point,i_state, i_k_point,i_k_task


      character(*), parameter :: func = 'cpscf_solver'
      
      logical first_iteration
     !shanghui------------------------------------------------------------------------

 
      real*8  change_of_first_order_DM
      real*8  time_start,time_end
      real*8  temp1, temp2
      integer  max_occ_number(n_spin)

      real*8  pulay_hessian(3,n_atoms,3,n_atoms)
      real*8  hellman_feynman_hessian(3,n_atoms,3,n_atoms)
      real*8  hellman_feynman_hessian_free_part(3,n_atoms,3,n_atoms)
      real*8  hellman_feynman_hessian_delta_part(3,n_atoms,3,n_atoms)
      real*8  hellman_feynman_hessian_temp(3,n_atoms)

     !shanghui------------------------------------------------------------------------



!-------------------shanghui begin add for DFPT------------------------------------
!----------------------(1)grid-------------------------------------------
       allocate(first_order_rho(3,n_atoms,n_full_points))
       allocate(first_order_potential(3,n_atoms,n_full_points))

       allocate(rho_free_gradient(3,n_atoms,n_full_points))
       allocate(v_free_gradient(3,n_atoms,n_full_points))

!----------------------(2)matrix-----------------------------------------
       allocate(first_order_S(3, n_atoms, n_Cbasis,n_Cbasis))
       allocate(first_order_S_complex(3,n_atoms, n_basis, n_basis,n_k_points_task))
       allocate(first_order_H(3, n_atoms, n_Cbasis,n_Cbasis))
       allocate(first_order_H_complex(3,n_atoms, n_basis, n_basis,n_k_points_task))

       allocate(first_order_U(3, n_atoms, n_basis,n_basis))
       allocate(first_order_U_complex(3, n_atoms, n_basis,n_basis,n_k_points_task))
       allocate(first_order_E(3, n_atoms, n_basis))

       allocate(density_matrix(n_Cbasis,n_Cbasis))
       allocate(first_order_density_matrix(3, n_atoms, n_Cbasis,n_Cbasis))
       allocate(old_first_order_density_matrix(3, n_atoms, n_Cbasis,n_Cbasis))
       allocate(energy_density_matrix(n_Cbasis,n_Cbasis))
       allocate(first_order_energy_density_matrix(3,n_atoms,n_Cbasis,n_Cbasis))
!-------------------shanghui end add for DFPT-----------------------------------


!---------------shanghui begin for vib-cal-------------
       lwork=9*n_atoms
       allocate( workspace(lwork) ) 
       allocate( eigenvalues(3*n_atoms) )
       allocate( mass_vector(3*n_atoms) )
       allocate( reduced_mass(3*n_atoms) )
       allocate( frequencies(3*n_atoms) )
       allocate( hessian(3*n_atoms,3*n_atoms) )   
!---------------shanghui end  for vib-cal-------------




      ! begin work
      number_of_loops = 0

      below_it_limit = (number_of_loops.lt.sc_iter_limit)
     
      first_iteration = .true.
  

   !--------------shanghui test first_order_S-------------------
     first_order_S=0.0d0 
     first_order_S_complex=(0.0d0,0.0d0) 
     call integrate_first_order_S_p0(partition_tab, l_shell_max, first_order_S)

     call  construct_first_order_S_p0(first_order_S,first_order_S_complex)

   !-------shanghui begin parallel------
    i_k_task = 0

    do i_k_point = 1,2  !n_k_points, 1
    if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
 
       i_k_task = i_k_task + 1

      write(use_unit,*) '************shanghui begain first_order_S(atom1,X)****************'
      write(use_unit,*) 'k=',i_k_point 
      do i_basis=1,n_basis
      write(use_unit,'(6f15.9)') (first_order_S_complex(1,1,i_basis,j_basis,i_k_task),j_basis=1,n_basis )
      enddo
      write(use_unit,*) '************shanghui end first_order_S****************'
      write(use_unit,*) ''

    endif
    enddo 
   !-------shanghui end parallel------

   !--------------shanghui end test first_order_S---------------   

!     find the max_occ_number
        do i_spin = 1, n_spin, 1
            max_occ_number(i_spin) = 0
            do i_state = n_states, 1, -1
             if (dabs(occ_numbers(i_state,i_spin,1)).gt.0.d0) then
              max_occ_number(i_spin) = i_state
              exit
             endif
            enddo
        enddo

     
    !-------shanghui begin parallel------
    if(myid.eq.0) then
     write(use_unit,*) '-----------------------------------------------'
     write(use_unit,*) 'shanghui test n_basis,max_occ:',n_basis,max_occ_number(1)
     write(use_unit,*) 'shanghui test n_max_compute_ham',n_max_compute_ham
     write(use_unit,*) 'shanghui test n_occ_atoms',n_occ_atoms
     write(use_unit,*) '-----------------------------------------------'
    endif
    !-------shanghui end parallel------
    


    density_matrix=0.0d0
    first_order_density_matrix=0.0d0
    old_first_order_density_matrix=0.0d0

    first_order_H=0.0d0
    first_order_H_complex=(0.0d0,0.0d0)
    first_order_U=0.0d0 
    first_order_U_complex=(0.0d0,0.0d0)
    first_order_E=0.0d0
    v_free_gradient=0.0d0
    rho_free_gradient=0.0d0
!---------------------begin calculate gradient_V_MP-------------------------------
        call  integrate_free_atom_sum_gradient_p0 &
              (partition_tab, rho_free_gradient,v_free_gradient) 

!        call update_hartree_potential_p1 &
!            ( hartree_partition_tab, free_rho_superpos, &
!            rho, rho_multipole_old, &
!            delta_v_hartree_part_at_zero, &
!            delta_v_hartree_deriv_l0_at_zero, &
!            multipole_moments, multipole_radius_sq, &
!            l_hartree_max_far_distance, rho_radial_multipole_old, &
!            outer_potential_radius &
!            )
!
!
!        call sum_up_whole_potential_p1 &
!                ( delta_v_hartree_part_at_zero, &
!                delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
!                partition_tab, rho, &
!                free_hartree_superpos, free_rho_superpos,  &
!                hartree_potential,  &
!                hartree_delta_energy, en_elec_delta, hartree_multipole_correction, &
!                pot_ion_embed, en_density_embed, &
!                multipole_forces, .true., multipole_radius_sq, &
!                l_hartree_max_far_distance, hellman_feynman_forces, &
!                energy_deriv_tress_tensor, rho_multipole_old, &
!                outer_potential_radius, v_free_gradient,rho_free_gradient)

!---------------------end calculate gradient_V_MP-------------------------------

      !----this part is used to make old_first_order_density_matrix-----------
      !----for just check evaluate_first_order_DM, we comment it--------------

     call evaluate_first_order_DM_p0(first_order_S_complex,  &
             KS_eigenvector_complex, KS_eigenvalue, occ_numbers,   &
             first_order_U_complex,density_matrix,old_first_order_density_matrix)

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



!---------------(1) Begin  update first_order_rho-----------------------------------
       call cpu_time(time_start)

       call evaluate_first_order_DM_p0(first_order_S_complex,  &
             KS_eigenvector_complex, KS_eigenvalue, occ_numbers,   &
             first_order_U_complex,density_matrix,first_order_density_matrix)

!         if(myid.eq.0) then
!        write(use_unit,*) '************shanghui begain first_order_dm(atom1,X)****************'
!         do i_basis=1,n_Cbasis
!         write(use_unit,'(6f15.9)') (first_order_density_matrix(1,1,i_basis,j_basis)   & 
!                              ,j_basis=1,n_Cbasis )
!         enddo
!        write(use_unit,*) '************shanghui enddo first_order_dm(atom1,X)****************'
!         endif

         !shanghui: we need to perform density matrix mixing here
         
          change_of_first_order_DM =0.0d0         

         do i_basis = 1,n_Cbasis
         do j_basis = 1,n_Cbasis

            do i_atom = 1,n_atoms
 
            do i_coord = 1,3
         change_of_first_order_DM =                &
         max( change_of_first_order_DM,             &
         dabs(first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)  &
         - old_first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)) )
         
         first_order_density_matrix(i_coord,i_atom,i_basis,j_basis) =       &
         (1.0d0-alpha)*old_first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)+  &
         alpha*first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)
     
         old_first_order_density_matrix(i_coord,i_atom,i_basis,j_basis) =   &
         first_order_density_matrix(i_coord,i_atom,i_basis,j_basis)
            enddo

            enddo
       
         enddo
         enddo

        !-------shanghui begin parallel------
        if(myid.eq.0) then
        write(use_unit,*) "((((((((((((((((((((((((((((((((("
        write(use_unit,*) change_of_first_order_DM
        write(use_unit,*) ")))))))))))))))))))))))))))))))))"
        endif 
        !-------shanghui end parallel--------


!-------------(1) end first-order-density update and mixing--------
        first_iteration = .false.

!------------(2)begain to calculate first_order_H-----------------
        !----in pbc case, we have used density matrix instead of orbitals, it looks neat. 
        call integrate_first_order_rho_p0(partition_tab, l_shell_max,  &
             density_matrix,first_order_density_matrix, &
             first_order_rho)


        do i_atom=1,n_atoms
           do i_coord=1,3


              first_order_rho(i_coord,i_atom,1:n_full_points)=  &           !@ 
              first_order_rho(i_coord,i_atom,1:n_full_points)+  &           !@
              rho_free_gradient(i_coord,i_atom,1:n_full_points)             !@

        call update_hartree_potential_shanghui_p0 &
            ( hartree_partition_tab,first_order_rho(i_coord,i_atom,1:n_full_points),& 
            delta_v_hartree_part_at_zero, &
            delta_v_hartree_deriv_l0_at_zero, &
            multipole_moments, multipole_radius_sq, &
            l_hartree_max_far_distance, &
            outer_potential_radius )


        !-------shanghui begin parallel------
        if(myid.eq.0) then
       write(use_unit,*) '{----shanghui in cpscf_solver.f:-----for first_order_hartree_potential:'
       write(use_unit,*) 'i_coord,i_atom',i_coord,i_atom,forces_on
       write(use_unit,*) '----shanghui in cpscf_solver.f:-----for first_order_hartree_potential}'
        endif 
        !-------shanghui end parallel--------


        call sum_up_whole_potential_shanghui_p0 &
                ( delta_v_hartree_part_at_zero, &
                delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
                partition_tab, first_order_rho(i_coord,i_atom,1:n_full_points), &
                first_order_potential(i_coord,i_atom,1:n_full_points),  & !<--------get first_order_DM_potential
                .false., multipole_radius_sq, &
                l_hartree_max_far_distance, &
                outer_potential_radius, & 
                hellman_feynman_hessian_temp) 

 
                first_order_potential(i_coord,i_atom, 1:n_full_points)= &    !@
                first_order_potential(i_coord,i_atom,1:n_full_points)   &     !@
                +(-v_free_gradient(i_coord,i_atom,1:n_full_points))       !@
 


           enddo
        enddo


     call  integrate_first_order_H_p0 &
           (hartree_potential,first_order_potential, rho, rho_gradient,&
           partition_tab, l_shell_max, en_xc, en_pot_xc,    &
           density_matrix,first_order_density_matrix, first_order_H &
           )
        
      i_k_task = 0
      
      do i_k_point = 1,n_k_points, 1
      if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
      
         i_k_task = i_k_task + 1   
       
      call  construct_first_order_H_p0(first_order_H,   & 
            first_order_H_complex(:,:,:,:,i_k_task),i_k_point)

       if(i_k_point.eq.1) then 

        !-------shanghui begin parallel------
       write(use_unit,*) '************shanghui begain first_order_H(atom1,X)****************'
       write(use_unit,*) 'k=',i_k_point
       do i_basis=1,n_basis
        write(use_unit,'(6f15.9)') (first_order_H_complex(1,1,i_basis,j_basis,i_k_task),j_basis=1,n_basis )
       enddo
       write(use_unit,*) '************shanghui end first_order_H(atom1,X)****************'
       write(use_unit,*) ''
        !-------shanghui end parallel--------

       endif 
  

     call evaluate_first_order_U_p0(first_order_H_complex(:,:,:,:,i_k_task),& 
          first_order_S_complex(:,:,:,:,i_k_task),  &
          KS_eigenvector_complex(:,:,:,i_k_task), KS_eigenvalue(:,:,i_k_point), &  
          occ_numbers(:,:,i_k_point),  &
          first_order_U_complex(:,:,:,:,i_k_task),first_order_E)

     endif
     enddo ! n_k_point


        call cpu_time(time_end)

        if(myid.eq.0) then 
        write(use_unit,*) '@@@@@@@@time for CPSCF@@@@@@@@@@@@@@@@@@',time_end-time_start
        endif


!--------------(3) end solve first_order_U problem----------------


! --------- Check convergence ----->>


!         check convergence of self-consistency loop
         

        converged = (change_of_first_order_DM.lt.1.0d-4).and.(number_of_loops.ne.1)

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

        ! current SCF loop ends here

! ----- Printing out time data -- >>


! << ---- end printing out data---------


!       this is the end of the self-consistency cycle.

  end do SCF_LOOP
! << ------ end self consistent cycle--------

      total_number_of_loops = total_number_of_loops + number_of_loops



!-------------(1) for hellman_feynman_hessian--------------

         hellman_feynman_hessian_free_part=0.0d0
         hellman_feynman_hessian_delta_part=0.0d0

        !------(1.1) free_part------------
        call integrate_free_atom_hessian_p0 (hellman_feynman_hessian_free_part)

 

        !------(1.2) delta_part------------
        call integrate_first_order_rho_p0(partition_tab, l_shell_max,  &
             density_matrix,first_order_density_matrix, &
             first_order_rho)

        do i_atom=1,n_atoms
           do i_coord=1,3
 
              first_order_rho(i_coord,i_atom,1:n_full_points)=  &           !@ 
              first_order_rho(i_coord,i_atom,1:n_full_points)+  &           !@
              rho_free_gradient(i_coord,i_atom,1:n_full_points)             !@
 
        call update_hartree_potential_shanghui_p0 &
            ( hartree_partition_tab,first_order_rho(i_coord,i_atom,1:n_full_points),& 
            delta_v_hartree_part_at_zero, &
            delta_v_hartree_deriv_l0_at_zero, &
            multipole_moments, multipole_radius_sq, &
            l_hartree_max_far_distance, &
            outer_potential_radius )
 
 
        !-------shanghui begin parallel------
        if(myid.eq.0) then
       write(use_unit,*) '{----shanghui for hellman_feynman_hessian_delta_part:'
       write(use_unit,*) 'i_coord,i_atom',i_coord,i_atom,forces_on
       write(use_unit,*) '{----shanghui end for hellman_feynman_hessian_delta_part:'
        endif 
        !-------shanghui end parallel--------

 
        call sum_up_whole_potential_shanghui_p0 &
                ( delta_v_hartree_part_at_zero, &
                delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
                partition_tab, first_order_rho(i_coord,i_atom,1:n_full_points), &
                first_order_potential(i_coord,i_atom,1:n_full_points),  & !<--------get first_order_DM_potential
                .true., multipole_radius_sq, &
                l_hartree_max_far_distance, &
                outer_potential_radius,&
                hellman_feynman_hessian_delta_part(1:3,1:n_atoms,i_coord,i_atom) ) 
 
 
                first_order_potential(i_coord,i_atom, 1:n_full_points)= &    !@
                first_order_potential(i_coord,i_atom,1:n_full_points)   &     !@
                +(-v_free_gradient(i_coord,i_atom,1:n_full_points))       !@

           enddo
        enddo
 

       do i_atom = 1, n_atoms
       do i_coord = 1, 3
           
          do j_atom = 1, n_atoms
          do j_coord = 1 , 3 

               hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) =  & 
               hellman_feynman_hessian_free_part(i_coord,i_atom,j_coord,j_atom) + & 
               hellman_feynman_hessian_delta_part(i_coord,i_atom,j_coord,j_atom)

          enddo  
          enddo 
       
       enddo 
       enddo  



!-------------(2) for pulay_hessian------------------------
         density_matrix=0.0d0
         first_order_density_matrix=0.0d0
         energy_density_matrix=0.0d0
         first_order_energy_density_matrix=0.0d0

     !----------prepare EDM(0),EDM(1),DM(0),DM(1) for Hessian----------
     !----------here is the PBC code-------->
   call  evaluate_first_zero_order_DM_EDM_p0( first_order_S_complex,& 
         first_order_H_complex,  &
         KS_eigenvector_complex, KS_eigenvalue, occ_numbers,  &
         first_order_U_complex,density_matrix,first_order_density_matrix, &
         energy_density_matrix,first_order_energy_density_matrix)

   call integrate_pulay_hessian_p0 &
     ( hartree_potential, first_order_potential, &
       rho, rho_gradient,  &
       partition_tab, l_shell_max,     &
       KS_eigenvector(:,:,1,1),       &
       density_matrix, first_order_density_matrix,  &
       energy_density_matrix, first_order_energy_density_matrix,  &
       max_occ_number, &
       pulay_hessian &
     )


!---------------shanghui begin for vib-cal-------------

    !-------shanghui begin parallel------
    if(myid.eq.0) then
      write(use_unit,*) ''
      write(use_unit,*) 'Total_hessian in cpscf.f90:'
      
      write(use_unit,*) ''     
      write(use_unit,*) 'hellman_feynman_hessian_free_part in cpscf.f90:'
      do i_atom = 1 , n_atoms
      do i_coord = 1 , 3

         write(use_unit,'(100E20.12)')  &
             (((hellman_feynman_hessian_free_part(i_coord,i_atom,j_coord,j_atom))* &
                 hartree / (bohr*bohr) , j_coord=1,3) ,j_atom=1,n_atoms)
      enddo 
      enddo

      write(use_unit,*) ''     
      write(use_unit,*) 'hellman_feynman_hessian_delta_part in cpscf.f90:'
      do i_atom = 1 , n_atoms
      do i_coord = 1 , 3

         write(use_unit,'(100E20.12)')  &
             (((hellman_feynman_hessian_delta_part(i_coord,i_atom,j_coord,j_atom))* &
                 hartree / (bohr*bohr) , j_coord=1,3) ,j_atom=1,n_atoms)
      enddo 
      enddo
      write(use_unit,*) ''     
      write(use_unit,*) 'pulay_hessian in cpscf.f90:'
      do i_atom = 1 , n_atoms
      do i_coord = 1 , 3

         write(use_unit,'(100E20.12)')  &
             (((pulay_hessian(i_coord,i_atom,j_coord,j_atom))* &
                 hartree / (bohr*bohr) , j_coord=1,3) ,j_atom=1,n_atoms)
      enddo 
      enddo

    endif

        !-------shanghui end parallel--------


      do i_atom = 1 , n_atoms
      do i_coord = 1 , 3

!        write(use_unit,'(100E20.12)')  &
!            (((hellman_feynman_hessian_free_part(i_coord,i_atom,j_coord,j_atom)+ &
!               hellman_feynman_hessian_delta_part(i_coord,i_atom,j_coord,j_atom)+ &
!              pulay_hessian(i_coord,i_atom,j_coord,j_atom))* &
!                hartree / (bohr*bohr) , j_coord=1,3) ,j_atom=1,n_atoms)


        do j_atom=1, n_atoms
        do j_coord =1, 3
        hessian(3*i_atom+i_coord-3,3*j_atom+j_coord-3)= &
        ( hellman_feynman_hessian_free_part(i_coord,i_atom,j_coord,j_atom)+ &
          hellman_feynman_hessian_delta_part(i_coord,i_atom,j_coord,j_atom)+ &
          pulay_hessian(i_coord,i_atom,j_coord,j_atom))* &
                 hartree / (bohr*bohr)
        enddo
        enddo

      enddo
      enddo



    do i_coord = 1, 3*n_atoms
       do i_coord_2 = 1, i_coord - 1 
          buf = (hessian(i_coord_2,i_coord)+hessian(i_coord,i_coord_2))/2d0
          hessian(i_coord_2,i_coord) = buf
          hessian(i_coord,i_coord_2) = buf
       end do
    end do

!--------------------shanghui make asr 2013.10.21----------------------

        !-------shanghui begin parallel------
        if(myid.eq.0) then
      write(use_unit,*) ''
      write(use_unit,*) 'Make Acoustic sum rule in cpscf.f90:'
        endif 
        !-------shanghui end parallel--------


      do i_coord = 1 , 3
      do j_coord =1, 3 
 
        do i_atom = 1 , n_atoms 
        hessian(3*i_atom+i_coord-3,3*i_atom+j_coord-3)=0.0d0
  
        do j_atom=1, n_atoms

        if(j_atom.ne.i_atom) then         
        hessian(3*i_atom+i_coord-3,3*i_atom+j_coord-3)=  & 
        hessian(3*i_atom+i_coord-3,3*i_atom+j_coord-3)-  &
        hessian(3*i_atom+i_coord-3,3*j_atom+j_coord-3)
        endif


        enddo 
        enddo 

      enddo
      enddo



 
!--------------------shanghui end make asr 2013.10.21----------------------


        !-------shanghui begin parallel------
        if(myid.eq.0) then

      do i_atom = 1 , n_atoms 
      do i_coord = 1 , 3

         write(use_unit,'(100E20.12)')  & 
             ((hessian(3*i_atom+i_coord-3,3*j_atom+j_coord-3) , j_coord=1,3) ,j_atom=1,n_atoms)
      enddo
      enddo

    write(use_unit,*) ''
    write(use_unit,*) 'Get frequencies in cpscf.f90:'
    hessian_factor   = const_eV/(const_u*const_Angstr*const_Angstr)


      do i_atom=1, n_atoms 
        do i_coord = 1, 3
        mass_vector(3*i_atom+i_coord-3) = 1.0d0/dsqrt(species_m(species(i_atom)))
        end do
      end do

    do i_coord = 1, 3*n_atoms, 1   
       hessian(:,i_coord) = hessian(:,i_coord)*mass_vector(i_coord)       
       hessian(i_coord,:) = hessian(i_coord,:)*hessian_factor*mass_vector(i_coord) 
   end do
    




    write(use_unit,*) 'Solving eigenvalue system for Hessian Matrix'
    call DSYEV('V','U',3*n_atoms,hessian,3*n_atoms,eigenvalues,workspace,lwork,info)
    write(use_unit,*) 'Done ... '
    write(use_unit,*)
    ! calculate the eigenvectors in cartesian coordinates ?
    do i_coord = 1, 3*n_atoms
       hessian(:,i_coord) = hessian(:,i_coord)*mass_vector(:)
    end do
    
    ! Renormalize eigenvectors for output - norm has units of 1/sqrt(mass) though 
    do i_coord = 1, 3*n_atoms
      reduced_mass(i_coord) = ddot(3*n_atoms, hessian(:, i_coord), 1, hessian(:, i_coord), 1)
      norm = sqrt(reduced_mass(i_coord))
      hessian(:,i_coord) = hessian(:,i_coord)/norm
    end do

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
        !-------shanghui end parallel--------

!---------------shanghui end for vib-cal-------------




!-------------------shanghui begin add for DFPT------------------------------------
!----------------------(1) grid------------------------------------
       deallocate(first_order_rho)
       deallocate(first_order_potential)

       deallocate(rho_free_gradient)
       deallocate(v_free_gradient)

!----------------------(2) matrix-----------------------------------------
       deallocate(first_order_S)
       deallocate(first_order_S_complex)
       deallocate(first_order_H)
       deallocate(first_order_H_complex)

       deallocate(first_order_U)
       deallocate(first_order_U_complex)
       deallocate(first_order_E)

       deallocate(density_matrix)
       deallocate(first_order_density_matrix)
       deallocate(old_first_order_density_matrix)
       deallocate(energy_density_matrix)
       deallocate(first_order_energy_density_matrix)
!-------------------shanghui end add for DFPT------------------------------------




!---------------shanghui begin for vib-cal-------------
       deallocate( workspace )
       deallocate( eigenvalues )
       deallocate( mass_vector )
       deallocate( reduced_mass )
       deallocate( frequencies )
       deallocate( hessian )
!---------------shanghui end for vib-cal-------------


    end subroutine cpscf_solver_phonon_p0
!******
