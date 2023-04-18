!****s* FHI-aims/cpscf_solver
!  NAME
!    cpscf_solver
!  SYNOPSIS

    subroutine cpscf_solver_reduce_memory &
    (converged)

!  PURPOSE
!  an cpscf process.
!  USES

      use constants, only: pi
      use dimensions
      use timing
      use physics
      use species_data
      use localorb_io
      use geometry
      use debugmanager, only: module_is_debugged, debugprint
      use runtime_choices, only : DFPT_sc_accuracy_dm, DFPT_mixing, &
          sc_iter_limit
      use mpi_tasks, only: myid
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


     !----shanghui begin for vib-cal-------------
     real*8, parameter:: const_u          = 1.66053886d-27         
     real*8, parameter:: const_eV         = 1.60217653d-19
     real*8, parameter:: const_c          = 299792458d0
     real*8, parameter:: const_Angstr     = 1d-10
     real*8, parameter:: gradient_dipole_factor = 4.80320577161 ! bring grad_dipole to D/Ang
     
     integer  lwork, info  
     real*8, allocatable :: workspace(:)
     real*8, allocatable :: eigenvalues(:)
     real*8, allocatable :: mass_vector(:)
     real*8, allocatable :: reduced_mass(:)
     real*8, allocatable :: frequencies(:)
     real*8, allocatable :: hessian(:,:)   
     real*8, allocatable :: gradient_dipole(:,:)
     real*8, allocatable :: gradient_dipole_internal(:,:) 
     real*8, allocatable :: infrared_intensity(:) 
     real*8   norm,buf,hessian_factor
     !----shanghui end for vib-cal----------------



!      real*8, parameter:: alpha=0.2d0


  !-------------------(1) define first-order------------------------------------
     real*8, allocatable :: first_order_density_matrix(:,:)
     real*8, allocatable :: first_order_energy_density_matrix(:,:)
     real*8, allocatable :: old_first_order_density_matrix(:,:)

     real*8, allocatable :: first_order_rho(:)
     real*8, allocatable :: first_order_potential(:)

     real*8, allocatable :: first_order_S(:,:)
     real*8, allocatable :: first_order_U(:,:)
     real*8, allocatable :: first_order_E(:)
     real*8, allocatable :: first_order_H(:,:)

     real*8, allocatable :: density_matrix(:,:) 
     real*8, allocatable :: energy_density_matrix(:,:) 
     real*8  change_of_first_order_DM
     real*8  time_start,time_end
     integer  max_occ_number(n_spin)


      ! counters

      integer :: i_spin
      integer :: j_atom, j_coord
      integer :: i_basis, j_basis 
      integer :: i_atom,  i_coord,i_coord_2, i_point,i_state


      character(*), parameter :: func = 'cpscf_solver'
      
      real*8, external :: ddot

     !-------------(2) begin allocate-----------------------------------------------------------

       allocate(first_order_density_matrix(n_basis,n_basis))
       allocate(first_order_energy_density_matrix( n_basis,n_basis))
       allocate(old_first_order_density_matrix( n_basis,n_basis))
       allocate(first_order_rho(n_full_points))
       allocate(first_order_potential(n_full_points))
      
       allocate(first_order_S( n_basis,n_basis))
       allocate(first_order_U( n_basis,n_basis))
       allocate(first_order_E( n_basis))
       allocate(first_order_H( n_basis,n_basis))

       allocate(density_matrix(n_basis, n_basis)) 
       allocate(energy_density_matrix(n_basis, n_basis)) 
 


       !----shanghui begin for vib-cal-------------
       lwork=9*n_atoms
       allocate( workspace(lwork) ) 
       allocate( eigenvalues(3*n_atoms) )
       allocate( mass_vector(3*n_atoms) )
       allocate( reduced_mass(3*n_atoms) )
       allocate( frequencies(3*n_atoms) )
       allocate( hessian(3*n_atoms,3*n_atoms) )   
       allocate( gradient_dipole(3*n_atoms,3) )   
       allocate( gradient_dipole_internal(3, n_atoms*3))
       allocate(infrared_intensity(n_atoms*3))
       !----shanghui end for vib-cal-------------



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

    call evaluate_zero_order_DM_reduce_memory( & 
         KS_eigenvector, occ_numbers, max_occ_number, density_matrix)

   !-------shanghui begin debug_mode------
    if (module_is_debugged("DFPT")) then
      write(info_str,'(A)') '************shanghui begain zero_order_DM****************'
      call localorb_info(info_str, use_unit,'(A)',OL_norm)

      do i_basis=1,n_basis 
       write(info_str,'(60f15.9)') (density_matrix(i_basis,j_basis),j_basis=1,n_basis )
       call localorb_info(info_str, use_unit,'(A)', OL_norm)
      enddo
      write(info_str,'(A)') '************shanghui end zero_order_DM****************'
      call localorb_info(info_str, use_unit,'(A)', OL_norm )
      write(info_str,'(A)') ' '
      call localorb_info(info_str, use_unit,'(A)', OL_norm )
    end if
   !-------shanghui end debug_mode------

!----------- becase DFPT is working for the second derivative of E, so we use j_atom,j_coord
   do j_atom = 1, n_atoms 
   do j_coord =1 ,3

      write(info_str,'(A)') ''
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
      write(info_str,'(A)') ''
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
      write(info_str,'(A)') "==========================================================================="
      call localorb_info(info_str, use_unit,'(A)', OL_norm  )
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
     
   !--------------shanghui test first_order_S-------------------
       call get_timestamps(time_first_order_S, clock_time_first_order_S)

       call integrate_first_order_S_reduce_memory(partition_tab, l_shell_max, first_order_S,j_atom,j_coord)

       call get_times(time_first_order_S, clock_time_first_order_S, &
        &              tot_time_first_order_S, tot_clock_time_first_order_S)
     
   !-------shanghui begin debug_mode------
    if (module_is_debugged("DFPT")) then
      write(info_str,'(A)') '************shanghui begain first_order_S(atom1,X)****************'
      call localorb_info(info_str, use_unit,'(A)',OL_norm)

      do i_basis=1,n_basis 
       write(info_str,'(60f15.9)') (first_order_S(i_basis,j_basis),j_basis=1,n_basis )
       call localorb_info(info_str, use_unit,'(A)', OL_norm)
      enddo
      write(info_str,'(A)') '************shanghui end first_order_S****************'
      call localorb_info(info_str, use_unit,'(A)', OL_norm )
      write(info_str,'(A)') ' '
      call localorb_info(info_str, use_unit,'(A)', OL_norm )
    end if
   !-------shanghui end debug_mode------

   !--------------shanghui end test first_order_S---------------   


     


    first_order_H=0.0d0
    first_order_U=0.0d0 
    first_order_E=0.0d0

       call evaluate_first_order_DM_reduce_memory(first_order_S,  &
             KS_eigenvector, KS_eigenvalue, occ_numbers,max_occ_number,   &
             first_order_U,old_first_order_density_matrix & 
             ) 

! ------------------------ self-consistency loop -------------->>
  SCF_LOOP: do while ( (.not.converged) .and.  &
  &                    below_it_limit )
        number_of_loops = number_of_loops + 1

        write(info_str,'(A)') ''
        call localorb_info(info_str, use_unit,'(A)', OL_norm  )
        write(info_str,'(A)') "--------------------------------------------------------------"
        call localorb_info(info_str, use_unit,'(A)', OL_norm  )
        
        write(info_str,'(2X,A,1X,I4,5X,A,1X,I4)') "For j_atom =", j_atom, "j_coord =",j_coord
        call localorb_info(info_str, use_unit,'(A)', OL_norm )
        write(info_str,'(10X,A,1X,I4)') "Begin CP-self-consistency iteration #", number_of_loops
        call localorb_info(info_str, use_unit,'(A)', OL_norm )
        write(info_str,*) 'DFPT_mixing = ', DFPT_mixing
        call localorb_info(info_str, use_unit,'(A)', OL_norm )
        write(info_str,*) 'DFPT_sc_accuracy_dm = ', DFPT_sc_accuracy_dm
        call localorb_info(info_str, use_unit,'(A)', OL_norm )
        
        write(info_str,'(A)') "------------------------------------------------------------"
        call localorb_info(info_str, use_unit,'(A)', OL_norm )


       call get_timestamps ( time_cpscf_loop, clock_time_cpscf_loop )

!---------------(1) Begin  update first_order_rho-----------------------------------
       call get_timestamps(time_first_order_density, clock_time_first_order_density)

       call evaluate_first_order_DM_reduce_memory(first_order_S,  &
             KS_eigenvector, KS_eigenvalue, occ_numbers,max_occ_number,   &
             first_order_U,first_order_density_matrix)

         !shanghui: we need to perform density matrix mixing here
    
          change_of_first_order_DM =0.0d0         

         do i_basis = 1,n_basis
         do j_basis = 1,n_basis

         change_of_first_order_DM =                &
         max( change_of_first_order_DM,             &
         dabs(first_order_density_matrix(i_basis,j_basis)  &
         - old_first_order_density_matrix(i_basis,j_basis)) )
         
         first_order_density_matrix(i_basis,j_basis) =       &
         (1.0d0-DFPT_mixing)*old_first_order_density_matrix(i_basis,j_basis)+  &
         DFPT_mixing*first_order_density_matrix(i_basis,j_basis)
     
         old_first_order_density_matrix(i_basis,j_basis) =   &
         first_order_density_matrix(i_basis,j_basis)
       
         enddo
         enddo

         !shanghui debug H1 by adding DM1 as input 
         !first_order_density_matrix=0.078623504

        !-------shanghui begin debug_mode------
       if (module_is_debugged("DFPT")) then
        write(info_str,'(A)') "((((((((((((((((((((((((((((((((("
        call localorb_info(info_str, use_unit,'(A)', OL_norm)
        write(info_str,*) change_of_first_order_DM
        call localorb_info(info_str, use_unit,'(A)', OL_norm)
        write(info_str,'(A)') ")))))))))))))))))))))))))))))))))"
        call localorb_info(info_str, use_unit,'(A)', OL_norm)
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


        call integrate_first_order_rho_reduce_memory(partition_tab, l_shell_max,  &
             KS_eigenvector(:,:,1,1),density_matrix,first_order_density_matrix,max_occ_number, &
             first_order_rho,j_atom,j_coord)

        call get_times(time_first_order_density, clock_time_first_order_density, &
        &              tot_time_first_order_density, tot_clock_time_first_order_density)

!------------(1) end first-order-density update and mixing--------


!------------(2)begain to calculate first_order_H-----------------

        call get_timestamps(time_first_order_potential, clock_time_first_order_potential)

        call update_hartree_potential_p2_shanghui &
            ( hartree_partition_tab,first_order_rho(1:n_full_points),& 
            delta_v_hartree_part_at_zero, &
            delta_v_hartree_deriv_l0_at_zero, &
            multipole_moments, multipole_radius_sq, &
            l_hartree_max_far_distance, &
            outer_potential_radius )

        call sum_up_whole_potential_p2_shanghui &
                ( delta_v_hartree_part_at_zero, &
                delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
                partition_tab, first_order_rho(1:n_full_points), &
                first_order_potential(1:n_full_points),  & !<--------get first_order_DM_potential
                .false., multipole_radius_sq, &
                l_hartree_max_far_distance, &
                outer_potential_radius)

        call get_times(time_first_order_potential, clock_time_first_order_potential, &
        &              tot_time_first_order_potential, tot_clock_time_first_order_potential)


        call get_timestamps(time_first_order_H, clock_time_first_order_H)

        call  integrate_first_order_H_reduce_memory &
             (hartree_potential,first_order_potential, rho, rho_gradient,&
              first_order_rho, & 
              partition_tab, l_shell_max, &
              first_order_density_matrix, density_matrix, & 
              first_order_H, j_atom, j_coord &
             )

        call get_times(time_first_order_H, clock_time_first_order_H, &
        &              tot_time_first_order_H, tot_clock_time_first_order_H)

        !-------shanghui begin debug_mode------

      if (module_is_debugged("DFPT")) then
        if(myid.eq.0) then
        write(info_str,'(A)') '************shanghui begain first_order_H(atom1,X)****************'
         call localorb_info(info_str, use_unit,'(A)', OL_norm)
         do i_basis=1,n_basis
         write(info_str,'(60f15.9)') (first_order_H(i_basis,j_basis),j_basis=1,n_basis )
         call localorb_info(info_str, use_unit,'(A)', OL_norm)
         enddo
        write(info_str,'(A)') '************shanghui end first_order_H****************'
        call localorb_info(info_str, use_unit,'(A)' , OL_norm )
        write(info_str,'(A)') ''
        call localorb_info(info_str, use_unit,'(A)', OL_norm)
        endif
      endif
        !-------shanghui end debug_mode--------

       !shanghui debug H1 by adding DM1 as input 
       !stop
! ------------(2) end to calculate first_order_H-----------------




!--------------(3) begin to calculate first_order_U-----------------
        call get_timestamps(time_Sternheimer, clock_time_Sternheimer)

        call evaluate_first_order_U_reduce_memory(first_order_H, first_order_S,  &
             KS_eigenvector, KS_eigenvalue, occ_numbers,  max_occ_number, &
             first_order_U,first_order_E)

        call get_times(time_Sternheimer, clock_time_Sternheimer, &
        &              tot_time_Sternheimer, tot_clock_time_Sternheimer)

!--------------(3) end solve first_order_U problem----------------


! --------- Check convergence ----->>

        converged = (change_of_first_order_DM.lt.DFPT_sc_accuracy_dm).and.(number_of_loops.gt.1) 

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
        call output_times(deffmt, "Solution of Sternheimer eqns.", &
        &                 time_Sternheimer, clock_time_Sternheimer, OL_norm)
        write(info_str,'(A)') &
        "------------------------------------------------------------"
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

! << ---- end printing out data---------



  end do SCF_LOOP
! << ------ end self consistent cycle--------

      !total_number_of_loops = total_number_of_loops + number_of_loops


   !---------begin Hessian calculation-------------
   call get_timestamps(time_Hessian, clock_time_Hessian)

   !------------this is the old loop version, which is slower-----------------
   !call  evaluate_first_zero_order_DM_EDM_reduce_memory( first_order_S,first_order_H,  &
   !      KS_eigenvector, KS_eigenvalue, occ_numbers, max_occ_number,  &
   !      first_order_U,density_matrix,first_order_density_matrix, &
   !      energy_density_matrix,first_order_energy_density_matrix)

   !------------this is the new dgemm version, which is faster-----------------
    call evaluate_zero_order_DM_reduce_memory( & 
         KS_eigenvector, occ_numbers, max_occ_number, density_matrix)
    call evaluate_zero_order_EDM_reduce_memory( & 
         KS_eigenvector, KS_eigenvalue, occ_numbers, max_occ_number, energy_density_matrix)
    call evaluate_first_order_DM_reduce_memory( & 
         first_order_S, KS_eigenvector, KS_eigenvalue, occ_numbers, max_occ_number, & 
         first_order_U, first_order_density_matrix)
    call evaluate_first_order_EDM_reduce_memory( & 
         first_order_S, first_order_H, KS_eigenvector, KS_eigenvalue, occ_numbers, max_occ_number, & 
         first_order_U, first_order_energy_density_matrix)
 

      if (module_is_debugged("DFPT")) then
        if(myid.eq.0) then
        write(info_str,'(A)') '************shanghui begain first_order_DM(atom1,X)****************'
         call localorb_info(info_str, use_unit,'(A)', OL_norm)
         do i_basis=1,n_basis
         write(info_str,'(60f15.9)') (first_order_density_matrix(i_basis,j_basis),j_basis=1,n_basis )
         call localorb_info(info_str, use_unit,'(A)', OL_norm)
         enddo
        write(info_str,'(A)') '************shanghui end first_order_DM****************'
        call localorb_info(info_str, use_unit,'(A)' , OL_norm )
        write(info_str,'(A)') ''
        call localorb_info(info_str, use_unit,'(A)', OL_norm)
        endif
      endif

   call integrate_hessian_reduce_memory &
     ( hartree_potential, first_order_potential, & 
       rho, rho_gradient,  &
       partition_tab,  l_shell_max,     &
       KS_eigenvector(:,:,1,1),       &
       density_matrix, first_order_density_matrix,  & 
       energy_density_matrix, first_order_energy_density_matrix,  &
       max_occ_number, &
       hessian, gradient_dipole, &
       j_atom,j_coord  & 
     )

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



   !---------end Hessian calculation-------------
   enddo !j_coord
   enddo !j_atom


!--------------------shanghui make asr 2013.10.21----------------------
    if(use_ASR.and.myid.eq.0) then

      if (module_is_debugged("DFPT")) then
         write(use_unit,*) ''
         write(use_unit,*) 'Total_hessian in cpscf.f90:'

         do i_atom = 1 , n_atoms 
         do i_coord = 1 , 3

         write(use_unit,'(100E20.12)')  & 
             ((  hessian(3*i_atom+i_coord-3,3*j_atom+j_coord-3) & 
               , j_coord=1,3) ,j_atom=1,n_atoms)
        enddo
        enddo 
      endif ! debug_mode


      write(use_unit,*) ''
      write(use_unit,*) 'Make Acoustic sum rule in cpscf.f90:'

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

   endif
!--------------------shanghui end make asr 2013.10.21----------------------


  if(myid.eq.0) then
      write(use_unit,*) ''
      write(use_unit,*) 'grdient_dipole in cpscf_solver.f90:'
      do i_atom = 1 , n_atoms
      do i_coord = 1 , 3
         write(use_unit,'(3E20.12)')  gradient_dipole(3*i_atom+i_coord-3,1:3)
      enddo
      enddo  
  endif


   !---------begin vib calculation,----shanghui begin parallel------
   if(myid.eq.0) then

        write(info_str,'(A)') ''
        call localorb_info(info_str, use_unit,'(A)', OL_norm)
        write(info_str,'(A)') 'Get frequencies in cpscf.f90:'
        call localorb_info(info_str, use_unit,'(A)', OL_norm)
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
    
      do i_coord = 1, 3*n_atoms
       do i_coord_2 = 1, i_coord - 1 
          buf = (hessian(i_coord_2,i_coord)+hessian(i_coord,i_coord_2))/2d0
          hessian(i_coord_2,i_coord) = buf
          hessian(i_coord,i_coord_2) = buf
       end do
      end do

      write(info_str,'(A)') 'Solving eigenvalue system for Hessian Matrix'
      call localorb_info(info_str, use_unit,'(A)', OL_norm)
      call DSYEV('V','U',3*n_atoms,hessian,3*n_atoms,eigenvalues,workspace,lwork,info)
      write(info_str,'(A)') 'Done ... '
      call localorb_info(info_str, use_unit,'(A)', OL_norm)
      write(info_str,'(A)') ''
      call localorb_info(info_str, use_unit,'(A)', OL_norm)
      
      ! calculate the eigenvectors in cartesian coordinates ?
      do i_coord = 1, 3*n_atoms
         hessian(:,i_coord) = hessian(:,i_coord)*mass_vector(:)
      end do


    !---------------------begin for IR intensity--------------------------------------   
    gradient_dipole(:,:) = gradient_dipole(:,:) * gradient_dipole_factor
    ! transform dipole derivative to internal coordinates via directional derivative
    ! d/dQ = d/dR * Q_normalized, where Q_normalized is nothing but displacement eigenvector
    ! hessian(:,:) not normalized anymore since mass was divided out
     do i_coord = 1, 3
       ! loop over modes
       do i_coord_2 = 1, 3*n_atoms
          gradient_dipole_internal(i_coord, i_coord_2) = ddot(3*n_atoms, gradient_dipole(:, i_coord), 1, hessian(:, i_coord_2), 1)
       end do
     end do
 
    ! get infrared intensities
     do i_coord = 1, 3*n_atoms, 1
       infrared_intensity(i_coord) = ddot(3, gradient_dipole_internal(:, i_coord), 1, gradient_dipole_internal(:, i_coord), 1)
     end do

    ! scale infrared intensities
     infrared_intensity(:) = infrared_intensity(:)  !* ir_factor
    !---------------------end for IR intensity--------------------------------------   

 
    ! Renormalize eigenvectors for output - norm has units of 1/sqrt(mass) though 
     do i_coord = 1, 3*n_atoms
      reduced_mass(i_coord) = ddot(3*n_atoms, hessian(:, i_coord), 1, hessian(:, i_coord), 1)
      norm = sqrt(reduced_mass(i_coord))
      hessian(:,i_coord) = hessian(:,i_coord)/norm
     end do

     write(info_str,'(A)') 'DFPT-Results: '
     call localorb_info(info_str, use_unit,'(A)', OL_norm)
     write(info_str,'(A)') ''
     call localorb_info(info_str, use_unit,'(A)', OL_norm)
     write(info_str,'(A)') 'List of all frequencies found:'
     call localorb_info(info_str, use_unit,'(A)', OL_norm)
     write(info_str,'(A13,A25,A27)') 'Mode number','Frequency [cm^(-1)]','IR-intensity [D^2/Ang^2]'
     call localorb_info(info_str, use_unit,'(A)', OL_norm)

     do i_coord = 1, 3*n_atoms, 1
       if (eigenvalues(i_coord).gt.0d0) then
          frequencies(i_coord) = sqrt(eigenvalues(i_coord))
       else if (eigenvalues(i_coord).lt.0d0) then
          frequencies(i_coord) = -sqrt(-eigenvalues(i_coord))
       else 
          frequencies(i_coord) = 0d0
       end if

       write(info_str, '(I13,F25.8,F27.8)') i_coord, frequencies(i_coord)/(200*pi*const_c),infrared_intensity(i_coord)
       call localorb_info(info_str, use_unit,'(A)', OL_norm)
     enddo

   endif 
   !-------shanghui end parallel--------


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


       deallocate(first_order_rho)
       deallocate(first_order_potential)

       !----shanghui begin for vib-cal-------------
       deallocate( workspace )
       deallocate( eigenvalues )
       deallocate( mass_vector )
       deallocate( reduced_mass )
       deallocate( frequencies )
       deallocate( hessian )
       deallocate( gradient_dipole)
       deallocate( gradient_dipole_internal)
       deallocate(infrared_intensity)

       !----shanghui end for vib-cal-------------



    end subroutine cpscf_solver_reduce_memory
!******
