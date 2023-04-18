!****s* FHI-aims/advance_KS_solution
!  NAME
!    advance_KS_solution
!  SYNOPSIS
subroutine advance_KS_solution(overlap_matrix, hamiltonian, n_electrons, &
              KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, &
              occ_numbers, chemical_potential, chemical_potential_spin)
!  PURPOSE
!    This subroutine is (should be) an advance version of the old
!    subroutine get_KS_orbitals_p0. It uses the newly-developed
!    ELSI (ELectronic Structure Infrastructure) package to treat
!    the Kohn-Sham eigenvalue problem.
!  AUTHOR
!    FHI-aims team
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is
!    subject to the terms and conditions of the respective license
!    agreement.
!  NOTE
!    Visit http://elsi-interchange.org for more information on ELSI.
!    For any issue with this subroutine, please let the developers
!    know by sending an email to elsi-team@duke.edu.
!
!  USES
   use dimensions
   use geometry
   use basis
   use runtime_choices
   use constraint
   use localorb_io
   use mixing_constraint
   use bfgs
   use constants
   use synchronize_mpi
   use scalapack_wrapper
   use separate_core_states
   use plus_u
   use force_occupation
   use elsi_wrapper, only: eh_scf,aims_elsi_set_mu_spin_degen,&
                           aims_elsi_compute_mu_and_occ
   use pbc_lists, only: k_weights

   implicit none

!  ARGUMENTS
   real*8 :: hamiltonian(n_hamiltonian_matrix_size,n_spin)
   real*8 :: overlap_matrix(n_hamiltonian_matrix_size)

   real*8, dimension(n_states,n_spin,n_k_points) :: occ_numbers
   real*8, dimension(n_states,n_spin,n_k_points) :: KS_eigenvalue

   real*8,     dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector
   complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex

   real*8 :: n_electrons
   real*8 :: chemical_potential

!  INPUTS
!   o overlap_matrix -- overlap matrix.
!   o hamiltonian -- Hamiltonian matrix.
!   o n_electrons -- number of electrons in system.
!
!  OUTPUT
!   o KS_eigenvalue -- Kohn-Sham eigenvalues.
!   o KS_eigenvector -- Kohn-Sham eigenvectors if real number eigenvectors are in use.
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors if complex number eigenvectors are in use.
!   o occ_number -- occupation weights of the eigenstates.
!   o chemical_potential -- value of chemical potential.
!
!  SOURCE

   real*8, dimension(n_states,n_region,n_spin) :: constraint_proj
   real*8, dimension(n_region,n_spin) :: constraint_fermi
   real*8, dimension(n_region,n_spin) :: electrons_in_region
   real*8, dimension(n_spin) :: constraint_fermitot

   integer, dimension(n_spin) :: n_eval
   integer, dimension(n_spin) :: n_homo
   integer :: info, n_spin_tmp

   real*8,dimension(n_states,2) :: temp_n_electrons_forced_state

   external :: f_external_constraint

   logical :: constraint_converged
   logical :: t_out = .true.

   real*8 :: avg_zero
   real*8 :: avg_potential

   character*120 :: info_str

   ! counters
   integer :: i_spin,j_spin,i_k_point
   integer :: i_basis,j_basis,k_basis,l_basis
   integer :: i_region
   integer :: i_state
   integer :: i_hamiltonian
   integer :: number_constraint_iter

   integer :: n_arguments
   real*8, dimension(:), allocatable :: arguments
   real*8 :: eps
   real*8 :: fmin

   real*8, dimension(:), allocatable :: gradient
   real*8, dimension(:), allocatable :: hessian
   real*8, dimension(:), allocatable :: work
   real*8 :: dfn
   real*8, dimension(:), allocatable :: xm
   real*8 :: hh

   integer :: mode
   integer :: maxfn
   integer :: iexit,i_k

   real*8, dimension(n_spin) :: n_electrons_in_spin
   real*8, dimension(n_spin) :: chemical_potential_spin
   real*8 :: diff_in_chem_pot

   real*8 :: spin_shift_fsm(n_spin)

   if(.not. (force_occupation_projector .and. start_force_occ)) then
      if(use_elsi_ev) then
         write(info_str,'(A)') ""
         call localorb_info(info_str,use_unit,'(A)',OL_norm)

         if(elsi_solver == 1) then
            if(use_scalapack) then
               write(info_str,'(2X,2A)') &
                  "Updating Kohn-Sham eigenvalues and eigenvectors using ELSI",&
                  " and the ELPA eigensolver."
            else
               write(info_str,'(2X,2A)') &
                  "Updating Kohn-Sham eigenvalues and eigenvectors using ELSI",&
                  " and the (modified) LAPACK eigensolver."
            endif
         elseif(elsi_solver == 5) then
            write(info_str,'(2X,2A)') &
               "Updating Kohn-Sham eigenvalues and eigenvectors using ELSI",&
               " and the SLEPc-SIPs eigensolver."
         endif

         call localorb_info(info_str,use_unit,'(A)',OL_norm)
      elseif(use_elsi_dm) then
         write(info_str,'(A)') ""
         call localorb_info(info_str,use_unit,'(A)',OL_norm)

         if(elsi_solver == 1) then
            write(info_str,'(2X,A)') &
               "Updating density matrix using ELSI and the ELPA eigensolver."
         elseif(elsi_solver == 2) then
            write(info_str,'(2X,2A)') &
               "Updating density matrix using ELSI and the libOMM density",&
               " matrix solver."
         elseif(elsi_solver == 3) then
            write(info_str,'(2X,2A)') &
               "Updating density matrix using ELSI and the PEXSI density",&
               " matrix solver."
         elseif(elsi_solver == 5) then
            write(info_str,'(2X,2A)') &
               "Updating density matrix using ELSI and the SLEPc-SIPs",&
               " eigensolver."
         elseif(elsi_solver == 6) then
            write(info_str,'(2X,2A)') &
               "Updating density matrix using ELSI and the NTPoly density",&
               " matrix solver."
         endif

         call localorb_info(info_str,use_unit,'(A)',OL_norm)
      else
         write(info_str,'(A)') ""
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
         write(info_str,'(2X,2A)') &
               "Updating Kohn-Sham eigenvalues, eigenvectors, occupation",&
               " numbers, and Fermi level."
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
      endif

      if(.not.use_cg) then
         occ_numbers = 0.0d0
      endif

      ! if locally constrained DFT requested, prepare for it
      if(use_constraint) then
         ! add constraint potential(s) from previous iterations to Hamiltonian
         ! In case of a constraint calculation we need to store the original hamiltonian
         ! because we need it for any inner iteration
         ! We also allocate a full work version of the ovlp matrix here.

         call allocate_hamiltonian_work(overlap_matrix)

         ! copy original hamiltonian to hamiltonian_work, to save it from
         ! destruction by LAPACK.
         do i_spin=1,n_spin,1
            do i_hamiltonian=1,n_basis*(n_basis+1)/2
               hamiltonian_work(i_hamiltonian,i_spin) = hamiltonian(i_hamiltonian,i_spin)
            enddo
         enddo

         ! If a constraint was requested, we may already know constraint potentials for different regions.
         ! These must be added to hamiltonian_work.

         if(n_active_regions.gt.1) then
            call add_constraint_potentials(overlap_matrix,hamiltonian)
         endif

         t_out = .false.

      endif

      ! Solve KS eigenvalue problem
      call solve_KS_eigen(overlap_matrix, hamiltonian, KS_eigenvalue, &
           KS_eigenvector, KS_eigenvector_complex)

      ! Obtain occupation numbers if using eigensolver
      if(.not. use_elsi_dm) then
         if(.not. fixed_spin_moment) then
            if(mu_method == 0) then ! zeroin method
               call get_occupation_numbers_p0(KS_eigenvalue, n_electrons, t_out, &
                                              occ_numbers, chemical_potential)
               chemical_potential_spin(1:n_spin) = chemical_potential
            elseif(mu_method == 1) then ! bisection method in ELSI
               write(info_str,'(2X,A)') ''
               call localorb_info(info_str,use_unit,'(A)',OL_norm)
               write(info_str,'(2X,A)') &
                  "Obtaining occupation numbers and chemical potential using ELSI."
               call localorb_info(info_str,use_unit,'(A)',OL_norm)

               if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
                 !(Rundong) For fully-relativistic cases, due to the mismatch in the definition of
                 ! occupation, I have to set n_spin=1 here, in order to use the existing code in aims.
                 call aims_elsi_compute_mu_and_occ(eh_scf, n_electrons, n_states, &
                       1, n_k_points, k_weights, KS_eigenvalue(:,1,:), occ_numbers(:,1,:), &
                       chemical_potential)
               else
                 call aims_elsi_compute_mu_and_occ(eh_scf, n_electrons, n_states, &
                       n_spin, n_k_points, k_weights, KS_eigenvalue, occ_numbers, &
                       chemical_potential)
               endif

               if(t_out) then
                  write(info_str,'(2X,A,F15.8,A)') "| Chemical potential (Fermi level):", &
                     chemical_potential*hartree, " eV"
                  call localorb_info(info_str,use_unit,'(A)',OL_norm)
               endif

               chemical_potential_spin(1:n_spin) = chemical_potential
            endif
         else
            do i_spin = 1, n_spin, 1
               ! determine occupation numbers for each channel directly.
               ! Note that for this to work always, MUST have the shapes
               ! of KS_eigenvalues, occ_numbers arrays declared explicitly
               ! in the subroutine header of get_KS_orbitals_p0 above,
               ! NOT by simply USEing the array shapes from module physics.f90
               ! directly. Reason is that the dimension n_states can change as
               ! the code progresses, in rare but possible cases.
               if(mu_method == 0) then ! zeroin method
                  call get_occupation_numbers_single_channel(&
                          KS_eigenvalue(1:n_states,i_spin,1:n_k_points),&
                          fixed_spin_moment_electrons(i_spin),t_out,&
                          occ_numbers(1:n_states,i_spin,1:n_k_points),&
                          chemical_potential_spin(i_spin),i_spin)
               elseif(mu_method == 1) then ! bisection method in ELSI
                  call aims_elsi_set_mu_spin_degen(eh_scf, 1.0d0)

                  call aims_elsi_compute_mu_and_occ(eh_scf, &
                          fixed_spin_moment_electrons(i_spin), n_states, 1, &
                          n_k_points, k_weights, &
                          KS_eigenvalue(1:n_states,i_spin,1:n_k_points), &
                          occ_numbers(1:n_states,i_spin,1:n_k_points), &
                          chemical_potential_spin(i_spin))
               endif
            enddo

            ! Calculate overall chemical potential
            diff_in_chem_pot = chemical_potential_spin(1) - chemical_potential_spin(2)

            write(info_str, '(2X,A)') &
                  'Fixed spin moments - chemical potentials:'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str, '(2X,A,F15.10,A)') &
                  '| Chemical potential, spin up: ', &
                  chemical_potential_spin(1)*hartree, &
                  " eV."
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str, '(2X,A,F15.10,A)') &
                  '| Chemical potential, spin dn: ', &
                  chemical_potential_spin(2)*hartree, &
                  " eV."
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str, '(2X,A,A,F15.10,A)') &
                  '| Fixed spin moment: Difference in chemical potentials between ', &
                  'spin channels:', diff_in_chem_pot*hartree, &
                  " eV."
            call localorb_info(info_str,use_unit,'(A)')
            call localorb_info('',use_unit,'(A)')

            spin_shift_fsm(1) = fixed_spin_moment_electrons(2)/n_electrons * &
                                diff_in_chem_pot

            chemical_potential = chemical_potential_spin(1) - spin_shift_fsm(1)

         endif ! fixed-spin moment or not
      endif ! use_elsi_ev

      ! At this point, if an electron number constraint was requested, attempt to enforce
      ! it using auxiliary potentials

      if(use_constraint) then
         ! count number of electrons in each region, each spin
         do i_spin = 1, n_spin, 1
            call get_electrons_per_region_v2(i_spin, KS_eigenvector, overlap_matrix, &
                                             occ_numbers, constraint_proj, electrons_in_region)
         enddo

         if(constraint_it_lim.gt.1) then
            ! if iterative determination of constraint_potential requested, check whether
            ! the correct electron numbers have already been reached

            constraint_converged = .true.

            do i_spin = 1,n_spin,1
               do i_region = 1,n_active_regions,1
                  constraint_converged = constraint_converged .and. &
                     (dabs(electrons_in_region(i_region,i_spin)- &
                     constraint_electrons(i_region,i_spin)).lt. &
                     constraint_precision)
               enddo
            enddo

            if(constraint_converged) then
               write(info_str,'(2X,A,A)') &
                     "Locally constrained DFT: ", &
                     "Electron counting constraint fulfilled on entry."
               call localorb_info(info_str,use_unit,'(A)')
            else
               write(info_str,'(2X,A,A)') &
                     "Locally constrained DFT: ", &
                     "Constraint on electron numbers not yet fulfilled."
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
                     "Entering inner self-consistency loop."
               call localorb_info(info_str,use_unit,'(A)')

               ! restore hamiltonian
               hamiltonian = hamiltonian_work
            endif
         endif

         ! Begin inner self-consistency loop to determine
         ! constraint potentials
         if((mixer_constraint.eq.1)) then
            call allocate_pulay_constraint( )
         endif

         number_constraint_iter = 0

         if(mixer_constraint.le.1) then
            ! If constraint already converged, still do a final adjustment of
            ! the eigenvalue levels etc, just for consistency's sake!!
            !
            ! This is a patch for a problem where significant shifts of the
            ! eigenvalues of spin-up vs spin-down against each other can arise
            ! even though the overall spin-up <-> spin-down occupancies
            ! do not change at all.
            if(constraint_converged) then
               n_electrons_in_spin = 0.d0

               do i_spin = 1,n_spin,1
                  do i_region = 1,n_active_regions,1
                     n_electrons_in_spin(i_spin) = n_electrons_in_spin(i_spin) + &
                                                   constraint_electrons(i_region,i_spin)
                  enddo

                  call get_occupation_numbers_v2(KS_eigenvalue(:,i_spin,1), &
                                                 n_electrons_in_spin(i_spin), t_out, &
                                                 occ_numbers(:,i_spin,1), &
                                                 chemical_potential_spin(i_spin))

                  ! The number of electrons in region 1 is now by definition correct.
                  do i_region = 1, n_active_regions, 1
                     electrons_in_region(i_region,i_spin) = constraint_electrons(i_region,i_spin)
                  enddo
               enddo

               if(n_spin.eq.2) then

                  diff_in_chem_pot = chemical_potential_spin(1) - chemical_potential_spin(2)

                  if(constraint_debug) then
                     call localorb_info('',use_unit,'(A)')
                     write(info_str, '(2X,A,A,F15.10,A)') &
                           'Difference in chemical potential between ', &
                           'spin channels:', diff_in_chem_pot*hartree, &
                           " eV."
                     call localorb_info(info_str,use_unit,'(A)')
                  endif

                  spin_shift(1) = n_electrons_in_spin(2)/n_electrons * diff_in_chem_pot
                  spin_shift(2) = - n_electrons_in_spin(1)/n_electrons * diff_in_chem_pot

                  do i_spin = 1,n_spin,1

                     ! adjust constraint potentials to reflect the spin shift
                     do i_region = 1,n_active_regions,1

                        constraint_potential(i_region,i_spin) =  &
                           constraint_potential(i_region,i_spin) - spin_shift(i_spin)
                     enddo

                     ! adjust the final eigenvalues to reflect the spin shift
                     do i_state = 1,n_states,1
                        KS_eigenvalue(i_state,i_spin,1) = KS_eigenvalue(i_state,i_spin,1) - &
                                                          spin_shift(i_spin)
                     enddo

                     chemical_potential_spin(i_spin) = chemical_potential_spin(i_spin) - &
                                                       spin_shift(i_spin)
                  enddo

                  ! Both spin channel chemical potentials must now be equal
                  ! simply by construction!
                  chemical_potential = chemical_potential_spin(1)

                  ! end adjustment of different spin channels
               endif
               ! end the patch if the constraint was already converged.
            endif

            do while ((.not.constraint_converged) .and. &
                     (number_constraint_iter.lt.constraint_it_lim))

               number_constraint_iter = number_constraint_iter+1

               if(constraint_debug) then
                  write(info_str,*)
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') &
                        "------------------------------------------------------"
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A,I5)') &
                        "Locally constrained DFT: Inner SCF iteration no. ", &
                        number_constraint_iter
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)')  &
                        "------------------------------------------------------"
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,*)
                  call localorb_info(info_str,use_unit,'(A)')
               endif

               ! store present constraint potential for later mixing
               constraint_pot_old = constraint_potential

               ! Determine region-specific Fermi levels which would be needed to
               ! satisfy the constraint in each region, based on the present
               ! KS eigenstates alone
               do i_region = 1,n_active_regions,1
                  do i_spin = 1,n_spin,1
                     call get_constraint_fermi(KS_eigenvalue(:,i_spin,1), &
                             constraint_electrons(i_region,i_spin), &
                             constraint_proj(:,i_region,i_spin), &
                             constraint_fermi(i_region,i_spin))
                  enddo
               enddo

               ! write results
               if(constraint_debug) then
                  write(info_str,*)
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)')  &
                        "Required Fermi levels in different subsystems:"
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A,5X,A,I2,5X,A,I5)')  &
                        "| Region  ", ("Spin ", i_spin, i_spin=1,n_spin,1)
                  call localorb_info(info_str,use_unit,'(A)')
                  do i_region = 1,n_active_regions,1
                     write(info_str,'(2X,A,I6,2X,F10.5,2X,F10.5)') &
                           "| ", &
                           constraint_region_label(i_region), &
                           (constraint_fermi(i_region,i_spin), i_spin=1,n_spin,1)
                     call localorb_info(info_str,use_unit,'(A)')
                  enddo
                  write(info_str,*)
                  call localorb_info(info_str,use_unit,'(A)')
               endif

               ! compute average zero level to which all potentials will be referenced
               ! variant 1: avg_zero = constraint_fermi(1,1)
               ! variant 2: avg_zero = chemical_potential
               ! variant 3:
               avg_zero = 0.d0
               do i_spin = 1,n_spin,1
                  do i_region = 1,n_active_regions,1
                     avg_zero = avg_zero + &
                        constraint_electrons(i_region,i_spin) * &
                        constraint_fermi(i_region,i_spin)
                  enddo
               enddo
               avg_zero = avg_zero / n_electrons

               ! To each region, we apply a constraint potential that would shift
               ! the local Fermi level to the current overall Fermi level.
               ! We arbitrarily set the constraint potential of the first region
               ! to zero ...
               do i_region=1,n_active_regions,1
                  do i_spin=1,n_spin,1
                     delta_constraint_potential(i_region,i_spin)= &
                        avg_zero - constraint_fermi(i_region,i_spin)
                  enddo
               enddo

               ! write results
               if(constraint_debug) then
                  write(info_str,*)
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A,A)')  &
                        "Suggested local potential change", &
                        " in different subsystems:"
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A,5X,A,I2,5X,A,I5)')  &
                        "| Region  ", ("Spin ", i_spin, i_spin=1,n_spin,1)
                  call localorb_info(info_str,use_unit,'(A)')
                  do i_region = 1,n_active_regions,1
                     write(info_str,'(2X,A,I6,2X,F10.5,2X,F10.5)') &
                           "| ", &
                           constraint_region_label(i_region), &
                           (delta_constraint_potential(i_region,i_spin), &
                           i_spin=1,n_spin,1)
                     call localorb_info(info_str,use_unit,'(A)')
                  enddo
                  write(info_str,*)
                  call localorb_info(info_str,use_unit,'(A)')
               endif

               ! Damp additionally needed constraint potential by simple linear mixing factor

               if((mixer_constraint.eq.0).or.(number_constraint_iter.le. &
                  ini_linear_mixing_constraint)) then

                  do i_region = 1,n_active_regions,1
                     do i_spin = 1,n_spin,1
                        constraint_potential(i_region,i_spin) = &
                           constraint_pot_old(i_region,i_spin) + &
                           constraint_mix(1) * &
                           delta_constraint_potential(i_region,i_spin)
                     enddo
                  enddo
               endif

               if(mixer_constraint.eq.1) then
                  if(number_constraint_iter.le.ini_linear_mixing_constraint) then
                     do i_spin = 1,n_spin,1
                        do i_region = 1,n_active_regions,1
                           delta_potential_mixing(i_region,i_spin) = &
                              delta_constraint_potential(i_region,i_spin) &
                              - delta_constraint_potential(1,i_spin)
                        enddo
                     enddo

                     call prepare_pulay_mixing_constraint(delta_potential_mixing)

                  else
                     do i_spin = 1,n_spin,1
                        do i_region = 1,n_active_regions,1
                           delta_potential_mixing(i_region,i_spin) = &
                              delta_constraint_potential(i_region,i_spin) &
                              - delta_constraint_potential(1,1)
                        enddo
                     enddo

                     call execute_pulay_mixing_constraint(number_constraint_iter, &
                                                          constraint_potential, &
                                                          delta_potential_mixing)
                  endif
               endif

               ! write results
               if(constraint_debug) then
                  write(info_str,*)
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A,A)') &
                        "Local constraint potentials in different subsystems after mixing:"
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A,5X,A,I2,5X,A,I5)') "| Region  ", ("Spin ", i_spin, i_spin=1,n_spin,1)
                  call localorb_info(info_str,use_unit,'(A)')
                  do i_region = 1,n_active_regions,1
                     write(info_str,'(2X,A,I6,2X,F10.5,2X,F10.5)')  "| ", &
                           constraint_region_label(i_region), &
                           (constraint_potential(i_region,i_spin),i_spin=1,n_spin,1)
                     call localorb_info(info_str,use_unit,'(A)')
                  enddo
                  write(info_str,*)
                  call localorb_info(info_str,use_unit,'(A)')
               endif

               ! Restore original hamiltonian ...
               do i_spin=1,n_spin,1
                  do i_hamiltonian=1, n_basis*(n_basis+1)/2
                     hamiltonian_work(i_hamiltonian,i_spin) = hamiltonian(i_hamiltonian,i_spin)
                  enddo
               enddo

               ! ... and add the current constraint potentials.
               call add_constraint_potentials(overlap_matrix,hamiltonian_work)

               ! dummy setting for output purposes
               ! constrained DFT is currently implemented only for the non-periodic case, anyway.
               i_k_point = 1

               ! solve for new KS eigenstates and -energies ...
               do i_spin = 1,n_spin,1

                  call improve_real_eigenfunctions(overlap_matrix, hamiltonian_work(:,i_spin), &
                                                   t_out, KS_eigenvalue(:,i_spin,1), &
                                                   KS_eigenvector(:,:,i_spin,1), i_k_point)
               enddo

               ! ... and obtain the updated global Fermi level and
               ! and global occupation numbers
               call get_occupation_numbers(KS_eigenvalue, n_electrons, t_out, &
                                           occ_numbers, chemical_potential)

               ! update constraint projectors and count number of electrons in each region, each spin
               do i_spin = 1,n_spin,1
                  call get_electrons_per_region_v2(i_spin, KS_eigenvector, overlap_matrix, &
                                                   occ_numbers, constraint_proj, electrons_in_region)
               enddo

               ! Finally, check for convergence:
               constraint_converged = .true.

               do i_spin = 1,n_spin,1
                  do i_region = 1,n_active_regions,1
                     constraint_converged = constraint_converged .and. &
                        (dabs(electrons_in_region(i_region,i_spin)- &
                        constraint_electrons(i_region,i_spin)).lt. &
                        constraint_precision)
                  enddo
               enddo

               if(constraint_debug) then
                  if(constraint_converged) then
                     write(info_str,'(2X,A,A)') &
                           "Electron counting constraint fulfilled - ", &
                           "inner scf cycle converged."
                     call localorb_info(info_str,use_unit,'(A)')
                  else
                     write(info_str,'(2X,A,A,A)') &
                           "Constraint on electron numbers not yet fulfilled."
                     call localorb_info(info_str,use_unit,'(A)')
                  endif
               endif
               ! end inner SCF loop.
            enddo
            ! endif (mixer_constraint.le.1)
         endif

         if(mixer_constraint.eq.2) then
            ! call BFGS implementation
            ! first, allocations for BFGS
            n_arguments = n_active_regions - 1

            if((n_arguments.gt.0).and.(.not.constraint_converged)) then
               ! This is where we actually run through the optimizer again ...

               if(.not.allocated(arguments)) then
                  allocate(arguments(n_arguments),stat=info)
                  call check_allocation(info,'arguments')
               endif
               if(.not.allocated(xm)) then
                  allocate( xm(3*n_arguments),stat=info)
                  call check_allocation(info,'xm')
               endif

               ! parameters for BFGS
               ! (consult the module bfgs.f)

               ! magnitude of the solution, xm > 0
               xm = 0.1d0

               ! step lenght for calculating the gradient
               hh = 1.0e-6

               ! tolerance of convergence with respect to the argumetns, x
               eps = constraint_precision

               ! mode of utilisation, mode=1 corresponds to no Hessian given
               mode = 1

               ! maximum number of function evaluations premitted
               maxfn = constraint_it_lim

               ! loop over spins, so that the spin channels don't mix here
               do i_spin = 1,n_spin,1

                  ! initialize number of required electrons per spin channel
                  n_electrons_in_spin(i_spin) = 0.0d0
                  do i_region = 1,n_active_regions,1
                     n_electrons_in_spin(i_spin) = n_electrons_in_spin(i_spin) + &
                                                   constraint_electrons(i_region,i_spin)
                  enddo

                  if(constraint_debug) then
                     write(info_str,'(2X,A,A,I5,A)') &
                           '------------------',  &
                           ' Iteration for spin channel ', i_spin, &
                           ' ------------------'
                     call localorb_info(info_str,use_unit,'(A)')
                     write(info_str,'(2X,A)') &
                           'Region               Constraint potential'
                     call localorb_info(info_str,use_unit,'(A)')
                     do i_region = 1,n_active_regions,1
                        write(info_str,'(2X,I5,20X,F15.10)') &
                              i_region, &
                              constraint_potential(i_region,i_spin)*hartree
                        call localorb_info(info_str,use_unit,'(A)')
                     enddo
                  endif

                  ! preload the chemical potential
                  chemical_potential_spin(i_spin) = chemical_potential

                  ! prepare the arguments according to the solution last iteration
                  do i_region = 1,n_arguments,1
                     arguments(i_region) = constraint_potential(i_region,i_spin) - &
                                           constraint_potential(n_active_regions,i_spin)
                  enddo

                  ! initialize fmin
                  fmin = 0.d0

                  ! dfn is the expected change in the target functional
                  dfn = 0.0d0
                  do i_region = 1,n_active_regions,1
                     dfn = dfn + (electrons_in_region(i_region,i_spin) - &
                           constraint_electrons(i_region,i_spin))**2
                  enddo

                  ! call bfgs_constraint_v2 from the bfgs.f module
                  call bfgs_constraint_v2(f_external_constraint, n_arguments, &
                                          arguments, maxfn, eps, dfn, xm, hh, &
                                          number_constraint_iter, hamiltonian, &
                                          overlap_matrix, KS_eigenvalue, KS_eigenvector, &
                                          occ_numbers, n_electrons_in_spin(i_spin), &
                                          chemical_potential_spin(i_spin), &
                                          constraint_proj, electrons_in_region, i_spin)
                  n_eval(i_spin) = number_constraint_iter

                  ! end loop over spin channels
               enddo

            elseif ((n_arguments.eq.0).or.constraint_converged) then
               ! Either: Only one region, but with a spin constraint.
               ! Or: The constraint was already converged upon entry.
               ! In either case, simply redetermine needed Fermi level
               ! shifts for each spin channel in a single shot once, but
               ! do not affect the balance within each individual spin
               ! channel.

               do i_spin = 1,n_spin,1

                  n_electrons_in_spin(i_spin) = 0.d0
                  do i_region = 1,n_active_regions,1
                     n_electrons_in_spin(i_spin) = &
                         n_electrons_in_spin(i_spin) + &
                         constraint_electrons(i_region,i_spin)
                  enddo

                  call get_occupation_numbers_v2(KS_eigenvalue(:,i_spin,1), &
                                                 n_electrons_in_spin(i_spin), &
                                                 t_out, occ_numbers(:,i_spin,1), &
                                                 chemical_potential_spin(i_spin))

                  ! Do one final count of the electrons in each region,
                  ! just to catch any possible residual deviations within
                  ! the constraint convergence criterion.
                  !
                  ! THIS IS IMPORTANT TO ENSURE THAT THE ENERGY CORRECTION
                  ! DUE TO EIGENVALUE SHIFTS IS EXACTLY ZERO ON AVERAGE, AS IT SHOULD BE.
                  if(n_active_regions.gt.1) then
                     call get_electrons_per_region_v2(i_spin, KS_eigenvector, &
                                                      overlap_matrix, occ_numbers, &
                                                      constraint_proj, electrons_in_region)
                  else
                     electrons_in_region(1,i_spin) = constraint_electrons(1,i_spin)
                  endif
               enddo

               n_eval = 0

            else
               write(use_unit,*) "* No constraint regions defined?"
               stop
            endif

            if(n_spin.eq.2) then

               ! Balance the fermi levels of the spin channels

               diff_in_chem_pot = chemical_potential_spin(1) - chemical_potential_spin(2)

               if(constraint_debug) then
                  call localorb_info('',use_unit,'(A)')
                  write(info_str, '(2X,A,A,F15.10,A)')  &
                        'Difference in chemical potential between ', &
                        'spin channels:', diff_in_chem_pot*hartree, &
                        " eV."
                  call localorb_info(info_str,use_unit,'(A)')
               endif

               ! Shift all constraint potentials between the two spin channels
               ! to create one uniform Fermi level for the overall system

               ! Yes, the weighted average requires that 1 <-> 2 are exchanged
               ! in the equations.
               spin_shift(1) = n_electrons_in_spin(2)/n_electrons * diff_in_chem_pot
               spin_shift(2) = - n_electrons_in_spin(1)/n_electrons * diff_in_chem_pot

               do i_spin = 1,n_spin,1

                  ! adjust constraint potentials to reflect the spin shift
                  if(n_active_regions.gt.1) then
                     ! in this case, an initial shift makes sense and must be readded here.
                     do i_region = 1,n_active_regions,1
                        constraint_potential(i_region,i_spin) = &
                           constraint_potential(i_region,i_spin) - spin_shift(i_spin)
                     enddo
                  else
                     ! plain fixed-spin-moment constraint - initial shift was not applied.
                     constraint_potential(1,i_spin) = - spin_shift(i_spin)
                  endif

                  ! adjust the final eigenvalues to reflect the spin shift
                  do i_state = 1,n_states,1
                     KS_eigenvalue(i_state,i_spin,1) = KS_eigenvalue(i_state,i_spin,1) - &
                                                       spin_shift(i_spin)
                  enddo

                  chemical_potential_spin(i_spin) = chemical_potential_spin(i_spin) - &
                                                    spin_shift(i_spin)

               enddo

               ! Both spin channel chemical potentials must now be equal
               ! simply by construction!
               chemical_potential = chemical_potential_spin(1)

               ! end adjustment of spin channels if needed
            endif

            ! deallocations
            if(allocated(arguments)) deallocate(arguments)
            if(allocated(gradient))  deallocate(gradient)
            if(allocated(hessian))   deallocate(hessian)
            if(allocated(work))      deallocate(work)
            if(allocated(xm))        deallocate(xm)

            ! Finally, check for convergence:
            constraint_converged = .true.

            do i_spin = 1,n_spin,1
               do i_region = 1,n_active_regions,1
                  constraint_converged = constraint_converged .and. &
                     (dabs(electrons_in_region(i_region,i_spin)- &
                     constraint_electrons(i_region,i_spin)).lt. &
                     constraint_precision)
               enddo
            enddo

            ! end if (mixer_constraint.eq.3)
         endif

         ! after end of loop, compute average auxiliary potential for all subsystems and shift all eigenvalues,
         ! constraint_potentials uniformly so as to reduce the energy correction to zero.

         avg_potential = 0.d0
         do i_spin = 1,n_spin,1
            do i_region = 1,n_active_regions,1
               avg_potential = avg_potential + &
                  electrons_in_region(i_region,i_spin) * &
                  constraint_potential(i_region,i_spin)
            enddo
         enddo
         avg_potential = avg_potential / n_electrons

         do i_region = 1,n_active_regions,1
            do i_spin = 1,n_spin,1
               constraint_potential(i_region,i_spin) = &
                  constraint_potential(i_region,i_spin) - avg_potential
            enddo
         enddo

         do i_spin = 1,n_spin,1
            do i_state = 1,n_states,1
               KS_eigenvalue(i_state,i_spin,1) = &
                  KS_eigenvalue(i_state,i_spin,1) - avg_potential
            enddo
         enddo

         chemical_potential = chemical_potential - avg_potential

         ! Write all important results ...
         write(info_str,*)
         call localorb_info(info_str,use_unit,'(A)')

         ! issue warning if the inner scf loop was exited without convergence.
         if(constraint_converged.and.(mixer_constraint.le.1)) then
            write(info_str,'(2X,A,I5,A)') "Locally constrained DFT: Converged after ", &
                  number_constraint_iter, " inner scf iterations."
            call localorb_info(info_str,use_unit,'(A)')

         elseif(constraint_converged) then
            do i_spin = 1,n_spin,1
               write(info_str,'(2X,A,I5,A,I5,A)')"Locally constrained DFT: Spin ", i_spin, &
                     " converged after ", n_eval(i_spin), " iterations."
               call localorb_info(info_str,use_unit,'(A)')
            enddo
         elseif(mixer_constraint.eq.2) then
            write(info_str,'(2X,A)')"Locally constrained DFT: NO CONVERGENCE."
            call localorb_info(info_str,use_unit,'(A)')

            do i_spin = 1,n_spin,1
               write(info_str,'(2X,A,I5,A,I5,A)') &
                     "Spin ", i_spin, &
                     ": ", &
                     n_eval(i_spin), " iterations."
               call localorb_info(info_str,use_unit,'(A)')
            enddo
         else
            write(info_str,'(2X,A,I5,A)') &
                  "Locally constrained DFT: NO CONVERGENCE after ", &
                  number_constraint_iter, " inner scf iterations."
            call localorb_info(info_str,use_unit,'(A)')
         endif

         write(info_str,'(2X,A)') "| "
         call localorb_info(info_str,use_unit,'(A)')

         write(info_str,'(2X,A,A)') "| Final constraint potentials in different subsystems:"
         call localorb_info(info_str,use_unit,'(A)')

         write(info_str,'(2X,A,10X,A,I2,12X,A,I2)')"|   Region  ", (" Spin ", i_spin, i_spin=1,n_spin,1)
         call localorb_info(info_str,use_unit,'(A)')

         do i_region = 1,n_active_regions,1

            write(info_str,'(2X,A,I6,2X,F15.10,A,2X,F15.10,A)')"|   ", constraint_region_label(i_region), &
                 (constraint_potential(i_region,i_spin)*hartree, " eV", i_spin=1,n_spin,1)
            call localorb_info(info_str,use_unit,'(A)')
         enddo

         write(info_str,'(2X,A)') "| "
         call localorb_info(info_str,use_unit,'(A)')

         write(info_str,'(2X,A)') "| Electron count in different subsystems:"
         call localorb_info(info_str,use_unit,'(A)')

         write(info_str,'(2X,A,10X,A,I2,13X,A,I2)')"|   Region   ", ("Spin ", i_spin, i_spin=1,n_spin,1)
         call localorb_info(info_str,use_unit,'(A)')

         do i_region = 1,n_active_regions,1

            write(info_str,'(2X,A,I6,5X,F15.10,5X,F15.10)')"|   ", constraint_region_label(i_region), &
                 (electrons_in_region(i_region,i_spin), i_spin=1,n_spin,1)
            call localorb_info(info_str,use_unit,'(A)')
         enddo

         write(info_str,'(2X,A)') "| "
         call localorb_info(info_str,use_unit,'(A)')

         write(info_str,'(2X,A,A)')  "| Final average constraint potential:"
         call localorb_info(info_str,use_unit,'(A)')

         write(info_str,'(2X,A,F22.10,A,F22.10,A)') "|   ", avg_potential, " Ha   ", avg_potential*hartree, " eV."
         call localorb_info(info_str,use_unit,'(A)')

         write(info_str,'(2X,A)') "| "
         call localorb_info(info_str,use_unit,'(A)')

         write(info_str,'(2X,A,A)')  "| Overall Fermi level:"
         call localorb_info(info_str,use_unit,'(A)')

         write(info_str,'(2X,A,F22.10,A,F22.10,A)')"|   ", chemical_potential, " Ha   ", chemical_potential*hartree, " eV."
         call localorb_info(info_str,use_unit,'(A)')

         write(info_str,*)
         call localorb_info(info_str,use_unit,'(A)')

         if(mixer_constraint.eq.1) then
            call cleanup_pulay_constraint( )
         endif

         do i_spin = 1,n_spin,1
            do i_region = 1,n_active_regions,1
               do i_state = 1,n_states,1
                  constraint_proj_final(i_state, i_region, i_spin) = &
                     constraint_proj(i_state, i_region, i_spin)
               enddo
            enddo
         enddo
         ! endif relates to everything to do with constrained DFT.
      endif

   else ! (.not. ((force_occupation_projector).and.start_force_occ))
      call get_KS_orbitals_forced_occ_p0(overlap_matrix, hamiltonian, &
              n_electrons, KS_eigenvalue, KS_eigenvector, &
              KS_eigenvector_complex, occ_numbers, chemical_potential)
   endif

end subroutine advance_KS_solution
!******
