!  NAME
!    get_kinetic_energy - computes T through the kinetic energy matrix 
!  SYNOPSIS
      subroutine get_kinetic_energy &
       (rho, &
       partition_tab, basis_l_max, &
       KS_eigenvector, KS_eigenvector_complex, &
       occ_numbers, kinetic_energy &
       )
!  PURPOSE
!    Computes the matrix elements of the KS kinetic energy, computes the density
!    matrix, multiplies them and finally computes the kinetic energy. This is a
!    post-processing subroutine, called in scf_solver.f90. 
!  USES
      use dimensions
      use runtime_choices
      use grids
      use geometry
      use basis
      use constants
      use localorb_io
      use constants
      use density_matrix_evaluation
      use energy_density
      use sym_base, only: evaluate_densmat_sym
      use mpi_tasks, only: myid, check_allocation
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    kinetic_energy 
!  SEE ALSO
!    FHI-aims CPC publication (see copyright notice below)
!  COPYRIGHT
!    Copyright by FHI-aims Team,
!    Theory Department, Fritz Haber Institute, Berlin (2008).
!
!    FHI-aims comes without warranty of any kind. We hope that it will
!    be
!    useful for your work, and if you contact us, will try to help solve
!    any problems as best as we can.
!
!    FHI-aims has been distributed to you including its source code, 
!    but the terms of the license agreement apply. In particular, do 
!    not redistribute FHI-aims. Instead, please refer all interested 
!    parties to
!
!      aims-coordinators@fhi-berlin.mpg.de
!
!    instead. If you use FHI-aims, cite
!
!      Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!      Xinguo Ren, Karsten Reuter, and Matthias Scheffler, 
!      "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!      Computer Physics Communications (2008), submitted.
!
!    (or preferably the final version of this article.)
!******
!  Any other comments (not for documentation) follow below.

!  Declaration of variables
      implicit none

     ! output
       real*8, dimension(n_spin) :: kinetic_energy

     ! imported variables
      real*8, dimension(n_spin, n_full_points) :: rho
      real*8, dimension( n_full_points) :: partition_tab
      integer :: basis_l_max (n_species)
      real*8,dimension(n_basis, n_states, n_spin,n_k_points_task) :: KS_eigenvector
      complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task):: KS_eigenvector_complex
      real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
      real*8, external :: ddot


     ! local variables
      integer :: i_point
      integer :: i_spin
      integer :: info
      real*8, dimension(n_hamiltonian_matrix_size, n_spin) :: kinetic
      real*8, dimension(n_hamiltonian_matrix_size) :: density_matrix_sparse
      real*8, dimension(n_hamiltonian_matrix_size, n_spin) :: density_matrix_sparse_tmp
      real*8, dimension(:,:), allocatable :: density_matrix
      real*8, dimension(:,:), allocatable :: density_matrix_con
      real*8, dimension(:,:), allocatable :: kinetic_con
      real*8, dimension(n_hamiltonian_matrix_size, n_spin) :: density_matrix_tmp
      integer :: i_basis_1, i_basis_2, i_index

       !allocations
       if(.not. allocated(density_matrix))then
        allocate(density_matrix(n_centers_basis_T, n_centers_basis_T),stat=info)
        call check_allocation(info, 'density_matrix                ') 
       end if
       if(.not. allocated(density_matrix_con))then
        allocate(density_matrix_con(n_centers_basis_T, n_centers_basis_T),stat=info)
        call check_allocation(info, 'density_matrix_con            ') 
       end if
       if(.not. allocated(kinetic_con))then
        allocate(kinetic_con(n_centers_basis_T, n_centers_basis_T),stat=info)
        call check_allocation(info, 'kinetic_con            ') 
       end if

      ! start calculation
      kinetic_energy = 0.d0
      density_matrix_sparse_tmp(:,:) = 0.d0
      density_matrix_tmp(:,:)=0.d0
      density_matrix_con(:,:)=0.d0
      kinetic_con(:,:)=0.d0

      ! Calculate <phi_i | T | phi_j> - spin loop is inside already
      call integrate_real_kinetic_matrix &
               (rho,  &
                partition_tab, basis_l_max, &
                kinetic &
               )

      do i_spin=1, n_spin, 1
           ! Calculate  sum_l (occ_numbers  c_il  c_jl)
	   if(use_symmetry_reduced_spg)then
	     call evaluate_densmat_sym &
                   ( KS_eigenvector, KS_eigenvector_complex, &
                     occ_numbers,  &
                     density_matrix, density_matrix_sparse, i_spin, .false. &
                     )
	   else
	     call evaluate_densmat &
                   ( KS_eigenvector, KS_eigenvector_complex, &
                     occ_numbers,  &
                     density_matrix, density_matrix_sparse, i_spin, .false. &
                     )
	   endif
           if (packed_matrix_format /= PM_none) then

              density_matrix_sparse_tmp(1:n_hamiltonian_matrix_size, i_spin) = density_matrix_sparse(1:n_hamiltonian_matrix_size)
              ! reconstruct density matrix
              call unpack_matrices(density_matrix_sparse_tmp(:,i_spin), density_matrix_con, n_centers_basis_T)
              ! reconstruct kinetic
              call unpack_matrices(kinetic(:,i_spin), kinetic_con, n_centers_basis_T)
              ! now calculate the trace... in this case it is the same as calculating the dot_product for each line and summing it... UGH!
              do i_basis_1=1, n_centers_basis_T, 1
                kinetic_energy(i_spin) = kinetic_energy(i_spin)+&
                   dot_product(density_matrix_con(:,i_basis_1), kinetic_con(:,i_basis_1))
              enddo
              ! CC: Do a pointwise assessment of the kinetic energy if
              ! kinetic energy densities are computed 
              if (flag_chetty_martin_energy_density) then
                do i_point=1,n_full_points
                  ! reconstruct kinetic
                  call unpack_matrices(ed_kinetic_batch(:,i_spin,i_point), kinetic_con, n_centers_basis_T)
                  ! now calculate the trace... in this case it is the same as calculating the dot_product for each line and summing it... UGH!
                  do i_basis_1=1, n_centers_basis_T, 1
                    ed_kinetic_energy_density(i_point) = ed_kinetic_energy_density(i_point) + &
                      & dot_product(density_matrix_con(:,i_basis_1), kinetic_con(:,i_basis_1))
                  enddo
                end do
              end if


!              Some line like this shooould work, but obviously I am blas stupid... and it doesn't.
!              kinetic_energy = kinetic_energy + ddot(n_centers_basis_T*n_centers_basis_T, density_matrix_con, 1.d0, kinetic_con, 1.d0)
!DEBUG
!           write(use_unit,*) 'density_matrix_sparse_tmp', density_matrix_sparse_tmp
!           write(use_unit,*) 'kinetic', kinetic
!           write(use_unit,*) 'density_matrix_con', density_matrix_con
!           write(use_unit,*) 'kinetic_con', kinetic_con
!END DEBUG


           else
                i_index = 0
                do i_basis_2 = 1, n_centers_basis_T, 1
                  do i_basis_1 = 1, i_basis_2, 1
                    i_index = i_index+1
                    if(i_basis_1 .eq. i_basis_2) then 
                      density_matrix_tmp(i_index, i_spin) = density_matrix(i_basis_1, i_basis_2)
                    else
                      density_matrix_tmp(i_index, i_spin) = 2.d0*density_matrix(i_basis_1, i_basis_2)
                    endif
                  enddo
                enddo
                kinetic_energy(i_spin) = dot_product(density_matrix_tmp(:,i_spin), kinetic(:,i_spin))
                ! CC: Do a pointwise assessment of the kinetic energy if
                ! kinetic energy densities are computed 
                if (flag_chetty_martin_energy_density) then
                  do i_point=1,n_full_points
                      ed_kinetic_energy_density(i_point) = ed_kinetic_energy_density(i_point) + &
                        & dot_product(density_matrix_tmp(:,i_spin), ed_kinetic_batch(:,i_spin,i_point))
                  end do
                end if
!DEBUG
!           if (myid==0) then
!             write(use_unit,*) 'density_matrix_tmp', density_matrix_tmp
!             write(use_unit,*) 'kinetic', kinetic
!           end if
!END DEBUG
           endif
      end do !i_spin

      if(myid==0) then
        if (n_spin.eq.1) then
         write(use_unit,'(2X,A,1X,F16.5,A,1X,F16.5,A)') &
         "| Kinetic energy of the system     :", &
            kinetic_energy , &
         " Ha", (kinetic_energy)*hartree, " eV" 
        else
         write(use_unit,'(2X,A,1X,F16.5,A,1X,F16.5,A)') &
         "| Kinetic energy of the system spin up      :", &
            kinetic_energy(1) , &
         " Ha", (kinetic_energy(1))*hartree, " eV"  
         write(use_unit,'(2X,A,1X,F16.5,A,1X,F16.5,A)') &
         "| Kinetic energy of the system spin down    :", &
            kinetic_energy(2) , &
         " Ha", (kinetic_energy(2))*hartree, " eV"  
        endif      
      endif

!     Return to initial defaults and deallocate everything
      if(allocated(density_matrix))      deallocate(density_matrix)
      if(allocated(density_matrix_con))      deallocate(density_matrix_con)
      if(allocated(kinetic_con))      deallocate(kinetic_con)

      end subroutine get_kinetic_energy
