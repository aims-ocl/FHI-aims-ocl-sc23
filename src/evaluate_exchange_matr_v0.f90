!****s* FHI-aims/evaluate_exchange_matr_v0
!  NAME
!   evaluate_exchange_matr_v0
!  SYNOPSIS

      subroutine evaluate_exchange_matr_v0 &
           ( KS_eigenvector, number_of_loops, &
              occ_numbers)

! PURPOSE
! construct the exchange matrix, in this version of the construction density matrix
! is used. This reduces the memory consumption, but slows down the calculation 
! in particular we the occupied states are few.
! 
! USES

      use dimensions
      use prodbas
      use hartree_fock
      use runtime_choices
      use mpi_tasks
      use synchronize_mpi
      use constants
      use localorb_io, only: use_unit

      implicit none

! ARGUMENTS

      integer :: number_of_loops
      real*8, dimension(n_basis,n_states,n_spin,n_k_points) :: &
                        KS_eigenvector
      real*8, dimension(n_states,n_spin,n_k_points) :: occ_numbers

!  INPUTS
!  o  number_of_loops -- integer number, the current number of the self-consistent loop
!  o  occ_numbers -- real array,
!       the occupation number of the electrons for each eigenstate and each spin
!  o  KS_eigenvector -- real array,
!            the eigenvector of the single-particle (KS/HF) self-consistent calculation
!  OUTPUTS
!  o  none
!     the exchange matrix (the "fock_matr" in the source code) is defined in MODULE
!     hartree_fock
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

!     auxiliary matrices for Level 3 Blas matrix multiplications

       real*8, dimension(:,:), allocatable ::  aux_eigenvector
       real*8, dimension(:,:), allocatable ::  density_matrix
       real*8, dimension(:,:), allocatable ::  ovlp_denmat
       real*8, dimension(:,:), allocatable ::  temp_ovlp_matr

!     counters

      integer :: i_state

      integer i_spin
      integer i_basis_1
      integer i_basis_2
      integer i_basis_3
      integer i_index

!     begin work


      if(myid.eq.0) then
        write(use_unit,'(2X,A)')"Evaluating the exchange matrix ... "
      endif

      allocate(aux_eigenvector(n_basis,n_homo_max))
      allocate(density_matrix(n_basis,n_basis))
      allocate(temp_ovlp_matr(n_basis,n_basis))
      allocate(ovlp_denmat(n_basis, n_basis))

      fock_matr(:,:,:) = 0.d0
      if (use_lc_wpbeh .and. hybrid_coeff .ne. 0.d0) then
      	fock_matr_SR(:,:,:) = 0.d0
      	if(myid.eq.0) then
		     write(use_unit,'(4X,A)')"for LR ... "
		   endif
      end if
		do i_spin = 1, n_spin
		  do i_state = 1, n_homo_max, 1
			  do i_basis_1 = 1, n_basis, 1
			    aux_eigenvector(i_basis_1, i_state) = &
			     KS_eigenvector(i_basis_1, i_state, i_spin, 1) * &
			     occ_numbers(i_state, i_spin, 1)*dble(n_spin)/2.d0
			  enddo
	 	  enddo

		  density_matrix (:,:) = 0
		  call dgemm('N', 'T', n_basis, n_basis, &
		       n_homo(i_spin), 1.d0, &
		       KS_eigenvector(:,:,i_spin,1), n_basis, &
		       aux_eigenvector(:,1:n_homo(i_spin)), n_basis, &
		       0.d0, &
		       density_matrix(:,:), n_basis &
		      )

		  if(use_density_matrix_hf .and. number_of_loops >= 0) then
		    call density_matrix_mixing &
		         (i_spin, number_of_loops, density_matrix)
		  endif
      
        do i_basis_1 = 1, n_loc_prodbas, 1
          i_index = 0
          temp_ovlp_matr(:,:) = 0.d0
          do i_basis_2 = 1, n_basis, 1
            do i_basis_3 = 1, i_basis_2, 1

              i_index = basis_nghr(i_basis_3, i_basis_2)
              if(i_index .gt. 0) then
	              temp_ovlp_matr(i_basis_3, i_basis_2) = &
                 		 ovlp_3fn(i_index, i_basis_1)
              endif


             enddo
           enddo

           ovlp_denmat(:,:) = 0.d0
           call dsymm('L', 'U', n_basis, n_basis, &
                       1.d0, &
                       temp_ovlp_matr, n_basis, &
                       density_matrix, n_basis, &
                       0.d0, &
                       ovlp_denmat, n_basis &
                      )

           call dsymm('R', 'U', n_basis, n_basis, &
                       1.d0, &
                       temp_ovlp_matr, n_basis, &
                       ovlp_denmat, n_basis, &
                       1.d0, &
                       fock_matr(:,:,i_spin), n_basis &
                      )
            
			  if (use_lc_wpbeh .and. hybrid_coeff .ne. 0.d0) then
			  	 if(myid.eq.0 .and. i_basis_1 .eq. 1 .and. i_spin .eq. 1) then
				   write(use_unit,'(4X,A)')"for SR ... "
				 endif
				 temp_ovlp_matr(:,:) = 0.d0
		       do i_basis_2 = 1, n_basis, 1
		         do i_basis_3 = 1, i_basis_2, 1
		           i_index = basis_nghr(i_basis_3, i_basis_2)
		           if(i_index .gt. 0) then
			           temp_ovlp_matr(i_basis_3, i_basis_2) = &
		              		 ovlp_3fn_SR(i_index, i_basis_1)
		           endif
		         enddo
		       enddo

		       ovlp_denmat(:,:) = 0.d0
		       call dsymm('L', 'U', n_basis, n_basis, &
		                   1.d0, &
		                   temp_ovlp_matr, n_basis, &
		                   density_matrix, n_basis, &
		                   0.d0, &
		                   ovlp_denmat, n_basis &
		                  )

		       call dsymm('R', 'U', n_basis, n_basis, &
		                   1.d0, &
		                   temp_ovlp_matr, n_basis, &
		                   ovlp_denmat, n_basis, &
		                   1.d0, &
		                   fock_matr_SR(:,:,i_spin), n_basis &
		                  )
			  end if
!   end of loop over i_basis_1
        enddo

       if(use_mpi) then
         call sync_matrix(fock_matr(:,:,i_spin), &
                   n_basis, n_basis)
         if (use_lc_wpbeh .and. hybrid_coeff .ne. 0.d0) then
         	call sync_matrix(fock_matr_SR(:,:,i_spin), &
                   n_basis, n_basis)
         end if
       endif

!       if(myid.eq.0) then
!       write(use_unit,*) "fock_matr: "
!       do i_basis_1 = 1, n_basis, 1
!         do i_basis_2 = 1, i_basis_1, 1
!             write(use_unit,'(2I6,f18.10)') i_basis_1, i_basis_2, &
!               fock_matr(i_basis_2, i_basis_1, i_spin)
!         enddo
!       enddo
!       endif
!! loop over i_spin
      enddo

      if(allocated(aux_eigenvector)) then
         deallocate(aux_eigenvector)
      endif
      if(allocated(density_matrix)) then
         deallocate(density_matrix)
      endif
      if(allocated(ovlp_denmat)) then
         deallocate(ovlp_denmat)
      endif
      if(allocated(temp_ovlp_matr)) then
         deallocate(temp_ovlp_matr)
      endif

      end subroutine evaluate_exchange_matr_v0
!---------------------------------------------------------------------
!******
