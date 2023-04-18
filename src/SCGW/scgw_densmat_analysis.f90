      subroutine scgw_densmat_analysis ()
      
      use dimensions
      use runtime_choices
      use species_data
      use physics
      use prodbas
      use hartree_fock
      use gw_para
      use constants
      use mpi_tasks
      use synchronize_mpi
      use numerical_utilities 
      use scgw_grid
      use poles_fit
      use scgw_allocations 
      use gt
      use localorb_io, only: use_unit

      implicit none      
      real*8       gev(n_basis)

!------------------------------------------------------------------
        if(.true.)then !diagonalize the density matrix
          if(myid.eq.0) write(use_unit,*) " --- Start Density Matrix ", &
                  "analysis --- "
 
          !get the density matrix from G
          if(.not. allocated (green_fn))then
             allocate(green_fn(n_basis,n_basis))
          endif
          green_fn(:,:) = 0.d0
          if (n_spin .eq.1)then
            green_fn(:,:) = green_fn_time(:,:,0,1)*2.d0
          elseif (n_spin .eq.2)then
            green_fn(:,:) = green_fn_time(:,:,0,1) + green_fn_time(:,:,0,2)
          endif

          !get the full overlap matrix
!          if(.not. allocated (full_hamiltonian))then
!             allocate(full_hamiltonian(n_basis,n_basis))
!          endif 
          if(.not. allocated (aux_ovlp))then
             allocate(aux_ovlp(n_basis,n_basis))
          endif
 
          aux_ovlp(:,:) = 0.d0
          call invert_overlap_matrix_2 (aux_ovlp)

!          i_index = 0
!          full_hamiltonian(:,:) = 0.d0
!          do j_basis = 1, n_basis, 1
!            do i_basis = 1, j_basis, 1
!             i_index = i_index+1
!             full_hamiltonian (j_basis,i_basis) = hamiltonian (i_index,1)
!             full_hamiltonian (i_basis,j_basis) = full_hamiltonian (j_basis,i_basis)
!            enddo
!          enddo

          if(.not.allocated(green_tmp))then
             allocate(green_tmp (n_basis,n_basis))
          endif
           green_tmp (:,:) = 0.d0
          call dgemm ('N','N',n_basis,n_basis, n_basis,&
             1.d0, aux_ovlp, n_basis, &
              green_fn , n_basis, 0.d0,&
             green_tmp, n_basis)
           green_fn(:,:) = 0.d0
          call dgemm ('N','N',n_basis,n_basis, n_basis,&
             1.d0, green_tmp , n_basis, &
             aux_ovlp, n_basis, 0.d0, &
             green_fn, n_basis)

!Test of the diagonalization of the Hamiltonian (it just obtains the KS eigenvalues if correct)
!           green_tmp (:,:) = 0.d0
!          call dgemm ('N','N',n_basis,n_basis, n_basis,&
!             1.d0, aux_ovlp, n_basis, &
!              full_hamiltonian, n_basis, 0.d0,&
!             green_tmp, n_basis)
!!          full_hamiltonian(:,:)= green_tmp(:,:)
!          call dgemm ('N','N',n_basis,n_basis, n_basis,&
!             1.d0, green_tmp , n_basis, &
!             aux_ovlp, n_basis, 0.d0, &
!             full_hamiltonian, n_basis)
          if(allocated(green_tmp))then
             deallocate(green_tmp)
          endif

          gev(:) = 0.d0
          call diagonalize_rmatrix (n_basis, green_fn, gev, .false.)
!          call diagonalize_rmatrix (n_basis, full_hamiltonian, gev, .false.)
         if(myid.eq.0)then 
            print *, "      -> Occupation of the natural orbitals <-"  
            do i_basis =1, n_basis, 1
              print *,  i_basis, gev (i_basis) !*hartree, KS_eigenvalue (i_basis, 1,1)*hartree
            enddo
          endif
          if(allocated (green_fn))then
             deallocate(green_fn)
          endif
          if( allocated (aux_ovlp))then
             deallocate(aux_ovlp)
          endif
 
          if(myid.eq.0) write(use_unit,*) " --- Finish Density Matrix ", &
                  "analysis --- "
        endif


      end subroutine scgw_densmat_analysis
