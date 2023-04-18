!****s* FHI-aims/extrapolar
!  NAME
!   get_green_freq 
!  SYNOPSIS

      subroutine  evaluate_ovlp_NAO_KS (&
         ovlp_NAO_KS)


      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
!      use hartree_fock
      use physics
!      use gw_para

      implicit none
!ARGUMENTS
      real*8 ovlp_NAO_KS (n_states, n_basis, n_spin)
      real*8 aux_overlap_matrix (n_basis, n_basis) 

!COUNTERS 
      integer i_state
      integer i_basis
      integer j_basis
      integer i_index
      integer i_spin



      i_index = 0
      ovlp_NAO_KS(:,:,:) = 0.d0
      aux_overlap_matrix(:,:) = 0.d0 

!      call sync_integrate_ovlp( &
!           overlap_matrix )
       
!      print *, overlap_matrix
      !unpack the overlap matrix 
      do i_basis = 1, n_basis, 1
        do j_basis= 1, i_basis, 1
      
         i_index= i_index+1
         aux_overlap_matrix (i_basis, j_basis) = &
         overlap_matrix(i_index)

         if(i_basis.ne.j_basis)then
           aux_overlap_matrix (j_basis, i_basis) = &
           aux_overlap_matrix (i_basis, j_basis)
         endif

       enddo
      enddo 
 
       
      do i_spin = 1, n_spin, 1
        do i_basis = 1, n_basis, 1
          do j_basis= 1, n_basis, 1
            do i_state= 1, n_states, 1
      
             ovlp_NAO_KS(i_state, i_basis,i_spin) = &
             ovlp_NAO_KS(i_state, i_basis,i_spin) +&
             KS_eigenvector(j_basis, i_state, i_spin, 1)*&
             aux_overlap_matrix(i_basis, j_basis)
      
           enddo
         enddo
        enddo 
      enddo

      return
      end subroutine evaluate_ovlp_NAO_KS
