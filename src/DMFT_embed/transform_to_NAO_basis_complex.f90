           subroutine transform_to_NAO_basis_complex(KS_matr,&
                                                    i_k_point,&
                                                     NAO_matr,&
                                                     KS_eigenvector_tmp)






     use dimensions
     use runtime_choices
     use species_data
     use pbc_lists
     use physics
     use prodbas
     use constants
     use mpi_tasks
     use synchronize_mpi
     use localized_basbas
     use timing

      implicit none

      integer i_k_point,i_k
      complex*16, dimension(n_basis,n_basis) :: NAO_matr
      complex*16, dimension(n_states,n_states) :: KS_matr 
      complex*16,dimension(:,:), allocatable :: aux_matr
      !complex*16,dimension(:,:,:), allocatable :: KS_eigenvector_tmp
      complex*16,dimension(n_basis,n_states,n_spin) :: KS_eigenvector_tmp

!      if(.not. allocated(KS_eigenvector_tmp))then
!        allocate(KS_eigenvector_tmp(n_basis,n_states,n_spin))
!      endif 

!!         i_k = 0

!           if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
!
!            i_k = i_k + 1
!
!             if(real_eigenvectors)then
!
!               KS_eigenvector_tmp(:,:,:) = DCMPLX(KS_eigenvector(:,:,:,i_k))
!              
!              else 
!         
!              KS_eigenvector_tmp(:,:,:) = KS_eigenvector_complex(:,:,:,i_k)
!
!             endif
!
!           else
!             ! zero temp. KS eigenvector on all other threads, prior to allreduce
!             KS_eigenvector_tmp(:,:,:) = (0.d0,0.d0)
!          end if
!           call sync_eigenvector_complex(KS_eigenvector_tmp)




       NAO_matr(:,:) = (0.d0,0.d0)

       if (.not. allocated(aux_matr))then
          allocate(aux_matr(n_states,n_basis))
       endif

        aux_matr(:,:) = (0.d0,0.d0)

      if(real_eigenvectors)then


           call zgemm ('N','T', n_states,&
            n_basis,n_states,(1.d0,0.d0),&
            KS_matr, n_states, &
            KS_eigenvector_tmp(:,:,1), &
            n_basis,(0.d0,0.d0), &
            aux_matr(:,:) , n_states )

    
           call zgemm ('N','N', n_basis,&
            n_basis,n_states,(1.d0,0.d0),&
            KS_eigenvector_tmp(:,:,1), n_basis, &
            aux_matr(:,:), &
            n_states,(0.d0,0.d0), &
            NAO_matr(:,:) , n_basis )


       else


           call zgemm ('N','C', n_states,&
            n_basis,n_states,(1.d0,0.d0),&
            KS_matr, n_states, &
            KS_eigenvector_tmp(:,:,1), &
            n_basis,(0.d0,0.d0), &
            aux_matr(:,:) , n_states )


           call zgemm ('N','N', n_basis,&
            n_basis,n_states,(1.d0,0.d0),&
            KS_eigenvector_tmp(:,:,1), n_basis, &
            aux_matr(:,:), &
            n_states,(0.d0,0.d0), &
            NAO_matr(:,:) , n_basis )


       endif
       if (allocated(aux_matr))then
          deallocate(aux_matr)
       endif

!      if(allocated(KS_eigenvector_tmp))then
!         deallocate(KS_eigenvector_tmp)
!      endif 

           end subroutine transform_to_NAO_basis_complex
