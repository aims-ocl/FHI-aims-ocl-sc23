           subroutine transform_to_KS_basis_complex(NAO_matr,&
                                                    i_k_point,&
                                                    KS_matr)






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

      integer i_k_point
      complex*16, dimension(n_basis,n_basis) :: NAO_matr
      complex*16, dimension(n_states,n_states) :: KS_matr 
      complex*16,dimension(:,:), allocatable :: aux_matr


        KS_matr(:,:) = (0.d0,0.d0)

       if (.not. allocated(aux_matr))then
          allocate(aux_matr(n_basis,n_states))
       endif

        aux_matr(:,:) = (0.d0,0.d0)

      if(real_eigenvectors)then

        call zgemm ('N','N',n_basis,n_states, n_basis,&
            (1.d0,0.d0),NAO_matr(:,:), n_basis, &
            !DCMPLX(KS_eigenvector(:,:,1,i_k_point)), n_basis, (0.d0,0.d0),&
            DCMPLX(KS_eigenvector(:,:,1,1)), n_basis, (0.d0,0.d0),&
            aux_matr(:,:),n_basis)


        call zgemm ('T','N',n_states,n_states, n_basis,&
            !(1.d0,0.d0),DCMPLX(KS_eigenvector(:,:,1,i_k_point)) , n_basis, &
            (1.d0,0.d0),DCMPLX(KS_eigenvector(:,:,1,1)) , n_basis, &
            aux_matr(:,:), n_basis, (0.d0,0.d0), &
            KS_matr(:,:) , n_states)

       else

        call zgemm ('N','N',n_basis,n_states, n_basis,&
            (1.d0,0.d0),NAO_matr(:,:), n_basis, &
            KS_eigenvector_complex(:,:,1,i_k_point), n_basis, (0.d0,0.d0),&
            aux_matr(:,:),n_basis)


        call zgemm ('C','N',n_states,n_states, n_basis,&
            (1.d0,0.d0),KS_eigenvector_complex(:,:,1,i_k_point) , n_basis, &
            aux_matr(:,:), n_basis, (0.d0,0.d0), &
            KS_matr(:,:) , n_states)

       endif
       if (allocated(aux_matr))then
          deallocate(aux_matr)
       endif


           end subroutine transform_to_KS_basis_complex
