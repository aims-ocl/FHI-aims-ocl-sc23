           subroutine multiply_ovlp_matr(matr_in,&
                                         ovlp_matr,&
                                         matr_out)






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

      complex*16, dimension(n_basis,n_basis) :: matr_in
      complex*16, dimension(n_basis,n_basis) :: matr_out
      complex*16, dimension(n_basis,n_basis) :: ovlp_matr
      complex*16,dimension(:,:), allocatable :: aux_matr


        matr_out(:,:) = (0.d0,0.d0)

       if (.not. allocated(aux_matr))then
          allocate(aux_matr(n_basis,n_basis))
       endif

        aux_matr(:,:) = (0.d0,0.d0)

        call zgemm ('N','N',n_basis,n_basis, n_basis,&
            (1.d0,0.d0),ovlp_matr(:,:), n_basis, &
            matr_in(:,:), n_basis, (0.d0,0.d0),&
            aux_matr(:,:),n_basis)


        call zgemm ('N','N',n_basis,n_basis, n_basis,&
            (1.d0,0.d0),aux_matr(:,:) , n_basis, &
            ovlp_matr(:,:), n_basis, (0.d0,0.d0), &
            matr_out(:,:) , n_basis)

       if (allocated(aux_matr))then
          deallocate(aux_matr)
       endif


           end subroutine multiply_ovlp_matr
