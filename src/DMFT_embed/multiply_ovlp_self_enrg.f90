  subroutine multiply_ovlp_self_enrg(self_enrg_gw, &
                                     full_ovlp_matrix_sqrt, &
                                     k_summed_overlap_matr, &
                                     self_enrg_gw_k)


  use dimensions
  use runtime_choices
  use species_data
  use pbc_lists
  use physics
  use prodbas
  use scgw_grid
  use constants
  use mpi_tasks
  use synchronize_mpi
  use hartree_fock
  use localized_basbas
  use gw_para
  use dmft_para
  use poles_fit 
  use gt
  use timing


 implicit none


 integer i_freq, i_k_point, i_basis, j_basis

! INPUT---------------------------------------------------------------------
 complex*16, dimension(n_basis,n_basis,n_k_points) :: full_ovlp_matrix_sqrt
 complex*16, dimension(n_basis,n_basis,nomega) :: self_enrg_gw
 complex*16, dimension(:,:,:), allocatable :: self_enrg_gw_temp

! OUTPUT---------------------------------------------------------------------
 complex*16, dimension(n_basis,n_basis,nomega,n_k_points) :: self_enrg_gw_k

! ARRAYS defined in this subroutine-----------------------------------------
 complex*16, allocatable :: aux_matr(:,:)
 !real*8, dimension(n_basis,n_basis) :: k_summed_overlap_matr
 complex*16, dimension(n_basis,n_basis) :: k_summed_overlap_matr


       if (.not. allocated(self_enrg_gw_temp))then
          allocate(self_enrg_gw_temp(n_basis,n_basis, nomega))
       endif
!do i_basis = 1, n_basis,1
!write(use_unit,*) 'self_energy GW', self_enrg_gw(i_basis,i_basis,1)
!enddo


     do i_freq = 1, nomega,1
       if (.not. allocated(aux_matr))then
          allocate(aux_matr(n_basis,n_basis))
       endif

 self_enrg_gw_temp(:,:,i_freq) = (0.d0,0.d0)
 aux_matr(:,:) = (0.d0,0.d0)

           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            k_summed_overlap_matr(:,:), n_basis, &
            self_enrg_gw(:,:,i_freq), &
            n_basis,(0.d0,0.d0), &
            aux_matr(:,:) , n_basis )

           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            aux_matr(:,:), n_basis, &
            k_summed_overlap_matr(:,:), &
            n_basis,(0.d0,0.d0), &
            self_enrg_gw_temp(:,:,i_freq) , n_basis )

       if ( allocated(aux_matr))then
          deallocate(aux_matr)
       endif

!self_enrg_gw(:,:,i_freq) = self_enrg_gw_temp(:,:,i_freq)
  enddo
!self_enrg_gw(:,:,:) = self_enrg_gw_temp(:,:,:)


 self_enrg_gw_k(:,:,:,:) = (0.d0,0.d0)


 do i_k_point = 1, n_k_points,1
    do i_freq = 1, nomega,1 

       if (.not. allocated(aux_matr))then
          allocate(aux_matr(n_basis,n_basis))
       endif

        aux_matr(:,:) = (0.d0,0.d0)


           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            full_ovlp_matrix_sqrt(:,:,i_k_point), n_basis, &
!            self_enrg_gw_temp(:,:,i_freq), &
            self_enrg_gw(:,:,i_freq), &
            !inv_omega_ovlp(:,:,i_freq,i_k_point), &
            n_basis,(0.d0,0.d0), &
            aux_matr(:,:) , n_basis )

           call zgemm ('N','N', n_basis,&
            n_basis,n_basis,(1.d0,0.d0),&
            aux_matr(:,:), n_basis, &
            full_ovlp_matrix_sqrt(:,:,i_k_point), &
            n_basis,(0.d0,0.d0), &
            self_enrg_gw_k(:,:,i_freq,i_k_point) , n_basis )


       if ( allocated(aux_matr))then
          deallocate(aux_matr)
       endif


    enddo
    enddo

end subroutine  multiply_ovlp_self_enrg
