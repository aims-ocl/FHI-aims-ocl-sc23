  subroutine get_NEW_embed_gf_freq(free_cluster_ham,&
                                   !k_summed_overlap_matr,&
                                   !inv_k_summed_overlap_matr,&
                                   hybrid_func,loc_self_enrg,&
                                   hartree_pot_LDA,&
                                   embed_gf, omega_term,&
                                   number_reiterations) ! (KS_eigenvalue)

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


! local parameters
   integer a, i_a, number_reiterations 
   integer i_matrix_size
   integer j_matrix_size
   integer j_basis
   integer i_basis, i_basis_1
   integer i_freq 
   integer i_k_point
   character*2 iter
   character*2 iter1
   character*17 filename_RE
   character*17 filename_IM
 
! subroutine actual parameters
   real*8, dimension(n_basis,n_basis), intent(in) :: free_cluster_ham
   real*8, dimension(n_basis,n_basis), intent(in) :: loc_self_enrg 
   complex*16, dimension(n_basis,n_basis,nomega) :: hybrid_func
!   real*8, dimension(n_basis,n_basis), intent(in) :: k_summed_overlap_matr
!   real*8, dimension(n_basis,n_basis), intent(in) :: inv_k_summed_overlap_matr
   complex*16, dimension(:,:,:), allocatable :: inv_embed_gf
   complex*16, dimension(n_basis,n_basis,nomega), intent(out) :: embed_gf
   !complex*16, dimension(n_basis,n_basis,nomega), intent(in) :: hybrid_func_temp
!   complex*16, dimension(:,:), allocatable :: inv_k_summed_overlap_matr


!   complex*16, dimension(:,:), allocatable:: inv_k_summed_overlap_matr_sqrt
!   complex*16, dimension(:,:), allocatable :: k_summed_overlap_matr_sqrt
!   complex*16, dimension(:,:,:), allocatable :: ovlp_omega
!   complex*16, dimension(:,:,:), allocatable :: ovlp_omega_ovlp
   complex*16, dimension(n_basis, n_basis, nomega):: omega_term
   real*8, dimension(n_basis,n_basis), intent(in) :: hartree_pot_LDA


   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 


   logical  :: output

 
!     if (.not. allocated (inv_k_summed_overlap_matr_sqrt)) then
!         allocate (inv_k_summed_overlap_matr_sqrt(n_basis,n_basis))
!     endif
!     if (.not. allocated (inv_k_summed_overlap_matr)) then
!         allocate (inv_k_summed_overlap_matr(n_basis,n_basis))
!     endif
!     if (.not. allocated (k_summed_overlap_matr_sqrt)) then
!         allocate (k_summed_overlap_matr_sqrt(n_basis,n_basis))
!     endif
     if (.not. allocated (inv_embed_gf)) then
         allocate (inv_embed_gf(n_basis,n_basis,nomega))
     endif




       if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif



!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------




        inv_embed_gf(:,:,:) = (0.d0,0.d0)
       do i_freq = 1, nomega, 1
          do i_matrix_size = 1, n_basis, 1
            do j_matrix_size = 1, n_basis, 1


                  inv_embed_gf(i_matrix_size,j_matrix_size,i_freq) = &
                  inv_embed_gf(i_matrix_size,j_matrix_size,i_freq) + &
!                  ((ovlp_omega_ovlp(i_matrix_size,j_matrix_size,i_freq))- &
                  ((omega_term(i_matrix_size,j_matrix_size,i_freq))- &
!                 ((((0.d0,1.d0)*omega(i_freq))+chemical_potential) - & 
                  free_cluster_ham(i_matrix_size,j_matrix_size)- &
                  !hartree_pot_LDA(i_matrix_size,j_matrix_size)- &
                  loc_self_enrg(i_matrix_size,j_matrix_size)- &
                  hybrid_func(i_matrix_size,j_matrix_size, i_freq))
!                  ovlp_hybrid_func(i_matrix_size,j_matrix_size, i_freq))


           enddo
          enddo
         enddo


      do i_basis_1 = 1, n_basis, 1
!             write(use_unit,*)'hartree_pot_LDA',  i_basis_1, &
!               hartree_pot_LDA(i_basis_1, i_basis_1)
      enddo





! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------






       !embed_gf(:,:,:) = inv_embed_gf(:,:,:)
       embed_gf(:,:,:) = inv_embed_gf(:,:,:)

! at each frequency point evaluate the LU factorization, and check for errors
     do i_freq = 1, nomega, 1

         !call zgetrf( n_basis, n_basis, embed_gf_temp(:,:,i_freq), &
         call zgetrf( n_basis, n_basis, embed_gf(:,:,i_freq), &
                    n_basis , ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then
            write(use_unit,*) " * Failure of LU decomposition at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
            stop
          endif
        endif
! this inverts the matrix in the LU factorized form
        !call zgetri( n_basis, embed_gf_temp(:,:,i_freq), n_basis, &
        call zgetri( n_basis, embed_gf(:,:,i_freq), n_basis, &
                 ipiv, work, n_basis, info)

      if (info.ne.0) then
         if(myid.eq.0)then
            write(use_unit,*) " * Failure of matrix inversion at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
          endif
        endif


      enddo

!      do i_basis_1 = 1, n_basis, 1
!            write(use_unit,*)'omega_term from NEW_gf_embed',  i_basis_1, &
!             omega_term(i_basis_1, i_basis_1,100)
!      enddo

!      do i_basis_1 = 1, n_basis, 1
!            write(use_unit,*)'embed_gf from NEW_gf_embed',  i_basis_1, &
!             embed_gf(i_basis_1, i_basis_1,100)
!      enddo

!      do i_basis_1 = 1, n_basis, 1
!            write(use_unit,*)'loc_self_enrg from NEW_gf_embed',  i_basis_1, &
!             loc_self_enrg(i_basis_1, i_basis_1)
!      enddo

!stop
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------




!endif !myid.eq.0
!  end subroutine get_NEW_embed_gf_freq 
!      do i_basis_1 = 1, n_basis, 1
!            write(use_unit,*)'embed_gf from NEW_gf_embed',  i_basis_1, &
!             embed_gf(i_basis_1, i_basis_1,100)
!      enddo

!stop
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------




!endif !myid.eq.0
  end subroutine get_NEW_embed_gf_freq 
