  subroutine get_embed_gf_freq_no_ovlp(free_cluster_ham,&
                               k_summed_overlap_matr,&
                               !inv_k_summed_overlap_matr, &
                               hybrid_func,loc_self_enrg,embed_gf,&
                               hartree_pot_LDA,&
                               number_reiterations,&
                               outter_loop_reiteration,& ! (KS_eigenvalue)
                               new_chemical_potential,& ! (KS_eigenvalue)
                               self_enrg_gw) ! (KS_eigenvalue)

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
   integer a, i_a 
   integer i_matrix_size
   integer j_matrix_size
   integer j_basis
   integer i_basis, i_basis_1
   integer i_freq 
   integer i_k_point
   integer   number_reiterations, outter_loop_reiteration
   character*2 iter
   character*2 iter1
   character*17 filename_RE
   character*17 filename_IM
 
! subroutine actual parameters
   real*8, dimension(n_basis,n_basis), intent(in) :: free_cluster_ham
   real*8, dimension(n_basis,n_basis), intent(in) :: loc_self_enrg 
   complex*16, dimension(n_basis,n_basis,nomega) :: hybrid_func
   real*8, dimension(n_basis,n_basis), intent(in) :: k_summed_overlap_matr
!   real*8, dimension(n_basis,n_basis), intent(in) :: inv_k_summed_overlap_matr
!   complex*16, dimension(n_basis,n_basis,nomega), intent(out) :: inv_embed_gf
   complex*16, dimension(n_basis,n_basis,nomega), intent(out) :: embed_gf
   !complex*16, dimension(n_basis,n_basis,nomega), intent(in) :: hybrid_func_temp
!   complex*16, dimension(:,:), allocatable :: inv_k_summed_overlap_matr


!   complex*16, dimension(:,:), allocatable:: inv_k_summed_overlap_matr_sqrt
!   complex*16, dimension(:,:), allocatable :: k_summed_overlap_matr_sqrt
   complex*16, dimension(:,:,:), allocatable :: inv_embed_gf 
   complex*16, dimension(n_basis,n_basis,nomega):: omega_term
   complex*16, dimension(n_basis,n_basis,nomega):: self_enrg_gw
   real*8, dimension(n_basis,n_basis) :: hartree_pot_LDA
   real*8 :: new_chemical_potential 


   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 


   logical  :: output

     if (.not. allocated (inv_embed_gf)) then
         allocate (inv_embed_gf(n_basis,n_basis,nomega))
     endif
 
!     if (.not. allocated (inv_k_summed_overlap_matr_sqrt)) then
!         allocate (inv_k_summed_overlap_matr_sqrt(n_basis,n_basis))
!     endif
!     if (.not. allocated (inv_k_summed_overlap_matr)) then
!         allocate (inv_k_summed_overlap_matr(n_basis,n_basis))
!     endif
!     if (.not. allocated (k_summed_overlap_matr_sqrt)) then
!         allocate (k_summed_overlap_matr_sqrt(n_basis,n_basis))
!     endif




       if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif



!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!         full_ovlp_matrix_sqrt(:,:,:) = full_ovlp_matrix(:,:,:)
!         inv_k_summed_overlap_matr_sqrt(:,:) = k_summed_overlap_matr(:,:) 

!!write(use_unit,*)'safe_minimum', safe_minimum
!do i_k_point = 1, n_k_points
!     call power_genmat_lapack_complex(n_basis,inv_k_summed_overlap_matr_sqrt(:,:), -0.5, &
!                             safe_minimum, 1.d-5 , '') 

!ovlp_hybrid_func_ovlp(:,:,:)=0.d0


!            k_summed_overlap_matr_sqrt(:,:), &
!            n_basis,0.d0, &
!            ovlp_hybrid_func_ovlp(:,:,i_freq) , n_basis )

!           call zgemm ('N','N', n_basis,&
!           n_basis,n_basis,1.d0,&
!            hybrid_func(:,:,i_freq), n_basis, &
!            k_summed_overlap_matr(:,:), &
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

          !if(number_reiterations.eq.0) then
          !    loc_self_enrg = on_site_xc_matr
          !endif
!          if(number_reiterations.eq.0) then
!              hartree_pot_LDA = 0.d0
!          endif

 
        inv_embed_gf(:,:,:) = (0.d0,0.d0)
       do i_freq = 1, nomega, 1 
          do i_matrix_size = 1, n_basis, 1 
            do j_matrix_size = 1, n_basis, 1 

                
                  inv_embed_gf(i_matrix_size,j_matrix_size,i_freq) = &
                  inv_embed_gf(i_matrix_size,j_matrix_size,i_freq) + &
                  !inv_embed_gf(i_matrix_size,j_matrix_size,i_freq) = &
                  !inv_embed_gf(i_matrix_size,j_matrix_size,i_freq) + &
                  (((0.d0,1.d0)*omega(i_freq)+chemical_potential)* &
                  k_summed_overlap_matr(i_matrix_size,j_matrix_size)- &
                  !((ovlp_omega_ovlp(i_matrix_size,j_matrix_size,i_freq))- &
                  !((omega_term(i_matrix_size,j_matrix_size,i_freq))- &
!                  free_cluster_ham(i_matrix_size,j_matrix_size)- &
                  !hartree_pot_LDA(i_matrix_size,j_matrix_size)- &
                  loc_self_enrg(i_matrix_size,j_matrix_size)- &
                  !self_enrg_gw(i_matrix_size,j_matrix_size,i_freq)- &
                  hybrid_func(i_matrix_size,j_matrix_size, i_freq))


           enddo
          enddo
         enddo

!       do i_basis_1 = 1, n_basis, 1
!            write(use_unit,*)'H and EXX from embed GF no_ovlp',  &
!               hartree_pot_LDA(i_basis_1, i_basis_1), loc_self_enrg(i_basis_1, i_basis_1)
!       enddo



!      do i_basis_1 = 1, n_basis, 1
!             write(use_unit,*)'loc_self_enrg from embed_GF',  i_basis_1, &
!               loc_self_enrg(i_basis_1, i_basis_1)
!       enddo

!      do i_basis_1 = 1, n_basis, 1
!             write(use_unit,*)'hartree_pot_LDA from embed_GF',  i_basis_1, &
!               hartree_pot_LDA(i_basis_1, i_basis_1)
!       enddo



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
!            write(use_unit,*)'embed_gf from embed_GF',  i_basis_1, &
!             embed_gf(i_basis_1, i_basis_1,100)
!      enddo



!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


         if(myid.eq.0)then

     !output = .false.
     output = .true.
      if(output) then
           ! do a = 1, n_freq, 1
!        do i_a = 1, n_basis, 1
          if( i_a.lt.10 ) then
             write(iter,'(A,I1)') "0",outter_loop_reiteration
          else
             write(iter,'(I2)') outter_loop_reiteration
          endif
          if( number_reiterations.lt.10 ) then
             write(iter1,'(A,I1)') "0", number_reiterations
          else
             write(iter1,'(I2)') number_reiterations
          endif

          filename_RE = "embGF_RE_"//iter1//"_"//iter//".dat"
          filename_IM = "embGF_IM_"//iter1//"_"//iter//".dat"
        !    do a = 1, n_freq, 1
          open(77, file=filename_RE)
          open(75, file=filename_IM)
           do a = 1, nomega, 1
           !do j_a = 1, n_basis, 1
          !  do a = 1, n_freq, 1
              write(77,*) omega(a), &
                       real(embed_gf(1,1,a)), (1.d0)/omega(a)
              write(75,*) omega(a), &
                      aimag(embed_gf(1,1,a)) !aimag(1.d0/((0.d0,1.d0)*omega(a)))
            ! enddo
            enddo
          close(77)
          close(75)
!            enddo
       ! enddo
      endif
      endif

     if (allocated (inv_embed_gf)) then
         deallocate (inv_embed_gf)
     endif





       if (allocated(ipiv))then
          deallocate(ipiv)
       endif

       if ( allocated(work))then
          deallocate(work)
       endif




!endif !myid.eq.0
  end subroutine get_embed_gf_freq_no_ovlp 
