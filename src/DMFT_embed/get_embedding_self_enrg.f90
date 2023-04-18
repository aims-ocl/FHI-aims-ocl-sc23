  subroutine get_embedding_self_enrg(free_cluster_ham,&
                                     inv_on_site_gf, &
                                     on_site_xc_matr,&
                                     k_summed_overlap_matr,&
                                     hybrid_func,&
                                     loc_self_enrg,&
                                     loc_self_enrg_GW,&
                                     hartree_pot_LDA,&
                                     number_reiterations,&
                                     new_chemical_potential_cluster) 

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
   integer  i_a 
   integer i_matrix_size
   integer j_matrix_size
   integer i_basis, i_basis_1
   integer j_basis
   integer i_freq 
   integer i_k_point
   character*2 iter
   character*17 filename_RE
   character*17 filename_IM
 
! subroutine actual parameters
   real*8, dimension(n_basis,n_basis) :: free_cluster_ham
   real*8, dimension(n_basis,n_basis) :: k_summed_overlap_matr
   complex*16, dimension(n_basis,n_basis, nomega) :: hybrid_func
   complex*16, dimension(n_basis,n_basis, nomega) :: loc_self_enrg_GW
   real*8, dimension(n_basis,n_basis):: loc_self_enrg
   real*8, dimension(n_basis,n_basis) :: on_site_xc_matr
   complex*16, dimension(n_basis,n_basis,nomega) :: inv_on_site_gf

   complex*16, dimension(:,:,:), allocatable :: inv_hybrid_func
   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 
   real*8, dimension(n_basis,n_basis):: hartree_pot_LDA
   real*8 :: new_chemical_potential_cluster

   integer number_reiterations

   logical  :: output




!----------------------------------------------------------------------------
!----------------------------------------------------------------------------


       if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif

       if (.not. allocated(inv_hybrid_func))then
          allocate(inv_hybrid_func(n_basis,n_basis,nomega))
       endif



!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------

          if(number_reiterations.eq.0) then
              loc_self_enrg(:,:) = on_site_xc_matr(:,:)
              hartree_pot_LDA(:,:) = 0.d0
              loc_self_enrg_GW(:,:,:) = (0.d0,0.d0)
          endif






         hybrid_func(:,:,:) = (0.d0,0.d0)
        do i_freq = 1, nomega, 1
          do i_matrix_size = 1, n_basis, 1
            do j_matrix_size = 1, n_basis, 1


                  hybrid_func(i_matrix_size,j_matrix_size, i_freq) = &
                  hybrid_func(i_matrix_size,j_matrix_size, i_freq) + &
                  (((0.d0,1.d0)*omega(i_freq)+new_chemical_potential_cluster)* &
                  k_summed_overlap_matr(i_matrix_size,j_matrix_size)-&!+ &
                  !free_cluster_ham(i_matrix_size,j_matrix_size)- &
                  !hartree_pot_LDA(i_matrix_size,j_matrix_size)- &
                  inv_on_site_gf(i_matrix_size,j_matrix_size,i_freq) - &
                  loc_self_enrg(i_matrix_size,j_matrix_size)-&
                  loc_self_enrg_GW(i_matrix_size,j_matrix_size,i_freq))!&


           enddo
          enddo
         enddo

! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------



      inv_hybrid_func(:,:,:) = hybrid_func(:,:,:)

! at each frequency point evaluate the LU factorization, and check for errors
     do i_freq = 1, nomega, 1

         call zgetrf( n_basis, n_basis, inv_hybrid_func(:,:,i_freq), &
                    n_basis , ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then
            write(use_unit,*) " * Failure of LU decomposition at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
           ! stop
          endif
        endif
! this inverts the matrix in the LU factorized form
        call zgetri(n_basis, inv_hybrid_func(:,:,i_freq), n_basis, &
                 ipiv, work, n_basis, info)

      if (info.ne.0) then
         if(myid.eq.0)then
            write(use_unit,*) " * Failure of matrix inversion at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
          endif
        endif


      enddo


         if(myid.eq.0)then

     output = .false.
      if(output) then
           ! do a = 1, n_freq, 1
        do i_a = 1, n_basis, 1
          if( i_a.lt.10 ) then
             write(iter,'(A,I1)') "0",i_a
          else
             write(iter,'(I2)') i_a
          endif
          filename_RE = "hybr_RE_"//iter//".dat"
          filename_IM = "hybr_IM_"//iter//".dat"
          open(77, file=filename_RE)
          open(75, file=filename_IM)
           do i_freq = 1, nomega, 1
              write(77,*) omega(i_freq), &
                       real(inv_hybrid_func(i_a,i_a,i_freq)),&
                       real(hybrid_func(i_a,i_a,i_freq))
              write(75,*) omega(i_freq), &
                      aimag(inv_hybrid_func(i_a,i_a,i_freq)),&
                      aimag(hybrid_func(i_a,i_a,i_freq))
            ! enddo
            enddo
          close(77)
          close(75)
            enddo
       ! enddo
      endif
      endif


       if (allocated(ipiv))then
          deallocate(ipiv)
       endif

       if ( allocated(work))then
          deallocate(work)
       endif


       if ( allocated(inv_hybrid_func))then
          deallocate(inv_hybrid_func)
       endif


!endif !myid.eq.0
  end subroutine get_embedding_self_enrg 
