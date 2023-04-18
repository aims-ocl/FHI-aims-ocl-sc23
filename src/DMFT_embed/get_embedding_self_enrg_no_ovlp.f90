  subroutine get_embedding_self_enrg_no_ovlp(free_cluster_ham,&
                                     inv_on_site_gf, &
                                     on_site_xc_matr,&
                                     k_summed_overlap_matr,&
!                                     inv_k_summed_overlap_matr,&
                                     hybrid_func,&
                                     loc_self_enrg,&
                                     loc_self_enrg_GW,&
                                     hartree_pot_LDA,number_reiterations,&
                                     !omega_term) ! (KS_eigenvalue)
                                     new_chemical_potential) ! (KS_eigenvalue)

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
!   real*8, dimension(n_basis,n_basis) :: inv_k_summed_overlap_matr
   !complex*16, dimension(n_basis,n_basis, nomega) :: hybrid_func_temp
   complex*16, dimension(n_basis,n_basis, nomega) :: hybrid_func
   complex*16, dimension(n_basis,n_basis, nomega) :: loc_self_enrg_GW
   !complex*16, dimension(n_basis,n_basis, nomega) :: inv_hybrid_func
   real*8, dimension(n_basis,n_basis):: loc_self_enrg
!   real*8, dimension(n_basis,n_basis):: loc_self_enrg
   real*8, dimension(n_basis,n_basis) :: on_site_xc_matr
!   real*8, dimension(n_basis,n_basis) :: on_site_xc_matr
   complex*16, dimension(n_basis,n_basis,nomega) :: inv_on_site_gf



!   real*8, dimension(:,:), allocatable :: inv_k_summed_overlap_matr_sqrt
!   complex*16, dimension(:,:), allocatable :: inv_k_summed_overlap_matr
!   complex*16, dimension(:,:), allocatable :: k_summed_overlap_matr_sqrt
   !complex*16, dimension(:,:,:), allocatable :: omega_term
   complex*16, dimension(n_basis,n_basis,nomega) :: omega_term
   complex*16, dimension(:,:,:), allocatable :: inv_hybrid_func
   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 
   real*8, dimension(n_basis,n_basis):: hartree_pot_LDA
   real*8 :: new_chemical_potential

   integer number_reiterations



   logical  :: output




!----------------------------------------------------------------------------
!----------------------------------------------------------------------------


!       if (.not. allocated(inv_k_summed_overlap_matr_sqrt))then
!          allocate(inv_k_summed_overlap_matr_sqrt(n_basis,n_basis))
!       endif
!       if (.not. allocated(inv_k_summed_overlap_matr))then
!          allocate(inv_k_summed_overlap_matr(n_basis,n_basis))
!       endif
!       if (.not. allocated(k_summed_overlap_matr_sqrt))then
!          allocate(k_summed_overlap_matr_sqrt(n_basis,n_basis))
!       endif


       if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif


       if (.not. allocated(inv_hybrid_func))then
          allocate(inv_hybrid_func(n_basis,n_basis,nomega))
       endif



!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!do i_basis = 1, n_basis,1

!write(use_unit,*) 'loc_self_enrg from embedding SE', loc_self_enrg(i_basis,i_basis)

!enddo


! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------

          if(number_reiterations.eq.0) then
              loc_self_enrg = on_site_xc_matr
!              hartree_pot_LDA = 0.d0
               loc_self_enrg_GW = (0.d0,0.d0)
               new_chemical_potential = chemical_potential
          endif


!do i_basis = 1, n_basis,1

!write(use_unit,*) 'hartree_pot_LDA from embedding SE', hartree_pot_LDA(i_basis,i_basis)
!write(use_unit,*) 'inv_on_site_gf from embedding SE', loc_self_enrg(i_basis,i_basis)

!enddo



!write(use_unit,*), hartree_pot_LDA(:,:), loc_self_enrg(:,:) 

!         hybrid_func_temp(:,:,:) = (0.d0,0.d0)
         hybrid_func(:,:,:) = (0.d0,0.d0)
        do i_freq = 1, nomega, 1
          do i_matrix_size = 1, n_basis, 1
            do j_matrix_size = 1, n_basis, 1


                  hybrid_func(i_matrix_size,j_matrix_size, i_freq) = &
                  hybrid_func(i_matrix_size,j_matrix_size, i_freq) + &
                  ((((0.d0,1.d0)*omega(i_freq))+chemical_potential)* &
                  k_summed_overlap_matr(i_matrix_size,j_matrix_size)- &
!                  free_cluster_ham(i_matrix_size,j_matrix_size)- &
                  inv_on_site_gf(i_matrix_size,j_matrix_size,i_freq) - &
                  loc_self_enrg(i_matrix_size,j_matrix_size))!-loc_self_enrg_GW(i_matrix_size,j_matrix_size,i_freq))!&


           enddo
          enddo
         enddo
!write(use_unit,*) "hartree_pot_LDA from embedding SE", hartree_pot_LDA(:,:)
!write(use_unit,*) "loc_self_enrg from embedding SE", loc_self_enrg(:,:)

! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------

! write(use_unit,*) 'loc_self_enrg from embedding SE', loc_self_enrg(:,:)
! write(use_unit,*) 'loc_self_enrg from embedding SE', hartree_pot_LDA(:,:)

!ovlp_hybrid_func(:,:,:)=0.d0


!      do i_basis_1 = 1, n_basis, 1
!         do i_basis_2 = 1, i_basis_1, 1
             !write(use_unit,*)'on_site_xc_pot from embeddinf self enrg',  i_basis_1, &
             !  on_site_xc_matr(i_basis_1, i_basis_1)
!             write(use_unit,*)'hybrid_func from embedding SE',  i_basis_1, &
!               hybrid_func(i_basis_1, i_basis_1,10)
!         enddo
!       enddo


!       do i_basis_1 = 1, n_basis, 1
!            write(use_unit,*)'H and EXX from embedding_no_ovlp',  &
!               hartree_pot_LDA(i_basis_1, i_basis_1), loc_self_enrg(i_basis_1, i_basis_1)
!       enddo

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



!  do i_a= 1, n_basis
!    write(use_unit,*)'cluster_hamiltoniani from embeddinf self-enrg', i_a, free_cluster_ham(i_a,i_a)
!  enddo
 ! do i_a= 1, n_basis
 !   write(use_unit,*)'k_summed_ovlp from embedding self-enrg', i_a, k_summed_overlap_matr(i_a,i_a), k_summed_overlap_matr(2,i_a)
 ! enddo


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
        !    do a = 1, n_freq, 1
          open(77, file=filename_RE)
          open(75, file=filename_IM)
           do i_freq = 1, nomega, 1
           !do j_a = 1, n_basis, 1
          !  do a = 1, n_freq, 1
              write(77,*) omega(i_freq), &
                       real(inv_hybrid_func(i_a,i_a,i_freq)), real(hybrid_func(i_a,i_a,i_freq))
              write(75,*) omega(i_freq), &
                      aimag(inv_hybrid_func(i_a,i_a,i_freq)), aimag(hybrid_func(i_a,i_a,i_freq))
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
  end subroutine get_embedding_self_enrg_no_ovlp 
