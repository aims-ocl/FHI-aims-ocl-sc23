  subroutine get_embed_gf_freq(free_cluster_ham,&
                               hybrid_func,&
                               loc_self_enrg,&
                               embed_gf,&
                               hartree_pot_LDA,&
                               number_reiterations,&
                               outter_loop_reiteration,&
                               self_enrg_gw,&
                               k_summed_overlap_matr,&
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
   real*8, dimension(n_basis,n_basis) :: k_summed_overlap_matr
   complex*16, dimension(n_basis,n_basis,nomega), intent(out) :: embed_gf

   complex*16, dimension(:,:), allocatable :: inv_embed_gf 
   complex*16, dimension(:,:), allocatable :: aux_matr 
   complex*16, dimension(n_basis,n_basis,nomega):: self_enrg_gw
   real*8, dimension(n_basis,n_basis) :: hartree_pot_LDA
   real*8  new_chemical_potential_cluster


   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 


   logical  :: output




!     if (.not. allocated (inv_embed_gf)) then
!         allocate (inv_embed_gf(n_basis,n_basis))
!     endif
 

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
 
!        inv_embed_gf(:,:,:) = (0.d0,0.d0)
       do i_freq = 1, nomega, 1 

     if (.not. allocated (inv_embed_gf)) then
         allocate (inv_embed_gf(n_basis,n_basis))
     endif

        inv_embed_gf(:,:) = (0.d0,0.d0)

          do i_matrix_size = 1, n_basis, 1 
            do j_matrix_size = 1, n_basis, 1 

                
                  inv_embed_gf(i_matrix_size,j_matrix_size) = &
                  inv_embed_gf(i_matrix_size,j_matrix_size) + &
                  (((0.d0,1.d0)*omega(i_freq) + &
                  new_chemical_potential_cluster)* &
                  k_summed_overlap_matr(i_matrix_size,j_matrix_size)-&!+ &
                  !free_cluster_ham(i_matrix_size,j_matrix_size)- &
                  !hartree_pot_LDA(i_matrix_size,j_matrix_size)- &
                  loc_self_enrg(i_matrix_size,j_matrix_size)- &
                  self_enrg_gw(i_matrix_size,j_matrix_size,i_freq)- &
                  hybrid_func(i_matrix_size,j_matrix_size, i_freq))


           enddo
          enddo
!         enddo


! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------------


       embed_gf(:,:,i_freq) = inv_embed_gf(:,:)

! at each frequency point evaluate the LU factorization, and check for errors
!     do i_freq = 1, nomega, 1

       if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif


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
        call zgetri( n_basis, embed_gf(:,:,i_freq), n_basis, &
                 ipiv, work, n_basis, info)

      if (info.ne.0) then
         if(myid.eq.0)then
            write(use_unit,*) " * Failure of matrix inversion at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
          endif
        endif

       if (allocated(ipiv))then
          deallocate(ipiv)
       endif

       if ( allocated(work))then
          deallocate(work)
       endif

     if (allocated (inv_embed_gf)) then
         deallocate (inv_embed_gf)
     endif

      enddo


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


         if(myid.eq.0)then

     output = .false.
     !output = .true.
      if(output) then
      if(number_reiterations.eq.0) then
           ! do a = 1, n_freq, 1
        do i_a = 1, n_basis, 1
          if( i_a.lt.10 ) then
             write(iter,'(A,I1)') "0",i_a
          else
             write(iter,'(I2)') i_a
          endif
!          if( number_reiterations.lt.10 ) then
          if( outter_loop_reiteration.lt.10 ) then
             write(iter1,'(A,I1)') "0", outter_loop_reiteration
          else
             write(iter1,'(I2)') outter_loop_reiteration
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
                       real(embed_gf(i_a,i_a,a)),&
                       real(self_enrg_gw(i_a,i_a,a)),a,omega(a),womega(a)
                       !real(self_enrg_gw(i_a,i_a,a))
              write(75,*) omega(a), &
                      aimag(embed_gf(i_a,i_a,a)),&
                      aimag(self_enrg_gw(i_a,i_a,a)) 
                      !aimag(self_enrg_gw(i_a,i_a,a)) !aimag(1.d0/((0.d0,1.d0)*omega(a)))
            ! enddo
            enddo
          close(77)
          close(75)
            !enddo
        enddo
      endif
      endif
      endif


!endif !myid.eq.0
  end subroutine get_embed_gf_freq 
