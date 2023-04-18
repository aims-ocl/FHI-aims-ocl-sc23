  subroutine check_norm_dmft_p2_cluster(inv_full_ovlp_matrix_sqrt,&
                                        full_ovlp_matrix_sqrt,&
                                        full_ovlp_matrix,&
                                        free_cluster_ham,&
                                        hybrid_func,& 
                                        loc_self_enrg,&
                                        hartree_pot_LDA,&
                                        elec_N,&
                                        diff_in_elec_N,&
                                        chemical_pot_i,&
                                        i_counter,&
                                        embed_part_number_new,&
                                        k_summed_overlap_matr,&
                                        self_energy_freq)


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
  use poles_fit 
  use gt
  use timing

   implicit none

! local parameters
   integer l,a,b,n, i_a, j_a 
   integer hamiltonian_size
   integer ovlp_matrix_size
   integer i_hamilton_size 
   integer inner_loop_reiterations, number_reiterations
   integer i_matrix_size
   integer j_matrix_size
   integer i_spin
   integer i_state
   integer j_basis
   integer k_basis
   integer i_basis
   integer i_basis_1
   integer i_freq
   integer i_k_point
   character*2 iter
   character*17 filename_RE
   character*17 filename_OUT
   character*17 filename_IM

! subroutine actual parameters

   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 

   complex*16, dimension(n_basis,n_basis, n_k_points), intent(in) :: inv_full_ovlp_matrix_sqrt 
   complex*16, dimension(n_basis,n_basis, n_k_points), intent(in) :: full_ovlp_matrix_sqrt 
   complex*16, dimension(n_basis,n_basis, n_k_points), intent(in) :: full_ovlp_matrix 
   complex*16, dimension(:,:,:), allocatable:: embed_gf
   complex*16, dimension(:,:,:), allocatable:: inv_embed_gf
   real*8, dimension(:,:,:), allocatable:: embed_gf_time

   real*8, dimension(n_basis,n_basis), intent(in) :: free_cluster_ham
   real*8, dimension(n_basis,n_basis), intent(in) :: k_summed_overlap_matr
   real*8, dimension(n_basis,n_basis), intent(in) :: loc_self_enrg
   complex*16, dimension(n_basis,n_basis,nomega) :: hybrid_func
   complex*16, dimension(n_basis,n_basis,nomega) :: self_energy_freq
   complex*16, dimension(:,:,:,:), allocatable:: omega_ovlp
   real*8, dimension(n_basis,n_basis) :: hartree_pot_LDA
   integer, intent(inout) :: i_counter


   real*8, intent(out) :: diff_in_elec_N 
   real*8, intent(out) :: embed_part_number_new 
   real*8, intent(in) :: elec_N
   real*8 :: embed_part_number
   real*8 ::chemical_pot_i 

      logical  :: output

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------


       
      if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif

      if (.not. allocated(embed_gf))then
          allocate(embed_gf(n_basis,n_basis,nomega))
       endif
      if (.not. allocated(embed_gf_time))then
          allocate(embed_gf_time(n_basis,n_basis,-ntau:ntau))
       endif
      if (.not. allocated(inv_embed_gf))then
          allocate(inv_embed_gf(n_basis,n_basis,nomega))
       endif




!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  i_counter = i_counter + 1




!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


        inv_embed_gf(:,:,:) = (0.d0,0.d0)

         do i_freq = 1, nomega, 1
          do i_matrix_size = 1, n_basis, 1
           do j_matrix_size = 1, n_basis, 1


                  inv_embed_gf(i_matrix_size,j_matrix_size,i_freq) = &
                  inv_embed_gf(i_matrix_size,j_matrix_size,i_freq) + &
                  (((0.d0,1.d0)*omega(i_freq) + chemical_pot_i)* &
                  k_summed_overlap_matr(i_matrix_size,j_matrix_size)-&!+&
                  !((omega_term(i_matrix_size,j_matrix_size,i_freq))+ &
!                  free_cluster_ham(i_matrix_size,j_matrix_size)- &
                  !hartree_pot_LDA(i_matrix_size,j_matrix_size)- &
                  loc_self_enrg(i_matrix_size,j_matrix_size)- &
                  self_energy_freq(i_matrix_size,j_matrix_size, i_freq)- &
                  hybrid_func(i_matrix_size,j_matrix_size, i_freq))


           enddo
          enddo
         enddo



!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

       embed_gf(:,:,:) = inv_embed_gf(:,:,:)

! at each frequency point evaluate the LU factorization, and check for errors
     do i_freq = 1, nomega, 1

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
      enddo

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


             call transform_G (embed_gf , n_basis, &
               n_basis, embed_gf_time(:,:,:))

 embed_part_number =0.d0
 do i_a = 1, n_basis,1

    embed_part_number= embed_part_number + (embed_gf_time(i_a,i_a,0))

 enddo



!       if(myid.eq.0) then
!

              if(myid.eq.0) write(use_unit,*) 'particle number',  embed_part_number, chemical_pot_i

 diff_in_elec_N = embed_part_number - elec_N
!              write(use_unit,*) 'diff in particle number',  diff_in_elec_N
!              write(use_unit,*) 'max_zeroin', max_zeroin 

embed_part_number_new = embed_part_number
!stop
!write(use_unit,*) "diff_electrons",i_counter, diff_in_elec_N

!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------


      if ( allocated(ipiv))then
          deallocate(ipiv)
       endif

       if ( allocated(work))then
          deallocate(work )
       endif

      if ( allocated(embed_gf))then
          deallocate(embed_gf)
       endif
      if ( allocated(inv_embed_gf))then
          deallocate(inv_embed_gf)
       endif

  end subroutine check_norm_dmft_p2_cluster 
