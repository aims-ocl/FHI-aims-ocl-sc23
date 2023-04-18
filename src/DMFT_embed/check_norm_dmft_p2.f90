  subroutine check_norm_dmft_p2(inv_full_ovlp_matrix_sqrt, &
                                full_ovlp_matrix_sqrt, &
                                full_ovlp_matrix, &
                                free_cluster_ham, &
                                loc_self_enrg,&
                                hartree_pot_LDA, &
                                elec_N, &
                                diff_in_elec_N, &
                                chemical_pot_i, &
                                i_counter, &
                                embed_part_number_new,&
                                on_site_xc_matr_no_ovlp,&
                                full_hamiltonian, &
                                on_site_PBE_xc_matr_no_ovlp, &
                                loc_self_enrg_GW)


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
   integer l,a,b,n, i_a, j_a,i_k 
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

   complex*16, dimension(n_basis,n_basis, n_k_points), intent(in) :: &
   inv_full_ovlp_matrix_sqrt 
   complex*16, dimension(n_basis,n_basis, n_k_points), intent(in) :: &
   full_ovlp_matrix_sqrt 
   complex*16, dimension(n_basis,n_basis, n_k_points), intent(in) :: &
   full_ovlp_matrix 
   complex*16, dimension(n_basis,n_basis,n_spin, n_k_points), intent(in) :: &
   full_hamiltonian 
   complex*16, dimension(:,:,:), allocatable:: embed_gf
   complex*16, dimension(:,:,:), allocatable:: inv_embed_gf
   real*8, dimension(:,:), allocatable:: embed_dens_matr_temp
   real*8, dimension(:,:), allocatable:: embed_dens_matr
   real*8, dimension(:,:), allocatable:: loc_self_enrg_temp

   real*8, dimension(n_basis,n_basis), intent(in) :: free_cluster_ham
   real*8, dimension(n_basis,n_basis), intent(in) :: loc_self_enrg
   complex*16, dimension(:,:,:,:), allocatable:: omega_ovlp
   real*8, dimension(n_basis,n_basis) :: hartree_pot_LDA
   real*8, dimension(n_basis,n_basis) :: on_site_xc_matr_no_ovlp
   real*8, dimension(n_basis,n_basis) :: on_site_PBE_xc_matr_no_ovlp
   integer, intent(inout) :: i_counter
   complex*16, dimension(:,:,:), allocatable :: gf_non_loc
   complex*16, dimension(:,:,:), allocatable :: gf_non_loc_temp
   complex*16, dimension(:,:,:), allocatable :: ovlp_gf_non_loc
   complex*16 , dimension (:,:),allocatable :: loc_self_enrg_no_orth
   complex*16 , dimension (:,:),allocatable :: hartree_pot_LDA_no_orth
   complex*16 , dimension (:,:), allocatable :: aux_no_orth
   complex*16 , dimension (:,:), allocatable :: aux_matr_1
   complex*16 , dimension (:,:,:), allocatable :: KS_eigenvector_tmp
   complex*16 , dimension (n_basis,n_basis,nomega) :: loc_self_enrg_GW

   real*8, intent(out) :: diff_in_elec_N 
   real*8, intent(out) :: embed_part_number_new 
   real*8, intent(in) :: elec_N
   real*8 :: embed_part_number
   real*8 ::chemical_pot_i 
   !real*8, dimension(n_basis,n_basis,-ntau:ntau,n_spin) :: embed_gf_time
   real*8, dimension(:,:,:), allocatable :: embed_gf_time

      logical  :: output

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------


      if (.not. allocated(embed_dens_matr))then
          allocate(embed_dens_matr(n_basis,n_basis))
       endif
      if (.not. allocated(embed_dens_matr_temp))then
          allocate(embed_dens_matr_temp(n_basis,n_basis))
       endif
       
      if (.not. allocated(ipiv))then
          !allocate(ipiv(n_basis))
          allocate(ipiv(n_states))
       endif

       if (.not. allocated(work))then
          !allocate(work (n_basis))
          allocate(work (n_states))
       endif

      if (.not. allocated(embed_gf))then
          allocate(embed_gf(n_basis,n_basis,nomega))
       endif
      if (.not. allocated(embed_gf_time))then
          allocate(embed_gf_time(n_basis,n_basis,-ntau:ntau))
       endif
!      if (.not. allocated(inv_embed_gf))then
!          allocate(inv_embed_gf(n_basis,n_basis,nomega))
!       endif



       if (.not. allocated(loc_self_enrg_temp))then
          allocate(loc_self_enrg_temp(n_basis,n_basis))
       endif

    loc_self_enrg_temp(:,:)= loc_self_enrg(:,:) - on_site_PBE_xc_matr_no_ovlp(:,:)

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
  i_counter = i_counter + 1

  embed_gf(:,:,:)=(0.d0,0.d0)
    i_k = 0

       do i_k_point=1,n_k_points,1


      if(.not. allocated(KS_eigenvector_tmp))then
        allocate(KS_eigenvector_tmp(n_basis,n_states,n_spin))
      endif

           if (myid .eq. MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then

            i_k = i_k + 1

             if(real_eigenvectors)then

               KS_eigenvector_tmp(:,:,:) = DCMPLX(KS_eigenvector(:,:,:,i_k))

              else

               KS_eigenvector_tmp(:,:,:) = KS_eigenvector_complex(:,:,:,i_k)

             endif

           else
             ! zero temp. KS eigenvector on all other threads, prior to allreduce
             KS_eigenvector_tmp(:,:,:) = (0.d0,0.d0)
          end if
           call sync_eigenvector_complex(KS_eigenvector_tmp)



      if (.not. allocated(inv_embed_gf))then
          allocate(inv_embed_gf(n_basis,n_basis,nomega))
       endif

  inv_embed_gf(:,:,:) = (0.d0,0.d0)


       do i_freq = 1, nomega, 1
          do i_basis = 1, n_basis, 1
            do j_basis = 1, n_basis, 1


                  inv_embed_gf(i_basis,j_basis,i_freq) = &
                  inv_embed_gf(i_basis,j_basis,i_freq) + &
                  (((0.d0,1.d0)*omega(i_freq)+DCMPLX(chemical_pot_i))* &
                  full_ovlp_matrix(i_basis,j_basis,i_k_point)- &
                  full_hamiltonian(i_basis,j_basis,1,i_k_point)+ &
                  DCMPLX(on_site_xc_matr_no_ovlp(i_basis,j_basis)) - &
                  DCMPLX(loc_self_enrg(i_basis,j_basis)) - &
                  loc_self_enrg_GW(i_basis,j_basis,i_freq))!- &



            enddo
           enddo
       enddo


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

       if (.not. allocated(gf_non_loc_temp))then
          allocate(gf_non_loc_temp(n_states,n_states,nomega))
       endif

gf_non_loc_temp(:,:,:) = (0.d0,0.d0)
do i_freq = 1, nomega

                call transform_to_KS_basis_complex(&
                                      inv_embed_gf(:,:,i_freq),&
                                      i_k_point,&
                                      gf_non_loc_temp(:,:,i_freq),&
                                      KS_eigenvector_tmp)

enddo

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


! at each frequency point evaluate the LU factorization, and check for errors
       do i_freq = 1, nomega, 1
!       do i_k_point = 1, n_k_points, 1
         call zgetrf( n_states, n_states, gf_non_loc_temp(:,:,i_freq), &
                    n_states , ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then
            write(use_unit,*) " * Failure of LU decomposition at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
  !          stop
          endif
        endif
! this inverts the matrix in the LU factorized form
        call zgetri(n_states, gf_non_loc_temp(:,:,i_freq), n_states, &
                 ipiv, work, n_states, info)

      if (info.ne.0) then
         if(myid.eq.0)then
            write(use_unit,*) " * Failure of matrix inversion at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
          endif
        endif
       enddo

       if ( allocated(inv_embed_gf))then
          deallocate(inv_embed_gf)
       endif

       if (.not. allocated(ovlp_gf_non_loc))then
          allocate(ovlp_gf_non_loc(n_basis,n_basis,nomega))
       endif
       if (.not. allocated(gf_non_loc))then
          allocate(gf_non_loc(n_basis,n_basis,nomega))
       endif


ovlp_gf_non_loc(:,:,:)=(0.d0,0.d0)
gf_non_loc(:,:,:)=(0.d0,0.d0)
do i_freq = 1, nomega

!
                call transform_to_NAO_basis_complex(&
                       gf_non_loc_temp(:,:,i_freq),i_k_point,&
                       gf_non_loc(:,:,i_freq),KS_eigenvector_tmp)




       embed_gf(:,:,i_freq) = embed_gf(:,:,i_freq) + gf_non_loc(:,:,i_freq) 



    enddo


       if (allocated(ovlp_gf_non_loc))then
          deallocate(ovlp_gf_non_loc)
       endif
       if (allocated(gf_non_loc_temp))then
          deallocate(gf_non_loc_temp)
       endif
       if (allocated(gf_non_loc))then
          deallocate(gf_non_loc)
       endif

      if(allocated(KS_eigenvector_tmp))then
         deallocate(KS_eigenvector_tmp)
      endif

  enddo


      


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
       embed_gf(:,:,:)=(1./(n_k_points))*embed_gf(:,:,:)
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


if(.false.) then
!if(.true.) then

         embed_dens_matr_temp(:,:) = 0.d0
        do i_freq = 1, nomega
            do j_matrix_size = 1, n_basis, 1
              do i_matrix_size = 1, n_basis, 1
                  embed_dens_matr_temp(i_matrix_size,j_matrix_size) = &
                  embed_dens_matr_temp(i_matrix_size,j_matrix_size) + &
                  womega(i_freq) * (&
                      real(embed_gf(i_matrix_size,j_matrix_size, i_freq))*&
                      cos((omega(i_freq))*(3.3523E-4)) &
                      -aimag(embed_gf(i_matrix_size,j_matrix_size, i_freq))*&
                      sin((omega(i_freq))*(3.3523E-4)) &
                 )
                 enddo
             enddo
          enddo
endif

!       do   i_spin = 1,n_spin,1

            call transform_G (embed_gf , n_basis, &
               n_basis, embed_gf_time(:,:,:))

!       enddo

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

 embed_dens_matr(:,:) = (1./(pi))*embed_dens_matr_temp(:,:) 

 embed_part_number =0.d0

 do i_a = 1, n_basis,1

    embed_part_number= embed_part_number + (embed_gf_time(i_a,i_a,0))

 enddo

if(myid.eq.0)then
write(use_unit,*) i_counter, chemical_pot_i, embed_part_number
endif


 diff_in_elec_N = embed_part_number - elec_N

embed_part_number_new = embed_part_number
!stop
!write(use_unit,*) "diff_electrons",i_counter, diff_in_elec_N
     !output = .true.
     output = .false.
      if(output) then

          if( i_counter.lt.10 ) then
             write(iter,'(A,I1)') "0",i_counter
          else
             write(iter,'(I2)') i_counter
          endif

          filename_RE = "ChemPot_RE_"//iter//".dat"
          filename_IM = "ChemPot_IM_"//iter//".dat"
        !    do a = 1, n_freq, 1
          open(77, file=filename_RE)
          open(75, file=filename_IM)
           do a = 1, nomega, 1
              write(77,*) omega(a), &
                       real(embed_gf(1,1,a))
              write(75,*) omega(a), &
                      aimag(embed_gf(1,1,a)) !aimag(1.d0/((0.d0,1.d0)*omega(a)))
            enddo
          close(77)
          close(75)
!            enddo
       ! enddo
      endif

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


      if ( allocated(embed_dens_matr))then
          deallocate(embed_dens_matr)
       endif
      if ( allocated(embed_dens_matr_temp))then
          deallocate(embed_dens_matr_temp)
       endif

      if ( allocated(ipiv))then
          deallocate(ipiv)
       endif

       if ( allocated(work))then
          deallocate(work )
       endif

      if ( allocated(embed_gf))then
          deallocate(embed_gf)
       endif
      !if ( allocated(inv_embed_gf))then
      !    deallocate(inv_embed_gf)
      ! endif

  end subroutine check_norm_dmft_p2 
