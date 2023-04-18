  subroutine get_NEW_on_site_gf_freq(full_hamiltonian,&
                                     full_ovlp_matrix,&
                                     on_site_gf,&
                                     on_site_xc_matr,&
                                     loc_self_enrg, &
                                     hartree_pot_LDA,& 
                                     number_reiterations,&
                                     scdmft_converged,&
                                     full_ovlp_matrix_sqrt, &
                                     loc_self_enrg_GW,&
                                     on_site_PBE_xc_matr,&
                                     dens_matr,&
                                     new_chemical_potential,&
                                     new_chemical_potential_cluster,&
                                     embed_part_number_LDA)!,&





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
   integer l,a,b,n, i_a, j_a, i_k,i
   integer hamiltonian_size
   integer ovlp_matrix_size
   integer i_basis, i_basis_1, i_basis_2, i_matrix_size, j_matrix_size
   integer j_basis
   integer i_k_point_local 
   integer k_basis
   integer i_count
   integer id_root 
   integer i_state
   integer i_freq
   integer i_k_point
   integer ham_size
   character*2 iter
   character*2 iter_1
   character*25 filename, filename_LDA
   character*20 filename_RE
   character*20 filename_IM
 
! subroutine actual parameters
   complex*16, dimension(n_basis,n_basis,n_k_points), intent(in) :: &
    full_ovlp_matrix 
   complex*16, dimension(n_basis,n_basis, n_spin,n_k_points), intent(in) &
   :: full_hamiltonian
   real*8, dimension(n_basis,n_basis),intent(in) :: loc_self_enrg 
   real*8, dimension(n_basis,n_basis),intent(in) :: on_site_xc_matr
   real*8  embed_part_number_LDA
   complex*16, dimension(n_basis,n_basis,nomega), intent(out) :: on_site_gf
   complex*16, dimension(:,:,:), allocatable :: gf_non_loc
   real*8, dimension(:,:,:), allocatable :: embed_gf_time
   complex*16, dimension(:,:), allocatable :: aux_matr

   complex*16, dimension(:,:,:), allocatable :: inv_gf_non_loc 
   complex*16, dimension(:,:,:), allocatable :: gf_non_loc_temp 
   complex*16, dimension(:,:), allocatable :: hamiltonian_GW 
   complex*16, dimension(:,:), allocatable :: hamiltonian_GW_temp 
   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 


   logical  :: output
   logical  :: scdmft_converged 
   logical  :: output_spectral
   integer number_reiterations
   integer max_reiteration 



      integer aux_nomega
      real*8  aux_omegamax
      real*8 , dimension (:), allocatable :: aux_omega
      real*8 , dimension (:), allocatable :: aux_womega
      real*8 , dimension (:), allocatable :: spectrum 
      real*8 , dimension (:), allocatable :: spectrum_k
      complex*16, dimension(:,:,:), allocatable:: green_fn_par
      real*8 , dimension (:), allocatable :: spectrum_LDA 
      real*8 , dimension (:), allocatable :: spectrum_k_LDA
      complex*16, dimension(:,:,:), allocatable:: green_fn_par_LDA
      real*8, dimension(n_basis,n_basis) :: hartree_pot_LDA
      complex*16 , dimension (n_basis,n_basis,nomega) :: loc_self_enrg_GW
!      complex*16, dimension(n_basis,n_basis,1,n_k_points) ::full_k_xc_matr


      complex*16, dimension(n_basis,n_basis,n_k_points) :: &
      full_ovlp_matrix_sqrt
      complex*16, dimension(:,:), allocatable:: diag_ham
      complex*16, dimension(:,:), allocatable:: diag_ham_k
      complex*16, dimension(:,:), allocatable:: diag_ham_k_LDA
      real*8 , dimension (:,:), allocatable :: hartree_pot_LDA_spec
      complex*16, dimension(:,:), allocatable :: ovlp_full_hamiltonian
      real*8, dimension(n_basis,n_basis) :: on_site_PBE_xc_matr
      complex*16, dimension(:,:), allocatable:: aux_matr_1
      complex*16, dimension(:,:,:), allocatable:: KS_eigenvector_tmp
      complex*16, dimension(:,:,:), allocatable:: DMFT_hamiltonian_k
      real*8 , dimension(n_basis,n_basis):: dens_matr
      real*8, dimension(:), allocatable :: DMFT_dos
      real*8, dimension(:), allocatable :: DMFT_dos_LDA
      real*8, dimension(:), allocatable :: DMFT_dos_k
      real*8, dimension(:), allocatable :: DMFT_dos_k_LDA


            real*8 new_chemical_potential
            real*8 new_chemical_potential_cluster


!------------------------------------------------------   


 

       if (.not. allocated(ipiv))then
          !allocate(ipiv(n_basis))
          allocate(ipiv(n_states))
       endif

       if (.not. allocated(work))then
          !allocate(work (n_basis))
          allocate(work (n_states))
       endif
       if (.not. allocated(hamiltonian_GW_temp))then
          allocate(hamiltonian_GW_temp(n_basis,n_basis))
       endif
       if (.not. allocated(hamiltonian_GW))then
          !allocate(hamiltonian_GW(n_states,n_states))
          allocate(hamiltonian_GW(n_basis,n_basis))
       endif
       if (.not. allocated(aux_matr))then
          allocate(aux_matr(n_basis,n_basis))
       endif
       if (.not. allocated(DMFT_hamiltonian_k))then
          !allocate(DMFT_hamiltonian_k(n_states,n_states,n_k_points))
          allocate(DMFT_hamiltonian_k(n_basis,n_basis,n_k_points))
       endif

    DMFT_hamiltonian_k(:,:,:) = (0.d0,0.d0)
    hamiltonian_GW(:,:) = (0.d0,0.d0)
    on_site_gf(:,:,:)=0.d0
    dens_matr (:,:) = 0.d0

    i_k =0

 do i_k_point = 1, n_k_points


      if(.not. allocated(KS_eigenvector_tmp))then
        allocate(KS_eigenvector_tmp(n_basis,n_states,n_spin))
      endif

           if (myid .eq. MOD(i_k_point, n_tasks) &
           .and. myid <= n_k_points ) then

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



!         do i_k_point = 1, n_k_points, 1 

       if (.not. allocated(inv_gf_non_loc))then
          allocate(inv_gf_non_loc(n_basis,n_basis,nomega))
       endif


        inv_gf_non_loc(:,:,:) = (0.d0,0.d0)


!----------------------------------------------------------------------------

       do i_freq = 1, nomega, 1 
          do j_basis = 1, n_basis, 1 
            do i_basis = 1, n_basis, 1 

                
                  inv_gf_non_loc(i_basis,j_basis,i_freq) = &
                  inv_gf_non_loc(i_basis,j_basis,i_freq) + &
                  (((0.d0,1.d0)*omega(i_freq)+new_chemical_potential)*&
                  full_ovlp_matrix(i_basis,j_basis,i_k_point)- &
                  full_hamiltonian(i_basis,j_basis,1,i_k_point)+ &
                  !DCMPLX(hartree_pot_LDA(i_basis,j_basis))+ &
                  DCMPLX(on_site_xc_matr(i_basis,j_basis))- &
                  DCMPLX(loc_self_enrg(i_basis,j_basis))- &
                  loc_self_enrg_GW(i_basis,j_basis,i_freq))



            enddo
           enddo
       enddo

!if(.false.)then
if(.true.)then

       if (.not. allocated(gf_non_loc_temp))then
          allocate(gf_non_loc_temp(n_states,n_states,nomega))
       endif

gf_non_loc_temp(:,:,:) = (0.d0,0.d0)
do i_freq = 1, nomega

                call transform_to_KS_basis_complex(&
                                      inv_gf_non_loc(:,:,i_freq),i_k_point,&
                                       gf_non_loc_temp(:,:,i_freq),&
                                       KS_eigenvector_tmp)

enddo

endif

!       if (.not. allocated(gf_non_loc_temp))then
!          allocate(gf_non_loc_temp(n_basis,n_basis,nomega))
!       endif
!          gf_non_loc_temp(:,:,:) = inv_gf_non_loc(:,:,:)
! at each frequency point evaluate the LU factorization, and check for errors
       do i_freq = 1, nomega, 1
!       write(use_unit,*) 'i_k_point= ', i_k_point 
         call zgetrf( n_states, n_states, gf_non_loc_temp(:,:,i_freq), &
                    n_states , ipiv, info )
!         call zgetrf( n_basis, n_basis, gf_non_loc_temp(:,:,i_freq), &
!                    n_basis , ipiv, info )

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
!        call zgetri(n_basis, gf_non_loc_temp(:,:,i_freq), n_basis, &
!                 ipiv, work, n_basis, info)

      if (info.ne.0) then
         if(myid.eq.0)then
            write(use_unit,*) " * Failure of matrix inversion at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
          endif
        endif
  
!       enddo
       enddo


       if ( allocated(inv_gf_non_loc))then
          deallocate(inv_gf_non_loc)
       endif
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

       if (.not. allocated(gf_non_loc))then
          allocate(gf_non_loc(n_basis,n_basis,nomega))
       endif


gf_non_loc(:,:,:)=(0.d0,0.d0)
!gf_non_loc_no_orth(:,:,i_k_point,:)=(0.d0,0.d0)
   do i_freq = 1, nomega


                call transform_to_NAO_basis_complex(&
                                       gf_non_loc_temp(:,:,i_freq),i_k_point,&
                                       gf_non_loc(:,:,i_freq),&
                                       KS_eigenvector_tmp)

on_site_gf(:,:,i_freq) = on_site_gf(:,:,i_freq) + gf_non_loc(:,:,i_freq)
!on_site_gf(:,:,i_freq) = on_site_gf(:,:,i_freq) + gf_non_loc_temp(:,:,i_freq)

  enddo

      if (allocated(gf_non_loc))then
          deallocate(gf_non_loc)
       endif

if(scdmft_converged)then

             if (.not. allocated(aux_matr_1))then
             allocate(aux_matr_1(n_basis,n_basis))
             endif

                 call multiply_ovlp_matr(DCMPLX(loc_self_enrg(:,:)),&
                                         full_ovlp_matrix_sqrt(:,:,i_k_point),&
                                         aux_matr_1(:,:))



        hamiltonian_GW(:,:) = (0.d0,0.d0)
        hamiltonian_GW_temp(:,:) = (0.d0,0.d0)
          do j_basis = 1, n_basis, 1
            do i_basis = 1, n_basis, 1

                 hamiltonian_GW_temp(i_basis,j_basis) = &
                 hamiltonian_GW_temp(i_basis,j_basis) + &
                 full_hamiltonian(i_basis,j_basis,1,i_k_point) - &
                 !DCMPLX(hartree_pot_LDA(:,:)) - &
                  DCMPLX(on_site_xc_matr(i_basis,j_basis))+ &
                  DCMPLX(loc_self_enrg(i_basis,j_basis))
                  !aux_matr_1(i_basis,j_basis)

            enddo
           enddo


      if (allocated(aux_matr_1))then
          deallocate(aux_matr_1)
       endif


!                 call multiply_ovlp_matr(hamiltonian_GW_temp(:,:),&
!                                         inv_full_ovlp_matrix_sqrt(:,:,i_k_point),&
!                                         hamiltonian_GW(:,:))

          do j_basis = 1, n_basis, 1
            do i_basis = 1, n_basis, 1
          !do j_basis = 1, n_states, 1
          !  do i_basis = 1, n_states, 1

                 DMFT_hamiltonian_k(i_basis,j_basis,i_k_point) = &
                 DMFT_hamiltonian_k(i_basis,j_basis,i_k_point)+ &
                 hamiltonian_GW_temp(i_basis,j_basis)

            enddo
           enddo




endif

       if ( allocated(gf_non_loc_temp))then
          deallocate(gf_non_loc_temp)
       endif


      if (allocated(KS_eigenvector_tmp))then
          deallocate(KS_eigenvector_tmp)
       endif


enddo
!stop
               if(scdmft_converged)then
!               if(.false.)then
                  call qp_spectrum_dmft (loc_self_enrg_GW(:,:,:), &
                                         DMFT_hamiltonian_k(:,:,:),&
                                         full_ovlp_matrix(:,:,:),&
                                         full_ovlp_matrix(:,:,:))
                                         !full_hamiltonian(:,:,1,:))
               endif



        on_site_gf(:,:,:)=(1./(n_k_points))*on_site_gf(:,:,:)


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
      if (.not. allocated( embed_gf_time))then
          allocate( embed_gf_time(n_basis,n_basis,-ntau:ntau))
       endif

             call transform_G (on_site_gf(:,:,:) , n_basis, &
               n_basis, embed_gf_time(:,:,:))

            dens_matr(:,:)=embed_gf_time(:,:,0)


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


       if ( allocated(gf_non_loc_temp))then
          deallocate(gf_non_loc_temp)
       endif



       if ( allocated(inv_gf_non_loc))then
          deallocate(inv_gf_non_loc)
       endif


       if (allocated(diag_ham_k))then
          deallocate(diag_ham_k)
       endif

       if ( allocated(diag_ham_k_LDA))then
          deallocate(diag_ham_k_LDA)
       endif






  end subroutine get_NEW_on_site_gf_freq
