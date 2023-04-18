  subroutine get_NEW_on_site_gf_freq_PBE0(full_hamiltonian,&
                                     full_ovlp_matrix,&
                                     on_site_gf,&
                                     on_site_xc_matr,&
                                     loc_self_enrg, &
                                     hartree_pot_LDA,& 
                                     number_reiterations,&
                                     scdmft_converged,&
                                     full_ovlp_matrix_sqrt, &
                                     on_site_PBE_xc_matr,&
                                     dens_matr,&
                                     new_chemical_potential,&
                                     new_chemical_potential_cluster,&
                                     embed_part_number_LDA)!!,&





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
   complex*16, dimension(n_basis,n_basis,n_k_points), intent(in) ::&
    full_ovlp_matrix 
   complex*16, dimension(n_basis,n_basis, n_spin,n_k_points), intent(in) ::&
    full_hamiltonian
   real*8, dimension(n_basis,n_basis),intent(in) :: loc_self_enrg 
   real*8, dimension(n_basis,n_basis),intent(in) :: on_site_xc_matr
   real*8  embed_part_number_LDA
   complex*16, dimension(n_basis,n_basis,nomega), intent(out) :: on_site_gf
   complex*16, dimension(:,:,:), allocatable :: gf_non_loc
   real*8, dimension(:,:,:), allocatable :: embed_gf_time

   complex*16, dimension(:,:,:), allocatable :: inv_gf_non_loc 
   complex*16, dimension(:,:,:), allocatable :: gf_non_loc_temp 
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


      complex*16, dimension(n_basis,n_basis,n_k_points) :: &
      full_ovlp_matrix_sqrt
      complex*16, dimension(:,:), allocatable:: diag_ham
      complex*16, dimension(:,:), allocatable:: diag_ham_k
      complex*16, dimension(:,:), allocatable:: diag_ham_k_LDA
      complex*16, dimension(:,:), allocatable:: DMFT_hamiltonian
      real*8 , dimension (:,:), allocatable :: hartree_pot_LDA_spec
      complex*16, dimension(:,:), allocatable :: ovlp_full_hamiltonian
      complex*16, dimension(:,:), allocatable :: full_hamiltonian_orth
      complex*16, dimension(:,:), allocatable :: ovlp_DMFT_hamiltonian
      real*8, dimension(n_basis,n_basis) :: on_site_PBE_xc_matr
      complex*16, dimension(:,:), allocatable:: aux_matr_1
      complex*16, dimension(:,:), allocatable:: loc_self_enrg_temp
      complex*16, dimension(:,:), allocatable:: on_site_xc_matr_temp
      complex*16, dimension(:,:,:), allocatable:: KS_eigenvector_tmp
      real*8 , dimension(n_basis,n_basis):: dens_matr
      real*8, dimension(:), allocatable :: DMFT_dos
      real*8, dimension(:), allocatable :: DMFT_dos_LDA
      real*8, dimension(:), allocatable :: DMFT_dos_k
      real*8, dimension(:), allocatable :: DMFT_dos_k_LDA


            real*8 new_chemical_potential
            real*8 new_chemical_potential_cluster


!------------------------------------------------------   


       if (.not. allocated(loc_self_enrg_temp))then
          allocate(loc_self_enrg_temp(n_basis,n_basis))
       endif

    on_site_gf(:,:,:)=0.d0
    dens_matr (:,:) = 0.d0
    loc_self_enrg_temp(:,:) = DCMPLX(loc_self_enrg(:,:)) - &
                              DCMPLX(on_site_PBE_xc_matr(:,:))
    i_k =0

 do i_k_point = 1, n_k_points


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



!----------------------------------------------------------------------------


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
                  DCMPLX(loc_self_enrg(i_basis,j_basis)))
                  !aux_matr_1(i_basis,j_basis)- &
                  !DCMPLX(on_site_PBE_xc_matr(i_basis,j_basis)))


            enddo
           enddo
       enddo

!       if (allocated(aux_matr_1))then
!          deallocate(aux_matr_1)
!       endif

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

       if ( allocated(inv_gf_non_loc))then
          deallocate(inv_gf_non_loc)
       endif


!       if (.not. allocated(gf_non_loc_temp))then
!          allocate(gf_non_loc_temp(n_basis,n_basis,nomega))
!       endif
!          gf_non_loc_temp(:,:,:) = inv_gf_non_loc(:,:,:)
! at each frequency point evaluate the LU factorization, and check for errors
       do i_freq = 1, nomega, 1

       if (.not. allocated(ipiv))then
          allocate(ipiv(n_states))
       endif

       if (.not. allocated(work))then
          allocate(work (n_states))
       endif


         call zgetrf( n_states, n_states, gf_non_loc_temp(:,:,i_freq), &
                    n_states , ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then
            write(use_unit,*) " * Failure of LU decomposition at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
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
  

       if ( allocated(ipiv))then
          deallocate(ipiv)
       endif

       if ( allocated(work))then
          deallocate(work)
       endif


       enddo


!       if ( allocated(inv_gf_non_loc))then
!          deallocate(inv_gf_non_loc)
!       endif
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



  enddo

      if(allocated(KS_eigenvector_tmp))then
        deallocate(KS_eigenvector_tmp)
      endif

  enddo



        on_site_gf(:,:,:)=(1./(n_k_points))*on_site_gf(:,:,:)


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
      if (.not. allocated( embed_gf_time))then
          allocate( embed_gf_time(n_basis,n_basis,-ntau:ntau))
       endif

             call transform_G (on_site_gf(:,:,:) , n_basis, &
               n_basis, embed_gf_time(:,:,:))

            dens_matr(:,:)=embed_gf_time(:,:,0)

      if (allocated( embed_gf_time))then
          deallocate( embed_gf_time)
       endif


output_spectral=.false.
!output_spectral=.true.
if(scdmft_converged)then
!if(.true.)then


!output_spectral=.true.
if (output_spectral)then
!if(scdmft_converged)then


       if (.not. allocated(diag_ham_k))then
          allocate(diag_ham_k(n_states, n_states))
       endif

       if (.not. allocated(diag_ham_k_LDA))then
          allocate(diag_ham_k_LDA(n_basis,n_basis))
       endif

       if (.not. allocated(diag_ham))then
          allocate(diag_ham(n_states, n_states))
       endif

       if (.not. allocated(DMFT_hamiltonian))then
          allocate(DMFT_hamiltonian(n_basis,n_basis))
       endif



       if(.not. allocated (DMFT_dos)) allocate(DMFT_dos(dos_n_en_points))
       if(.not. allocated (DMFT_dos_LDA)) allocate(DMFT_dos_LDA(dos_n_en_points))
       if(.not. allocated (DMFT_dos_k_LDA)) allocate(DMFT_dos_k_LDA(dos_n_en_points))
       if(.not. allocated (DMFT_dos_k)) allocate(DMFT_dos_k(dos_n_en_points))



      if (myid.eq.0) then
        write(use_unit,'(A)')
        write(use_unit,'(A)')"--------------------------------------------"
        write(use_unit,'(10X,A)') "Computation of the Spectral-functions using"
        write(use_unit,'(10X,A)') "analytic continuation starts ..."
      endif






!-----------------------------------------------------------------------------------------
      aux_nomega = 15000
      aux_omegamax = 8.50
      allocate(aux_omega(aux_nomega))
      allocate(aux_womega(aux_nomega))
      if(.not. allocated (spectrum_k)) allocate(spectrum_k(-aux_nomega:aux_nomega))
      if(.not. allocated (spectrum)) allocate(spectrum(-aux_nomega:aux_nomega))
      if(.not. allocated (spectrum_k_LDA)) allocate(spectrum_k_LDA(-aux_nomega:aux_nomega))
      if(.not. allocated (spectrum_LDA)) allocate(spectrum_LDA(-aux_nomega:aux_nomega))

              if(.not.allocated(green_fn_par))then
                !allocate(green_fn_par(n_max_par, n_states, n_states))
                allocate(green_fn_par(n_max_par, n_basis, n_basis))
              endif
              if(.not.allocated(green_fn_par_LDA))then
                !allocate(green_fn_par(n_max_par, n_states, n_states))
                allocate(green_fn_par_LDA(n_max_par, n_basis, n_basis))
              endif

      call gauleg(0d0, aux_omegamax, aux_omega, aux_womega, aux_nomega)




                 if (number_reiterations.lt.10)then
                   write(iter,'(A,I1)') "0", number_reiterations
                 else
                   write(iter,'(I2)') number_reiterations
                 endif

                 !if(n_spin .eq. 1 )then
                   filename = "sp_ImG_scDMFT.dat"
                   filename_LDA = "sp_ImG_LDA.dat"

spectrum(:) =0.d0
spectrum_LDA(:) =0.d0
diag_ham(:,:) =(0.d0,0.d0)
green_fn_par(:,:,:) = 0.d0
green_fn_par_LDA(:,:,:) = 0.d0
DMFT_dos(:) = 0.d0
DMFT_dos_LDA(:) = 0.d0
i_k = 0
do i_k_point = 1, n_k_points

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

DMFT_hamiltonian(:,:) = (0.d0,0.d0)
ovlp_full_hamiltonian(:,:)=(0.d0,0.d0)
full_hamiltonian_orth(:,:)=(0.d0,0.d0)
ovlp_DMFT_hamiltonian(:,:)=(0.d0,0.d0)
!DMFT_hamiltonian_orth(:,:)=(0.d0,0.d0)
spectrum_k(:) = 0.d0
spectrum_k_LDA(:) = 0.d0


if(.false.)then
          do i_basis = 1, n_basis, 1
            do j_basis = 1, n_basis, 1

                  DMFT_hamiltonian(i_basis,j_basis) = &
                   DMFT_hamiltonian(i_basis,j_basis) + &
                  (full_hamiltonian(i_basis,j_basis,1,i_k_point)-&
                  DCMPLX(on_site_xc_matr(i_basis,j_basis)) +&!- &
                  !(hartree_pot_LDA(i_basis,j_basis)) + &
                   DCMPLX(loc_self_enrg(i_basis,j_basis)))

            enddo
           enddo
endif


!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------



       if (.not. allocated(aux_matr_1))then
          allocate(aux_matr_1(n_basis,n_basis))
       endif


               call multiply_ovlp_matr(&
                                        loc_self_enrg_temp(:,:),&
                                       full_ovlp_matrix_sqrt(:,:,i_k_point),&
                                        aux_matr_1(:,:))

          do i_basis = 1, n_basis, 1
            do j_basis = 1, n_basis, 1

                  DMFT_hamiltonian(i_basis,j_basis) = &
                   DMFT_hamiltonian(i_basis,j_basis) + &
                  (full_hamiltonian(i_basis,j_basis,1,i_k_point)-&
                  DCMPLX(on_site_xc_matr(i_basis,j_basis)) +&!- &
                  !(on_site_xc_matr_temp(i_basis,j_basis)) +&!- &
                  ! (loc_self_enrg_temp(i_basis,j_basis)))
                    aux_matr_1(i_basis,j_basis)+&
                   DCMPLX(on_site_PBE_xc_matr(i_basis,j_basis)))
                  !DCMPLX(on_site_xc_matr(i_basis,j_basis)) - &
                  !(hartree_pot_LDA(i_basis,j_basis)) + &
                  !DCMPLX(loc_self_enrg(i_basis,j_basis)))

            enddo
           enddo

       if (allocated(aux_matr_1))then
          deallocate(aux_matr_1)
       endif


      if (myid.eq.0) then
        write(use_unit,'(A)')
        write(use_unit,'(A)')"------------------------------------------------------"
        write(use_unit,*) "Diagonalizing the Hamiltonians for k-point", i_k_point
      endif




!call spectrum_dmft ( DMFT_hamiltonian_orth(:,:),&
call spectrum_dmft ( DMFT_hamiltonian(:,:),&
                    spectrum_k(:), aux_omegamax, aux_omega, &
                    aux_womega, aux_nomega, green_fn_par, i_k_point,&
                    diag_ham_k(:,:),DMFT_dos_k, full_ovlp_matrix(:,:,i_k_point))





DMFT_dos(:) = DMFT_dos(:) + DMFT_dos_k(:)!
spectrum(:) = spectrum(:) + spectrum_k(:)!
diag_ham(:,:) = diag_ham(:,:) + diag_ham_k(:,:)


      if( allocated(KS_eigenvector_tmp))then
        deallocate(KS_eigenvector_tmp)
      endif

enddo
       if( allocated (DMFT_dos_k_LDA)) deallocate(DMFT_dos_k_LDA)
       if( allocated (DMFT_dos_k)) deallocate(DMFT_dos_k)


       if ( allocated(diag_ham_k))then
          deallocate(diag_ham_k)
       endif

       if ( allocated(diag_ham_k_LDA))then
          deallocate(diag_ham_k_LDA)
       endif
       if ( allocated(DMFT_hamiltonian))then
          deallocate(DMFT_hamiltonian)
       endif
       if ( allocated(ovlp_full_hamiltonian))then
          deallocate(ovlp_full_hamiltonian)
       endif
       if (allocated(full_hamiltonian_orth))then
          deallocate(full_hamiltonian_orth)
       endif
       if ( allocated(ovlp_DMFT_hamiltonian))then
          deallocate(ovlp_DMFT_hamiltonian)
       endif

           spectrum = (1./(n_k_points))&
                      *spectrum
           spectrum_LDA = (1./(n_k_points))&
                      *spectrum_LDA
           diag_ham(:,:)= 1./(n_k_points)&
                      *diag_ham(:,:)
         if(myid.eq.0)then

           do i_basis =1 , n_states
!do j_basis =1 , n_basis
              write(use_unit,*) 'diag_ham_dmft',i_basis,real(diag_ham(i_basis,i_basis))
!enddo
           enddo
          
         endif

       if ( allocated(diag_ham))then
          deallocate(diag_ham)
       endif


call get_spectrum_dmft(anacon_type, &
                       green_fn_par, &
                       n_max_par, &
                       omega, &
                       nomega, &
                       filename, &
                       n_basis, &
                       spectrum, &
                       aux_omegamax, &
                       aux_omega, &
                       aux_womega, &
                       aux_nomega, &
                       DMFT_dos_LDA,DMFT_dos,&
                       new_chemical_potential)


              if(allocated(green_fn_par))then
                deallocate(green_fn_par)
              endif
              if(allocated(green_fn_par_LDA))then
                deallocate(green_fn_par_LDA)
              endif
              if( allocated (DMFT_dos_LDA)) deallocate(DMFT_dos_LDA)
              if( allocated (DMFT_dos)) deallocate(DMFT_dos)

!stop
endif
endif
!stop
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
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

       if (allocated(loc_self_enrg_temp))then
          deallocate(loc_self_enrg_temp)
       endif





  end subroutine get_NEW_on_site_gf_freq_PBE0
