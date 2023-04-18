  !subroutine test_loc_gf(overlap_matrix, hamiltonian, inv_gf_loc) ! (KS_eigenvalue)
  subroutine get_on_site_gf_freq_no_ovlp(full_hamiltonian,&
                                 full_ovlp_matrix,&
                                 inv_on_site_gf,&
                                 on_site_gf,&
                                 on_site_xc_matr,&!gf_non_loc,&
                                 loc_self_enrg,&
                                 number_reiterations,&
                                 max_reiteration,&!inv_k_summed_overlap_matr,&
                                 hartree_pot_LDA,&
                                 on_site_lda_hamiltonian,&
                                 on_site_xc_matr_no_ovlp,&
                                 full_ovlp_matrix_sqrt,&
                                 inv_full_ovlp_matrix_sqrt,&
                                 !loc_self_enrg_no_orth,&
                                 non_loc_self_enrg_GW,&
                                 hartree_pot_LDA_no_orth,& ! (KS_eigenvalue)
                                 !dens_matr_LDA, coeff_dens_matr_LDA_coeff,& ! (KS_eigenvalue)
                                 on_site_PBE_xc_matr,& ! (KS_eigenvalue)
                                 on_site_PBE_xc_matr_no_ovlp,& ! (KS_eigenvalue)
                                 inv_full_ovlp_matrix,&
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
  use poles_fit 
  use gt
  use timing

   implicit none


! local parameters
   integer l,a,b,n, i_a, j_a, i_k,i 
   integer hamiltonian_size
   integer ovlp_matrix_size
   integer i_basis, i_basis_1, i_matrix_size, j_matrix_size
   integer j_basis
   integer k_basis
   integer i_count
   integer id_root
   integer i_k_point_local 
   integer i_state
   integer i_freq
   integer i_k_point
   integer ham_size
   character*2 iter
   character*2 iter_1
   character*15 filename
   character*20 filename_RE
   character*20 filename_IM
 
! subroutine actual parameters
   complex*16, dimension(n_basis,n_basis,n_k_points) :: full_ovlp_matrix 
   complex*16, dimension(n_basis,n_basis,n_k_points) :: inv_full_ovlp_matrix 
!   complex*16, dimension(n_basis,n_basis,n_k_points) :: inv_full_ovlp_matrix 
   complex*16, dimension(n_basis,n_basis,n_k_points) :: full_ovlp_matrix_sqrt 
   complex*16, dimension(n_basis,n_basis,n_k_points) :: inv_full_ovlp_matrix_sqrt 
   !complex*16, dimension(n_basis,n_basis,n_k_points) :: loc_self_enrg_no_orth 
   complex*16, dimension(n_basis,n_basis,nomega,n_k_points) :: non_loc_self_enrg_GW 
!   complex*16, dimension(n_basis,n_basis,nomega) :: non_loc_self_enrg_GW 
   !complex*16, dimension(n_basis,n_basis,nomega) :: non_loc_self_enrg_GW 
   complex*16, dimension(n_basis,n_basis,n_k_points) :: hartree_pot_LDA_no_orth
   complex*16, dimension(n_basis,n_basis, n_spin,n_k_points) :: full_hamiltonian
   real*8, dimension(n_basis,n_basis) :: on_site_lda_hamiltonian 
   real*8, dimension(n_basis,n_basis) :: on_site_xc_matr
   real*8, dimension(n_basis,n_basis) :: loc_self_enrg 
   real*8, dimension(n_basis,n_basis) :: inv_k_summed_overlap_matr
   real*8, dimension(n_basis,n_basis) :: k_summed_overlap_matr 
   complex*16, dimension(n_basis,n_basis,nomega) :: inv_on_site_gf
   complex*16, dimension(n_basis,n_basis,nomega) :: inv_gf_test
   complex*16, dimension(n_basis,n_basis,nomega) :: on_site_gf
!   complex*16, dimension(n_basis,n_basis,nomega) :: on_site_gf_LDA
   !complex*16, dimension(n_basis,n_basis,n_k_points,nomega) :: gf_non_loc_LDA
   complex*16, dimension(:,:,:), allocatable :: gf_non_loc
   !complex*16, dimension(n_basis,n_basis,n_k_points,nomega) :: gf_non_loc_broad



   complex*16, dimension(:,:,:), allocatable :: inv_gf_non_loc 
   !complex*16, dimension(:,:,:,:), allocatable :: inv_gf_non_loc_temp
   complex*16, dimension(:,:,:), allocatable :: gf_non_loc_temp 
   complex*16, dimension(:,:,:), allocatable :: ovlp_gf_non_loc 
   integer info
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 
   real*8, dimension(n_basis,n_basis) :: on_site_xc_matr_no_ovlp
   !complex*16, dimension(n_basis,n_basis) :: on_site_xc_matr_no_ovlp
   complex*16, dimension(n_basis,n_basis,n_k_points) ::full_free_hamiltonian
!   complex*16, dimension(n_basis,n_basis,1,n_k_points) ::full_k_xc_matr


   logical  :: output
!   logical  :: scdmft_converged
   logical  :: output_spectral
   integer number_reiterations
   integer max_reiteration 



      integer aux_nomega
      real*8  aux_omegamax
      real*8  new_chemical_potential
      real*8 , dimension (:), allocatable :: aux_omega
      real*8 , dimension (:), allocatable :: aux_womega
      real*8 , dimension (:), allocatable :: spectrum 
      real*8 , dimension (:,:), allocatable :: spectrum_k
!      real*8 , dimension (:,:), allocatable :: on_site_xc_matr_no_ovlp
      complex*16, dimension(:,:,:), allocatable:: green_fn_par
      real*8, dimension(n_basis,n_basis) :: hartree_pot_LDA

            real*8 new_chem_pot
      complex*16, dimension(:,:,:), allocatable:: aux_ovlp
      complex*16, dimension(:,:,:), allocatable:: aux_new_full_ovlp_matrix
      complex*16, dimension(:,:,:), allocatable:: new_full_ovlp_matrix
      complex*16, dimension(:,:,:), allocatable:: diag_ham_k
      real*8, dimension(n_basis,n_basis) :: on_site_PBE_xc_matr
      real*8, dimension(n_basis,n_basis) :: on_site_PBE_xc_matr_no_ovlp
!   complex*16, dimension(:,:), allocatable:: dens_matr_k_LDA
   !complex*16, dimension(n_states,n_states,n_k_points):: coeff_dens_matr_LDA_coeff
!   complex*16, dimension(n_basis,n_basis,n_k_points):: coeff_dens_matr_LDA_coeff
   complex*16, dimension(:,:), allocatable:: coeff_product
   !real*8 , dimension(n_states,n_states):: dens_matr_LDA
!   real*8 , dimension(n_basis,n_basis):: dens_matr_LDA

!-----------------------------------------------------   
 
       
!       if (.not. allocated(ovlp_gf_non_loc))then
!          allocate(ovlp_gf_non_loc(n_basis,n_basis,n_k_points,nomega))
!       endif
!       if (.not. allocated(gf_non_loc_temp))then
!          allocate(gf_non_loc_temp(n_basis,n_basis,n_k_points,nomega))
!       endif
!       if (.not. allocated(gf_non_loc)then
!          allocate(gf_non_loc(n_basis,n_basis,nomega))
!       endif
!       if (.not. allocated(inv_gf_non_loc))then
!          allocate(inv_gf_non_loc(n_basis,n_basis,nomega))
!       endif

!       if (.not. allocated(myid_tmp))then
!          allocate(myid_tmp(n_k_points))
!       endif



       if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif


!          if(number_reiterations.eq.0) then
!              hartree_pot_LDA_no_orth(:,:,:) = (0.d0,0.d0)
!!              loc_self_enrg(:,:) = on_site_xc_matr_no_ovlp(:,:)
!           do i_k_point = 1, n_k_points
!              loc_self_enrg_no_orth(:,:,i_k_point) = DCMPLX(on_site_xc_matr_no_ovlp(:,:))
!           enddo
!          endif

!chemical_potential = -0.30830507
!chemical_potential = -0.10830507
!       do i_basis_1 = 1, n_basis, 1
!            write(use_unit,*)'from on-site GF',  &
!              hartree_pot_LDA(i_basis_1, i_basis_1), loc_self_enrg(i_basis_1, i_basis_1)!, on_site_xc_matr(i_basis_1, i_basis_1), on_site_PBE_xc_matr(i_basis_1, i_basis_1)
!       enddo

!     do i_k_point = 1, n_k_points
!        if (use_scalapack) then
!write(use_unit,*) 'new_chemical_potential', new_chemical_potential, DCMPLX(new_chemical_potential)
!           do i=0, n_tasks-1
!                  if (i_k_point == i*n_k_points/n_tasks + 1) then
!                     myid_tmp(i_k_point) = i
!                 exit
!                  endif
!           enddo
!
!        else

!           myid_tmp(i_k_point) = MOD(i_k_point, n_tasks)
!
!       endif
!   enddo

!    dens_matr_LDA(:,:) = 0.d0

    on_site_gf(:,:,:)=(0.d0,0.d0)
    inv_on_site_gf(:,:,:)=(0.d0,0.d0)
    i_k =0
  do i_k_point = 1, n_k_points, 1 


       if (.not. allocated(inv_gf_non_loc))then
          allocate(inv_gf_non_loc(n_basis,n_basis,nomega))
       endif



        inv_gf_non_loc(:,:,:) = (0.d0,0.d0)

!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
          if(number_reiterations.eq.0) then

       do i_freq = 1, nomega, 1 
          do i_basis = 1, n_basis, 1 
            do j_basis = 1, n_basis, 1 

                
                  inv_gf_non_loc(i_basis,j_basis,i_freq) = &
                  inv_gf_non_loc(i_basis,j_basis,i_freq) + &
                  (((0.d0,1.d0)*omega(i_freq)+DCMPLX(chemical_potential))* &
                  full_ovlp_matrix(i_basis,j_basis,i_k_point)- &
                  full_hamiltonian(i_basis,j_basis,1,i_k_point))



            enddo
           enddo
       enddo

    else
 
       do i_freq = 1, nomega, 1
          do i_basis = 1, n_basis, 1
            do j_basis = 1, n_basis, 1


                  inv_gf_non_loc(i_basis,j_basis,i_freq) = &
                  inv_gf_non_loc(i_basis,j_basis,i_freq) + &
                  (((0.d0,1.d0)*omega(i_freq)+DCMPLX(new_chemical_potential))* &
!                  (((0.d0,1.d0)*omega(i_freq)+DCMPLX(chemical_potential))* &
                  full_ovlp_matrix(i_basis,j_basis,i_k_point)- &
                  full_hamiltonian(i_basis,j_basis,1,i_k_point)+ &
                  DCMPLX((hartree_pot_LDA(i_basis,j_basis))) + &
                  DCMPLX(on_site_xc_matr_no_ovlp(i_basis,j_basis)) - &
                  DCMPLX(loc_self_enrg(i_basis,j_basis))-non_loc_self_enrg_GW(i_basis,j_basis,i_freq,i_k_point))!- &
!                  DCMPLX(loc_self_enrg(i_basis,j_basis))-non_loc_self_enrg_GW(i_basis,j_basis,i_freq))!- &
!                  (hartree_pot_LDA_no_orth(i_basis,j_basis,i_k_point)) + &
!                  DCMPLX(on_site_xc_matr_no_ovlp(i_basis,j_basis)) - &
!                  loc_self_enrg_no_orth(i_basis,j_basis,i_k_point)- &
!                  DCMPLX(on_site_PBE_xc_matr_no_ovlp(i_basis,j_basis)) )


            enddo
           enddo
       enddo

endif

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------




       if (.not. allocated(gf_non_loc_temp))then
          allocate(gf_non_loc_temp(n_basis,n_basis,nomega))
       endif
!       if (.not. allocated(gf_non_loc_no_ovlp_temp))then
!          allocate(gf_non_loc_no_ovlp_temp(n_basis,n_basis,nomega))
!       endif

      gf_non_loc_temp(:,:,:) = inv_gf_non_loc(:,:,:)


! at each frequency point evaluate the LU factorization, and check for errors
       do i_freq = 1, nomega, 1
!       do i_k_point = 1, n_k_points, 1
         call zgetrf( n_basis, n_basis, gf_non_loc_temp(:,:,i_freq), &
                    n_basis , ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then
            write(use_unit,*) " * Failure of LU decomposition at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
  !          stop
          endif
        endif
! this inverts the matrix in the LU factorized form
        call zgetri(n_basis, gf_non_loc_temp(:,:,i_freq), n_basis, &
                 ipiv, work, n_basis, info)

      if (info.ne.0) then
         if(myid.eq.0)then
            write(use_unit,*) " * Failure of matrix inversion at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
          endif
        endif
  
      enddo
       if ( allocated(inv_gf_non_loc))then
          deallocate(inv_gf_non_loc)
       endif


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


       if (.not. allocated(ovlp_gf_non_loc))then
          allocate(ovlp_gf_non_loc(n_basis,n_basis,nomega))
       endif
       if (.not. allocated(gf_non_loc))then
          allocate(gf_non_loc(n_basis,n_basis,nomega))
       endif


!do i_freq = 1, nomega


       on_site_gf(:,:,:) = on_site_gf(:,:,:) + gf_non_loc_temp(:,:,:)!*k_weights(i_k_point) 


!enddo
!       on_site_gf_no_ovlp(:,:,:) = on_site_gf_no_ovlp(:,:,:) + gf_non_loc_no_ovlp_temp(:,:,:)!*k_weights(i_k_point) 




       if (allocated(ovlp_gf_non_loc))then
          deallocate(ovlp_gf_non_loc)
       endif
       if (allocated(gf_non_loc_temp))then
          deallocate(gf_non_loc_temp)
       endif
!do i_basis =1 , n_states

!if(myid.eq.0) write(use_unit,*) 'coeff_dens_matr_LDA_coeff',i_basis,i_k_point, coeff_dens_matr_LDA_coeff(i_basis,i_basis,i_k_point) 

!enddo
 enddo

       !dens_matr_LDA(:,:)= dens_matr_LDA(:,:)*(1./n_k_points)


    !   if(use_mpi) then
    !     call sync_matrix(dens_matr_LDA(:,:), &
    !               n_states, n_states)
    !   endif
       !dens_matr_LDA(:,:)= dens_matr_LDA(:,:)*(1./n_k_points)

       on_site_gf(:,:,:)=(1./(n_k_points))*on_site_gf(:,:,:)


!          if(number_reiterations.eq.0) then
!             on_site_gf_LDA(:,:,:) = on_site_gf(:,:,:)
!          endif


       if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif


        inv_on_site_gf(:,:,:) = on_site_gf(:,:,:)


! at each frequency point evaluate the LU factorization, and check for errors
     do i_freq = 1, nomega, 1

         call zgetrf( n_basis, n_basis, inv_on_site_gf(:,:,i_freq), &
                    n_basis , ipiv, info )

        if (info.ne.0) then
          if (myid.eq.0)then
            write(use_unit,*) " * Failure of LU decomposition at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
  !          stop
          endif
        endif
! this inverts the matrix in the LU factorized form
        call zgetri(n_basis, inv_on_site_gf(:,:,i_freq), n_basis, &
                 ipiv, work, n_basis, info)

      if (info.ne.0) then
         if(myid.eq.0)then
            write(use_unit,*) " * Failure of matrix inversion at",&
                            " frequency point " , i_freq
            write(use_unit,*) " * Error info = ", info
          endif
        endif

      enddo
       if ( allocated(ipiv))then
          deallocate(ipiv)
       endif

       if (allocated(work))then
          deallocate(work )
       endif



!do i_basis = 1, n_basis
!write(use_unit,*)  inv_on_site_gf(i_basis,i_basis,100), on_site_gf(i_basis,i_basis,100)  
!enddo
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

!if(number_reiterations.eq.1) stop


       if (allocated(gf_non_loc))then
          deallocate(gf_non_loc)
       endif
       if (allocated(ovlp_gf_non_loc))then
          deallocate(ovlp_gf_non_loc)
       endif
       if (allocated(gf_non_loc_temp))then
          deallocate(gf_non_loc_temp)
       endif
       if ( allocated(inv_gf_non_loc))then
          deallocate(inv_gf_non_loc)
       endif

       if ( allocated(ipiv))then
          deallocate(ipiv)
       endif

       if (allocated(work))then
          deallocate(work )
       endif





         if(myid.eq.0)then

     output = .false.
      if(output) then
           ! do a = 1, n_freq, 1
          if( i_a.lt.10 ) then
             write(iter_1,'(A,I1)') "0",number_reiterations
          else
             write(iter_1,'(I2)') number_reiterations
          endif
         do i_a = 1, n_basis, 1
          if( i_a.lt.10 ) then
             write(iter,'(A,I1)') "0",i_a
          else
             write(iter,'(I2)') i_a
          endif

          filename_RE = "gree_fr_"//iter_1//"RE_"//iter//".dat"
          filename_IM = "gree_fr_"//iter_1//"IM_"//iter//".dat"
        !    do a = 1, n_freq, 1
          open(77, file=filename_RE)
          open(75, file=filename_IM)
           do a = 1, nomega, 1
           !do j_a = 1, n_basis, 1
          !  do a = 1, n_freq, 1
              write(77,*) omega(a), &
                       real(on_site_gf(i_a,i_a,a))
              write(75,*) omega(a), &
                      aimag(on_site_gf(i_a,i_a,a))
            ! enddo
            enddo
          close(77)
          close(75)
        enddo
      endif




    endif

  end subroutine get_on_site_gf_freq_no_ovlp 
