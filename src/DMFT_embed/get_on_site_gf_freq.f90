  subroutine get_on_site_gf_freq(full_hamiltonian,&
                                 full_ovlp_matrix,&
                                 inv_on_site_gf,&
                                 on_site_gf,&
                                 on_site_xc_matr,&!gf_non_loc,&
                                 loc_self_enrg,&
                                 number_reiterations,&
                                 max_reiteration,&!inv_k_summed_overlap_matr,&
                                 hartree_pot_LDA,&
                                 full_ovlp_matrix_sqrt,&
                                 !inv_full_ovlp_matrix_sqrt,&
                                 loc_self_enrg_GW,&
                                 dens_matr_LDA,& 
                                 on_site_PBE_xc_matr,&
                                 !inv_full_ovlp_matrix,& ! (KS_eigenvalue)
                                 new_chemical_potential)!,&
                                 !full_k_xc_matr) ! (KS_eigenvalue)




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
   integer l,a,b,n, i_a, j_a, i_k,i,i_q 
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
   integer info
 
! subroutine actual parameters
   complex*16, dimension(n_basis,n_basis,n_k_points) :: full_ovlp_matrix 
   complex*16, dimension(n_basis,n_basis,n_k_points) :: full_ovlp_matrix_sqrt 
   complex*16, dimension(n_basis,n_basis,nomega) :: loc_self_enrg_GW 
   complex*16, dimension(n_basis,n_basis, n_spin,n_k_points) ::&
   full_hamiltonian
   real*8, dimension(n_basis,n_basis) :: on_site_xc_matr
   real*8, dimension(n_basis,n_basis) :: loc_self_enrg 
   real*8, dimension(n_basis,n_basis) :: k_summed_overlap_matr 
   real*8, dimension(:,:,:), allocatable :: embed_gf_time 
   complex*16, dimension(n_basis,n_basis,nomega) :: inv_on_site_gf
   complex*16, dimension(n_basis,n_basis,nomega) :: on_site_gf
   complex*16, dimension(:,:,:), allocatable :: gf_non_loc
   complex*16, dimension(:,:,:), allocatable :: inv_gf_non_loc 
   complex*16, dimension(:,:,:), allocatable :: gf_non_loc_temp 
   complex*16, dimension(:), allocatable :: ipiv
   complex*16, dimension(:), allocatable :: work 
   complex*16, dimension(n_basis,n_basis,n_k_points) ::full_free_hamiltonian
!   complex*16, dimension(n_basis,n_basis,1,n_k_points) ::full_k_xc_matr


   logical  :: output
   logical  :: output_spectral
   integer number_reiterations
   integer max_reiteration 



   real*8, dimension(n_basis,n_basis) :: hartree_pot_LDA
   real*8  new_chemical_potential
   real*8 new_chem_pot
   real*8, dimension(n_basis,n_basis) :: on_site_PBE_xc_matr
   real*8 , dimension(n_basis,n_basis):: dens_matr_LDA
   complex*16, dimension(:,:), allocatable:: aux_matr_1
   complex*16, dimension(:,:,:), allocatable:: KS_eigenvector_tmp
!-----------------------------------------------------   
 
       
          if(number_reiterations.eq.0) then
              hartree_pot_LDA(:,:) = 0.d0
              loc_self_enrg(:,:) = on_site_xc_matr(:,:)
          endif


   if(number_reiterations.eq.0) then
    dens_matr_LDA(:,:) = 0.d0
   endif

    on_site_gf(:,:,:)=(0.d0,0.d0)
    i_k =0

  do i_k_point = 1, n_k_points, 1 

      if(.not. allocated(KS_eigenvector_tmp))then
        allocate(KS_eigenvector_tmp(n_basis,n_states,n_spin))
      endif

           if (myid .eq. MOD(i_k_point, n_tasks) .and. &
           myid <= n_k_points ) then

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
                  full_hamiltonian(i_basis,j_basis,1,i_k_point))!+ &


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
                  full_ovlp_matrix(i_basis,j_basis,i_k_point)- &
                  full_hamiltonian(i_basis,j_basis,1,i_k_point)+ &
                  !DCMPLX(hartree_pot_LDA(i_basis,j_basis)) - &
                  DCMPLX(on_site_xc_matr(i_basis,j_basis)) - &
                  DCMPLX(loc_self_enrg(i_basis,j_basis)) - &
                  loc_self_enrg_GW(i_basis,j_basis,i_freq))

            enddo
           enddo
       enddo

endif

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
if(.true.)then
!if(.false.)then

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
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------


!! at each frequency point evaluate the LU factorization, and check for errors
       do i_freq = 1, nomega, 1

       if (.not. allocated(ipiv))then
          allocate(ipiv(n_states))
!          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_states))
!          allocate(work (n_basis))
       endif

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

       if (allocated(ipiv))then
          deallocate(ipiv)
       endif

       if (allocated(work))then
          deallocate(work)
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
!-----------------------------------------------------------------------------------------


       if (.not. allocated(gf_non_loc))then
          allocate(gf_non_loc(n_basis,n_basis,nomega))
       endif


gf_non_loc(:,:,:)=(0.d0,0.d0)
do i_freq = 1, nomega



                call transform_to_NAO_basis_complex(&
                                      gf_non_loc_temp(:,:,i_freq),&
                                      i_k_point,&
                                      gf_non_loc(:,:,i_freq), &
                                      KS_eigenvector_tmp)



        on_site_gf(:,:,i_freq) = on_site_gf(:,:,i_freq) + &
                                 gf_non_loc(:,:,i_freq)!*k_weights(i_k_point) 


enddo




       if (allocated(gf_non_loc_temp))then
          deallocate(gf_non_loc_temp)
       endif

      if(allocated(KS_eigenvector_tmp))then
        deallocate(KS_eigenvector_tmp)
      endif

 enddo

   if(number_reiterations.eq.0) then
       dens_matr_LDA(:,:)= dens_matr_LDA(:,:)*(1./n_k_points)
   endif

       if (allocated(gf_non_loc_temp))then
          deallocate(gf_non_loc_temp)
       endif


       on_site_gf(:,:,:)=(1./(n_k_points))*on_site_gf(:,:,:)

  if(number_reiterations.eq.0) then

      if (.not. allocated( embed_gf_time))then
          allocate( embed_gf_time(n_basis,n_basis,-ntau:ntau))
       endif

             call transform_G ( on_site_gf(:,:,:) , n_basis, &
               n_basis, embed_gf_time(:,:,:))
            dens_matr_LDA(:,:)=embed_gf_time(:,:,0)


  endif


        inv_on_site_gf(:,:,:) = on_site_gf(:,:,:)

! at each frequency point evaluate the LU factorization, and check for errors
     do i_freq = 1, nomega, 1

       if (.not. allocated(ipiv))then
          allocate(ipiv(n_basis))
       endif

       if (.not. allocated(work))then
          allocate(work (n_basis))
       endif

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
  
       if (allocated(ipiv))then
          deallocate(ipiv)
       endif

       if (allocated(work))then
          deallocate(work)
       endif

      enddo

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

       if (allocated(gf_non_loc))then
          deallocate(gf_non_loc)
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
    ! output = .true.
      if(output) then
           ! do a = 1, n_freq, 1
!        do i_a = 1, n_basis, 1
!          if( i_a.lt.10 ) then
!             write(iter,'(A,I1)') "0",outter_loop_reiteration
!          else
!             write(iter,'(I2)') outter_loop_reiteration
!          endif
          if( number_reiterations.lt.10 ) then
             write(iter,'(A,I1)') "0", number_reiterations
          else
             write(iter,'(I2)') number_reiterations
          endif

          filename_RE = "on_site_RE_"//iter//".dat"
          filename_IM = "on_site_IM_"//iter//".dat"
        !    do a = 1, n_freq, 1
          open(77, file=filename_RE)
          open(75, file=filename_IM)
           do a = 1, nomega, 1
              write(77,*) omega(a), &
                       real(on_site_gf(1,1,a))
              write(75,*) omega(a), &
                      aimag(on_site_gf(1,1,a)) !aimag(1.d0/((0.d0,1.d0)*omega(a)))
            enddo
          close(77)
          close(75)
!            enddo
       ! enddo
      endif
!      endif
!      endif


    endif

  end subroutine get_on_site_gf_freq 
