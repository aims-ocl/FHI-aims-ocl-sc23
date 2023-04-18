       subroutine get_qp_spectrum_dmft            &
           (anacon_type,n_max_par, &
            n_low_state,n_high_state, &
            n_freq, omega, &
            sigma_par_loc, DMFT_ham_KS, &
            ovlp_matr, &
            filename )


              use constants, only: one_over_sqrt2, pi
              use dimensions
              use runtime_choices
              use species_data
              use physics
              use pbc_lists
              use mpi_tasks, only: myid
              use localorb_io, only: OL_high, output_priority
!              use gw_para
 
              implicit none


!  ARGUMENTS

      integer :: anacon_type
      integer :: n_max_par
      integer :: n_low_state, n_high_state
      integer :: n_freq
      integer :: i_freq, i_freq_1, i_freq_2
      integer :: sign_freq

      real*8 :: omega(n_freq)
      !complex*16 :: DMFT_ham_KS(n_states,n_states,n_k_points)
      complex*16 :: DMFT_ham_KS(n_basis,n_basis,n_k_points)
      complex*16 :: ovlp_matr(n_basis,n_basis,n_k_points)
      character*17 filename
!      complex*16 :: sigma_par_loc(n_max_par,n_states,n_states)
      complex*16 :: sigma_par_loc(n_max_par,n_states,n_states,n_k_points)
!      complex*16 :: sigma_par_loc(n_max_par,n_basis,n_basis,n_k_points)
!      complex*16 :: sigma_par_loc(n_max_par,n_basis,n_basis)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8  :: e_diff
      real*8  :: de
      integer :: n,l
      real*8  :: en
      real*8  :: w
      real*8  :: w1
      real*8  :: w2
      real*8  :: mu
      real*8  :: delta_mu
      integer :: output_priority_old
 
      real*8  :: ims, res

      complex*16  :: selfe, dselfe

      integer :: i_state,i_state_1
      integer :: i, i_count
      integer :: i_spin, i_k_point
      integer :: i_e
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
      real*8, dimension(:,:), allocatable::  spectral_function
      real*8, dimension(:), allocatable::  summed_spectral_function
      real*8, dimension(:), allocatable::  dos_dmft

      integer aux_nomega, n_nonsingular
      real*8  aux_omegamax
      real*8 , dimension (:), allocatable :: aux_omega
      real*8 , dimension (:), allocatable :: aux_womega
      real*8 , dimension (:,:), allocatable :: dmft_eigenvalues
      complex*16 , dimension (:,:,:), allocatable :: dmft_eigenvectors
      complex*16 , dimension (:), allocatable :: ham_triangle
      complex*16 , dimension (:), allocatable :: ovlp_triangle


      aux_nomega = 15000
      aux_omegamax = 6.50
      allocate(aux_omega(aux_nomega))
      allocate(aux_womega(aux_nomega))
      allocate(spectral_function(-aux_nomega:aux_nomega,n_k_points))
      allocate(summed_spectral_function(-aux_nomega:aux_nomega))
      allocate(dos_dmft(-aux_nomega:aux_nomega))
      allocate(dmft_eigenvalues(n_states,n_k_points))
      allocate(dmft_eigenvectors(n_states,n_states,n_k_points))
      !allocate(dmft_eigenvalues(n_basis,n_k_points))
      !allocate(dmft_eigenvectors(n_basis,n_basis,n_k_points))
      !allocate(summed_spectral_function(dos_n_en_points))

      call gauleg(0d0, aux_omegamax, aux_omega, aux_womega, aux_nomega)

!      if (n_spin.ne.1)then
!        write(use_unit,*) " Quasi particle spectrum ",&
!                      " implemented only for spin unpolarized"
!      endif

      
        if(.not.allocated(spectral_function))then
          allocate(spectral_function(n_omega,n_k_points))
        endif


        if(myid.eq.0 )then
             open(86,file=filename)
        endif       

           do i_k_point = 1, n_k_points,1


            allocate(ham_triangle   (n_basis*(n_basis+1)/2))
            allocate(ovlp_triangle   (n_basis*(n_basis+1)/2))


                    n = 0
                        do l = 1, n_basis

                           ham_triangle(n+1:n+l) = DMFT_ham_KS(1:l,l,i_k_point)
                           ham_triangle(n+1:n+l-1) = &
                           dconjg(DMFT_ham_KS(l,1:l-1,i_k_point))

                           ovlp_triangle(n+1:n+l) = ovlp_matr(1:l,l,i_k_point)
                           ovlp_triangle(n+1:n+l-1) = &
                           dconjg(ovlp_matr(l,1:l-1,i_k_point))

                           n = n+l
                        enddo



                output_priority_old = output_priority
                output_priority = OL_high

                   call improve_complex_eigenfunctions &
                        ( ovlp_triangle, &
                          ham_triangle(:),  &
                          dmft_eigenvalues(:,i_k_point), &
                          dmft_eigenvectors (:,:,1),1)

                output_priority = output_priority_old

            deallocate(ham_triangle)
            deallocate(ovlp_triangle)


           enddo
      
         spectral_function (:,:) = 0.d0
         summed_spectral_function(:)=0.d0


       do sign_freq = -1, 1, 2
        do i_freq = 1, aux_nomega, 1

         w = sign_freq*aux_omega(i_freq)

         do i_state = 1, n_states,1
!         do i_state = 1, n_basis,1
           do i_k_point = 1, n_k_points,1

!if(.false.) &
                                     call get_real_selfenergy (&
                                                   anacon_type,&
                                                       n_freq ,&
                                                        omega ,&
                                               dcmplx(w, 0.d0),&
                                                    n_max_par ,&
          sigma_par_loc(1:n_max_par,i_state,i_state,i_k_point),&
                                                          selfe)
!                       sigma_par_loc(1:n_max_par,i_state,i_state),selfe)


           res = real(selfe)
           ims = aimag(selfe)

!           do i_k_point = 1, n_k_points,1
              spectral_function(sign_freq*i_freq,i_k_point)  = &
              spectral_function(sign_freq*i_freq,i_k_point)    &
          + 1.d0/pi * abs(ims) &
             /((w  - dmft_eigenvalues(i_state,i_k_point)-&
               (res))**2 &
             + (ims)**2 )

                      summed_spectral_function(sign_freq*i_freq) = &
                      summed_spectral_function(sign_freq*i_freq)+ &
                      spectral_function(sign_freq*i_freq,i_k_point)
          enddo
        enddo

       enddo  ! i_freq
      enddo !sign_freq

        summed_spectral_function(:)=summed_spectral_function(:)*&
                                    (1.d0/n_k_points)*(1.d0/n_states)




       de= (dos_high_energy - dos_low_energy)/dble(dos_n_en_points-1)
       dos_dmft(:) =0d0

     !if(.false.) then
      if(.true.) then
       do sign_freq = -1, 1, 2
        do i_freq_1 = 1, aux_nomega, 1

         w1 = sign_freq*aux_omega(i_freq_1)
       if(.true.) then
             do i_freq_2 = 1, aux_nomega, 1

          w2 = sign_freq*aux_omega(i_freq_2)


                 dos_dmft (sign_freq*i_freq_1) =  dos_dmft (sign_freq*i_freq_1) + &
                      (1/sqrt(pi))*(one_over_sqrt2/dos_alpha)&
                      *aux_womega(i_freq_2) &
                      *(exp((-1d0)*(w1-w2)&
                      **2*(one_over_sqrt2/dos_alpha)**2))&
                      * summed_spectral_function(sign_freq*i_freq_2)!/(2d0*dE)

           enddo
        endif 


       if(.false.) then
             do i_state = 1, n_states, 1
             do i_k_point = 1, n_k_points, 1

                 dos_dmft (sign_freq*i_freq_1) =  dos_dmft (sign_freq*i_freq_1) + &
                      (1/sqrt(pi))*(one_over_sqrt2/dos_alpha)&
                      *aux_womega(i_freq_1) &
                      *(exp((-1d0)*(w1-dmft_eigenvalues(i_state,i_k_point)&!+&
                      )**2*(one_over_sqrt2/dos_alpha)**2)&
                      )!/(2d0*dE)

             enddo
            enddo
       endif

           enddo
          enddo
        endif

     if(.true.) then
       do sign_freq = -1, 1, 2
        do i_freq = 1, aux_nomega, 1

         w = sign_freq*aux_omega(i_freq)

         if(myid.eq.0)then
             write(86,*) real(w)*hartree, dos_dmft(sign_freq*i_freq)
           endif 
         enddo
        enddo
     endif


     if(.false.) then
!     if(.true.) then
        do i_e = 1, dos_n_en_points ,1
           en = dos_low_energy + dble(i_e-1)*de
         if(myid.eq.0)then
             write(86,*) en, summed_spectral_function(i_e) 
           endif
        enddo
     endif

        if(myid.eq.0)then
             close(86)
        endif
 
!      enddo!i_spin

      end subroutine get_qp_spectrum_dmft


