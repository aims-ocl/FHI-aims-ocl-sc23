      subroutine get_spectrum_k_dmft (anacon_type,&
                                      green_fn_par,&
                                      n_max_par, &
                                      omega, nomega,&
                                      n_matrix, &
                                      spectrum,&
                                      aux_omegamax,&
                                      aux_omega,&
                                      aux_womega,&
                                      aux_nomega)

!get the spectrum from the (analytically continued) diagonalized Green's function


      use dimensions
      use runtime_choices
      use species_data
      use constants
      use mpi_tasks
      use synchronize_mpi
      use prodbas
      use physics

!ARGUMENTS
      implicit none
      integer n_matrix
      integer n_max_par
      integer nomega
      integer anacon_type
      complex*16 green_fn_par (n_max_par, n_matrix, n_matrix)
      complex*16 inv_overlap_matrix(n_basis,n_basis)
      real*8 omega(nomega)
      real*8 new_chem_pot

!internal
!      complex*16 diagonal_green (n_matrix)
      complex*16, dimension(:,:),allocatable :: aux_green 
      complex*16 denominator
      complex*16 numerator
!      real*8     w
      complex*16 w 
      complex*16 xdata (n_max_par)
      complex*16 gtmp
      real*8 n_particles
      real*8 n_particles_w
      real*8 n_particles_w_1
      real*8 n_particles_w_2
!the local grid
      integer aux_nomega
      real*8  aux_omegamax
      real*8 , dimension (aux_nomega) :: aux_omega
      real*8 , dimension (aux_nomega) :: aux_womega
      real*8,  dimension (-aux_nomega:aux_nomega) ::  spectrum

!counter
      integer :: i_basis
      integer :: j_basis
      integer :: i_freq
      integer :: i_par
      integer :: i_state
      integer :: sign_freq
      integer :: n_step
      integer :: i_dat
      integer i_index
      real*8 e_fermi
      integer n_homo
 
!define a grid on which the spectrum is calculated
      e_fermi = chemical_potential
      if(.not. allocated(aux_green)) allocate(aux_green(n_matrix, n_matrix))

      do i_state = 1, n_states
        if(occ_numbers(i_state,1,1).gt.1.d-6) then
          n_homo=i_state
        endif
      enddo 


      spectrum (:) = 0.d0

      if (anacon_type.eq.0)then
       do sign_freq = -1, 1, 2
        do i_freq = 1, aux_nomega, 1
         do i_basis=1, n_matrix, 1
!         do j_basis=1, n_matrix, 1

         j_basis = i_basis         
         numerator = (0.d0,0.d0)
         denominator = (0.d0,0.d0)
         w = sign_freq*aux_omega(i_freq) 

         if(mod(n_max_par,2) .eq.0) then

            do i_par = n_max_par/2, 2, -1 
              numerator = (numerator + &
               green_fn_par (i_par,i_basis,j_basis))&
              *(w)
            enddo

            do i_par = n_max_par, n_max_par/2+1, -1 
              denominator = (denominator + &
               green_fn_par (i_par,i_basis, j_basis))&
              *(w)
            enddo

         else

            do i_par = n_max_par/2+1, 2, -1 
              numerator = (numerator + &
               green_fn_par (i_par,i_basis,j_basis))&
              *(w)
            enddo

            do i_par = n_max_par, n_max_par/2+2, -1 
              denominator = (denominator + &
               green_fn_par (i_par,i_basis,j_basis))&
              *(w)
            enddo

         endif
         
         numerator = numerator +  green_fn_par (1,i_basis,j_basis)
         denominator = denominator + dcmplx(1.d0,0.d0)
         aux_green (i_basis, j_basis ) = numerator/denominator

!        enddo !i_basis
        enddo !j_basis

        i_index = 0 
       do i_basis = 1, n_matrix , 1
!        do j_basis = 1, n_matrix , 1
          j_basis = i_basis

          i_index = i_index+1

          spectrum (sign_freq*i_freq) = spectrum (sign_freq*i_freq) -&
             (aimag (aux_green (j_basis,i_basis)))*&
             !2./pi
             1./pi
 
        enddo



       enddo  ! i_freq
      enddo !sign_freq

             endif

! the integral of the spectrum should give the number of particle!
! let's check
      if(.false.)then
         n_particles = 0.d0
!       do sign_freq = -1, 1, 2
         do i_freq= 1, aux_nomega, 1
          if( aux_omega(i_freq) .lt. chemical_potential )then
            n_particles = n_particles + spectrum(i_freq)*aux_womega(i_freq)
          else
            n_particles = n_particles +&! 
              spectrum(-i_freq)*aux_womega(i_freq)
          endif

         enddo
!        enddo


         n_particles_w = 0.d0
         n_particles_w_1 = 0.d0
         n_particles_w_2 = 0.d0
       do sign_freq = -1, 1, 2
         do i_freq= 1, aux_nomega, 1
          if( sign_freq*aux_omega(i_freq)*hartree .le. -130 )then
            n_particles_w = n_particles_w + &
            spectrum(sign_freq*i_freq)*aux_womega(i_freq)
          endif 
         enddo
         do i_freq= 1, aux_nomega, 1

          if ( sign_freq*aux_omega(i_freq)*hartree .le. -80 ) then
            n_particles_w_1 = n_particles_w_1 + &
            spectrum(sign_freq*i_freq)*aux_womega(i_freq)
          endif
         enddo
         do i_freq= 1, aux_nomega, 1
          if ( sign_freq*aux_omega(i_freq)*hartree .le. 0) then !chemical_potential*hartree) then
            n_particles_w_2 = n_particles_w_2 + spectrum(sign_freq*i_freq) &
                              *aux_womega(i_freq)

          endif

         enddo
        enddo
       
      endif
!!stop
      return

      end subroutine get_spectrum_k_dmft
